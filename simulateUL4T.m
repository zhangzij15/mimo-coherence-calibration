function result = simulateUL4T(cfg)
%SIMULATEUL4T Paired 2T or 4T PUSCH comparison on a common receiver noise sweep.
if nargin<1, cfg = configUL4T(); end
cfg = prepareFixedArrayConfig(cfg);
ue = ul4tUEArray(cfg);
nTx = numel(ue.antenna_indices);
cfg.ul.num_tx = nTx;
cfg = validateULConfig(cfg);
ul = cfg.ul;
identityTx = eye(nTx);
nCandidate = size(ue.candidate_positions_m,1);
schedule = ul4tSchedule(cfg);
driftEnabled = ul.coherence_enabled && strcmp(ul.coherence_model,'drift');
normalizeTxPower = ul.hold_tx_power && ~ul.normalize_pusch_energy;
useToolbox = strcmp(cfg.cdl.backend,'5g-toolbox');
if useToolbox
    required = {'nrCDLChannel','nrCarrierConfig','nrPUSCHConfig', ...
        'nrPUSCHDMRS','nrPUSCHDMRSIndices'};
    for j = 1:numel(required)
        assert(~isempty(which(required{j})), ...
            'The optional 5G Toolbox backend requires %s. Use base-cdl for Base MATLAB.',required{j});
    end
end
pattern = fixedArrayChannel('pattern',cfg);

% Normal CP and a full 14-symbol NR PUSCH allocation. Two front-loaded DMRS.
N = cfg.N;
F = cfg.ofdm.num_tones;
r = ul.rank;
[pilot,groups] = buildDMRS(cfg,0,0);
dmrsTone = find(any(any(abs(pilot)>0,3),2));
linkTable = ul4tLinkMapping('prepare',cfg,groups);
[pilotBySlot,rotationBySlot] = slotPilots(cfg,schedule,dmrsTone,0);
firstSlot = 1+mod(schedule.ul_slot(1),schedule.slots_per_frame);
pilotDesign = reshape(pilotBySlot{firstSlot},[],r);
pilotDesign = bsxfun(@times,pilotDesign,conj(rotationBySlot{firstSlot}));

nc = numel(ul.cases);
calibrationMode = {ul.cases.calibration};
calibrated = ~strcmp(calibrationMode,'none');
feedback = strcmp(calibrationMode,'feedback');
oracle = strcmp(calibrationMode,'oracle');
hasCalibration = any(calibrated);
if hasCalibration
    calPilot = ul4tCalibrationPilots(cfg);
    [calDMRS,calRotation] = slotPilots(cfg,schedule,dmrsTone,calPilot.prefix_symbols);
    calLinkTable = ul4tLinkMapping('prepare',cfg,groups,calPilot.prefix_symbols);
    calRows = bsxfun(@plus,(1:N).',N*(calPilot.tone-1));
    calIndices = bsxfun(@plus,calRows,N*F*(calPilot.port-1));
    caseIndices = find(calibrated);
end
book = cell(1,nc);
tpmi = cell(1,nc);
labels = cell(1,nc);
if nTx==2
    nrCodebook = @nrUL2TCodebook;
else
    nrCodebook = @nrUL4TCodebook;
end
for c = 1:nc
    mode = ul.cases(c);
    labels{c} = upper(mode.precoder);
    if strcmp(mode.precoder,'nr')
        [book{c},tpmi{c}] = nrCodebook(r,mode.sets,ul.codebook_unit_power);
        category = mode.sets;
        if ischar(category), category = {category}; end
        labels{c} = ['NR (' strjoin(category,'+') ')'];
    end
    if feedback(c)
        labels{c} = [labels{c} ' / RF correction'];
    elseif oracle(c)
        labels{c} = [labels{c} ' / oracle RF correction'];
    elseif calibrated(c)
        labels{c} = [labels{c} ' / precoder reselection'];
    elseif driftEnabled && mode.coherence_error
        labels{c} = [labels{c} ' / RF drift'];
    elseif ul.coherence_enabled && mode.coherence_error
        labels{c} = [labels{c} ' / coherence error'];
    else
        labels{c} = [labels{c} ' / ideal coherence'];
    end
end
noiseVar = cfg.ul.rx_noise_power;
ns = numel(noiseVar);
nu = numel(schedule.ul_ms);
nSrs = numel(schedule.srs_ms);
dims = [cfg.MC nu ns nc];
result.good_bits = zeros(dims);
result.mcs = -ones(dims);
result.eff_sinr_dB = zeros(dims);
% A common EESM beta supports paired gains without MCS-dependent beta jumps.
result.fixed_beta_eff_sinr_dB = zeros(dims);
result.fixed_beta_linear = 1;
result.dmrs_nmse = zeros(dims);
result.dmrs_filter_cache = struct('enabled',ul.dmrs_filter_cache_enabled, ...
    'hits',zeros(cfg.MC,1),'misses',zeros(cfg.MC,1));
result.precoder_tpmi = -ones(cfg.MC,nSrs,nc);
result.srs_propagation_tx_correlation = zeros(cfg.MC,nSrs);
result.precoder = complex(zeros(nTx,r,cfg.MC,nSrs,nc));
result.srs_initial_phase_rad = zeros(nTx,cfg.MC,nSrs,nc);
result.coherence_gain = complex(ones(nTx,cfg.MC,nu));
result.pusch_tx_power = zeros(cfg.MC,nu,nc);
result.pusch_srs_reference_energy = zeros(cfg.MC,nu,nc);
result.pusch_reference_energy = zeros(cfg.MC,nu,nc);
result.pusch_effective_energy = zeros(cfg.MC,nu,nc);
result.ue_position_m = zeros(cfg.MC,3);
if hasCalibration
    result.calibration.pilot = calPilot;
    result.calibration.case_indices = caseIndices;
    result.calibration.data_re_per_rb = calLinkTable.data_re_per_rb;
    matchedBits = zeros(dims);
    calDims = [nTx cfg.MC nu ns numel(caseIndices)];
    result.calibration.estimated_gain = complex(ones(calDims));
    result.calibration.innovation = complex(ones(calDims));
    result.calibration.tx_correction = complex(ones(calDims));
    result.calibration.tpmi = -ones(cfg.MC,nu,ns,numel(caseIndices));
    result.calibration.precoder = complex(zeros(nTx,r,cfg.MC,nu,ns,numel(caseIndices)));
    result.calibration.tx_power = zeros(cfg.MC,nu,ns,numel(caseIndices));
    result.calibration.rx_energy = zeros(cfg.MC,nu,ns,numel(caseIndices));
    result.calibration.estimation_phase_error_deg = zeros(calDims);
    result.calibration.estimation_amplitude_error_dB = zeros(calDims);
    result.calibration.estimation_error_definition = ...
        'True/estimated SRS-relative RF gain after each is divided by its TX1 gain; not RAN4 per-port drift.';
    calPowerDims = [cfg.MC nu ns numel(caseIndices)];
    result.calibration.pilot_epre_min = NaN(calPowerDims);
    result.calibration.pilot_epre_max = NaN(calPowerDims);
    result.calibration.pilot_total_energy = NaN(calPowerDims);
    result.calibration.pilot_power_definition = ...
        'Digital probe RE power after correction and normalization, before unknown RF gain; Oracle does not simulate a calibration probe.';
end
if driftEnabled
    % Sample every slot, including downlink and idle slots, on one RF trajectory.
    time_ms = (0:round(ul.duration_ms/schedule.slot_ms)-1)*schedule.slot_ms;
    srsIndex = 1+round(schedule.srs_ms/schedule.slot_ms);
    puschIndex = 1+schedule.ul_slot;
    referenceIndex = srsIndex(schedule.srs_id);
    result.coherence_drift.time_ms = time_ms;
    result.coherence_drift.srs_sample_index = srsIndex;
    result.coherence_drift.pusch_sample_index = puschIndex;
    result.coherence_drift.phase_deg = zeros(nTx,cfg.MC,numel(time_ms));
    result.coherence_drift.amplitude_dB = zeros(nTx,cfg.MC,numel(time_ms));
end
previous = rng;
restoreRng = onCleanup(@() rng(previous));
mcOffset = 0:cfg.MC-1;
geometrySeed = mod(double(cfg.seed)+mcOffset,2^32);
initialPhaseSeed = mod(double(cfg.seed)+2^28+mcOffset,2^32);
mainNoiseSeed = mod(double(cfg.seed)+2^29+mcOffset,2^32);
result.randomization = struct('version','mc_streams_v1', ...
    'geometry_seed',geometrySeed,'initial_phase_seed',initialPhaseSeed, ...
    'main_noise_seed',mainNoiseSeed);
for mc = 1:cfg.MC
    % Earlier MC rank/layout/case draw counts must not change this scene or motion.
    rng(geometrySeed(mc),'twister');
    if useToolbox
        [model,scene] = ul4tToolboxChannel(cfg,pattern,schedule,mc,ue);
    else
        [model,scene] = ul4tNativeCDL(cfg,pattern,schedule,mc,ue);
    end
    result.channel_diagnostics(mc) = model.diagnostics;
    result.ue_position_m(mc,:) = scene.positions_m;
    agingEnabled = ul.channel_aging_enabled && model.max_doppler_Hz>0;
    basis = model.delayBasis;
    B = size(basis,2);
    pilotBasis = repmat(basis(dmrsTone,:),2,1);
    X = complex(zeros(size(pilotDesign,1),B*r));
    for layer = 1:r
        cols = (layer-1)*B+(1:B);
        X(:,cols) = bsxfun(@times,pilotBasis,pilotDesign(:,layer));
    end
    % Priors use the configured CDL PDP and the SRS-selected precoder.
    % No noiseless PUSCH coefficients enter this Wiener filter.
    XHX = X'*X;
    % pRx/X/XHX are fixed only within this MC. Retain one K per SNR.
    dmrsCache = struct('power',{cell(1,ns)},'filter',{cell(1,ns)}, ...
        'hits',0,'misses',0);
    if hasCalibration
        mainRng = rng;
        rng(mod(double(cfg.seed)+2^30+mc-1,2^32),'twister');
        calibrationNoise = (randn(N,calPilot.num_re,nu)+ ...
            1i*randn(N,calPilot.num_re,nu))/sqrt(2);
        rng(mainRng);
    end
    % Keep the main RNG sequence paired with the legacy uniform simulation.
    phaseError = (2*rand(nCandidate,nu)-1)*ul.coherence_deg;
    ampError = (2*rand(nCandidate,nu)-1)*ul.coherence_dB;
    phaseError = phaseError(ue.antenna_indices,:);
    ampError = ampError(ue.antenna_indices,:);
    nominalPhaseDeg = [];
    if driftEnabled
        mainRng = rng;
        rng(mod(double(cfg.seed)+2^31+mc-1,2^32),'twister');
        if isfield(ul.coherence_drift,'bounded')
            bounded = ul.coherence_drift.bounded;
            if strcmp(bounded.mode,'load')
                % State at a slot boundary depends on previous TX activity.
                bounded.txActive = false(size(time_ms));
                bounded.txActive([srsIndex puschIndex]) = true;
            end
            if isfield(bounded,'initial_offset_deg') || isfield(bounded,'initial_jitter_deg')
                [phaseDeg,ampDb,nominalPhaseDeg] = ul4tBoundedRfDrift( ...
                    nCandidate,time_ms,ue.antenna_indices,bounded);
            else
                [phaseDeg,ampDb] = ul4tBoundedRfDrift( ...
                    nCandidate,time_ms,ue.antenna_indices,bounded);
            end
        else
            [phaseDeg,ampDb] = ul4tCoherenceDrift(nCandidate,time_ms,ul.coherence_drift);
        end
        rng(mainRng);
        phaseDeg = phaseDeg(ue.antenna_indices,:);
        ampDb = ampDb(ue.antenna_indices,:);
        result.coherence_drift.phase_deg(:,mc,:) = reshape(phaseDeg,nTx,1,[]);
        result.coherence_drift.amplitude_dB(:,mc,:) = reshape(ampDb,nTx,1,[]);
        phaseError = phaseDeg(:,puschIndex)-phaseDeg(:,referenceIndex);
        ampError = ampDb(:,puschIndex)-ampDb(:,referenceIndex);
    end
    gains = 10.^(ampError/20).*exp(1i*pi/180*phaseError);
    if ~ul.coherence_enabled, gains = ones(nTx,nu); end
    result.coherence_gain(:,mc,:) = reshape(gains,nTx,1,nu);
    % One hardware phase per physical location and MC, shared by every case.
    rng(initialPhaseSeed(mc),'twister');
    initialPhase = 2*pi*rand(nCandidate,1);
    if ~isempty(nominalPhaseDeg)
        % Use the configured phase at the SRS reference. For load mode this
        % is the cold phase when the first SRS precedes any TX activity.
        initialPhase = (pi/180)*nominalPhaseDeg;
    end
    initialPhase = repmat(initialPhase(ue.antenna_indices),1,nc);
    initialGain = exp(1i*initialPhase);
    result.srs_initial_phase_rad(:,mc,:,:) = repmat(reshape(initialPhase,nTx,1,1,nc),1,1,nSrs,1);
    Hsrs = cell(nSrs,1);
    WwithInitialPhase = cell(nSrs,nc);
    for s = 1:nSrs
        Hsrs{s} = reshape(ul4tResponse(model,schedule.srs_ms(s),identityTx),[],nTx);
        propagationGram = (Hsrs{s}'*Hsrs{s})/F;
        columnNorms = sqrt(max(real(diag(propagationGram)),realmin));
        normalizedGram = abs(propagationGram)./(columnNorms*columnNorms.');
        result.srs_propagation_tx_correlation(mc,s) = mean(normalizedGram(~eye(nTx)));
        for c = 1:nc
            % Ideal sounding includes the fixed initial phase for this Case.
            portPhase = initialGain(:,c);
            gram = propagationGram.*(conj(portPhase)*portPhase.');
            [selectedW,id] = selectPrecoder(gram,r,book{c},tpmi{c},ul.selection_noise_power);
            result.precoder_tpmi(mc,s,c) = id;
            result.precoder(:,:,mc,s,c) = selectedW;
            WwithInitialPhase{s,c} = bsxfun(@times,portPhase,selectedW);
        end
    end
    filter = cell(ns,nc);
    srsEffective = cell(1,nc);
    srsNorm = zeros(1,nc);
    previousSrs = 0;
    % PUSCH noise starts independently of geometry, RF and hardware phase draws.
    rng(mainNoiseSeed(mc),'twister');
    for u = 1:nu
        s = schedule.srs_id(u);
        if s~=previousSrs
            if hasCalibration
                trackedGain = complex(ones(nTx,ns,numel(caseIndices)));
                correction = complex(ones(nTx,ns,numel(caseIndices)));
                srsPilot = Hsrs{s}(calIndices);
                srsGram = (Hsrs{s}'*Hsrs{s})/F;
                if agingEnabled
                    pathPilot = ul4tCalibrationPathBasis(model,Hsrs{s},calIndices,calPilot.port,N,F);
                end
            end
            for c = 1:nc
                nominalW = WwithInitialPhase{s,c};
                srsEffective{c} = reshape(Hsrs{s}*nominalW,N,F,r);
                srsNorm(c) = norm(srsEffective{c}(:));
                if ~calibrated(c)
                    for i = 1:ns
                        [filter{i,c},dmrsCache] = cachedDmrsFilter(dmrsCache,i, ...
                            model,nominalW,X,XHX,noiseVar(i),ul.dmrs_filter_cache_enabled);
                    end
                end
            end
            previousSrs = s;
        end
        % Reuse the SRS precoder and evaluate propagation at the PUSCH time.
        useCurrentChannel = agingEnabled && schedule.srs_age_ms(u)>0;
        Hpusch = Hsrs{s};
        if useCurrentChannel
            Hpusch = reshape(ul4tResponse(model,schedule.ul_ms(u),identityTx),[],nTx);
        end
        g = gains(:,u);
        slot = 1+mod(schedule.ul_slot(u),schedule.slots_per_frame);
        P = reshape(pilotBySlot{slot},1,numel(dmrsTone),2,r);
        rot = rotationBySlot{slot};
        nPilot = numel(rot);
        z = (randn(N,nPilot)+1i*randn(N,nPilot))/sqrt(2);
        if hasCalibration
            puschPilot = Hpusch(calIndices);
            calP = reshape(calDMRS{slot},1,numel(dmrsTone),2,r);
            calRot = calRotation{slot};
        end
        for c = 1:nc
            nominalW = WwithInitialPhase{s,c};
            txPower = norm(nominalW,'fro')^2;
            G = srsEffective{c};
            referenceNorm = srsNorm(c);
            if useCurrentChannel
                G = reshape(Hpusch*nominalW,N,F,r);
                referenceNorm = norm(G(:));
            end
            result.pusch_srs_reference_energy(mc,u,c) = srsNorm(c)^2;
            result.pusch_reference_energy(mc,u,c) = referenceNorm^2;
            if calibrated(c)
                k = find(caseIndices==c,1);
                for i = 1:ns
                    gainBefore = trackedGain(:,i,k);
                    correctionBefore = correction(:,i,k);
                    w = result.precoder(:,:,mc,s,c);
                    id = result.precoder_tpmi(mc,s,c);
                    if oracle(c)
                        % Oracle: exact relative RF gain at this PUSCH, referenced
                        % to the latest SRS. A common complex scalar is irrelevant.
                        gainAfter = g./g(1);
                        innovation = gainAfter./gainBefore;
                        trackedGain(:,i,k) = gainAfter;
                        correction(:,i,k) = 1./gainAfter;
                    else
                        % The previous effective channel and probe weights predict this RS.
                        probe = correctionBefore(calPilot.port);
                        probe = calPilot.amplitude*probe/norm(probe)*sqrt(calPilot.num_re);
                        probePower = abs(probe).^2;
                        result.calibration.pilot_epre_min(mc,u,i,k) = min(probePower);
                        result.calibration.pilot_epre_max(mc,u,i,k) = max(probePower);
                        result.calibration.pilot_total_energy(mc,u,i,k) = sum(probePower);
                        pilotWeight = reshape(initialGain(calPilot.port,c),1,[]).*reshape(probe,1,[]);
                        referenceGain = reshape(gainBefore(calPilot.port),1,[]);
                        reference = bsxfun(@times,srsPilot,pilotWeight.*referenceGain);
                        received = bsxfun(@times,puschPilot,pilotWeight.*reshape(g(calPilot.port),1,[])) ...
                            +sqrt(noiseVar(i))*calibrationNoise(:,:,u);
                        if agingEnabled
                            gainAfter = ul4tTrackCalibration(gainBefore,received,pathPilot,pilotWeight, ...
                                calPilot.port,noiseVar(i),ul.calibration.innovation_variance);
                            innovation = gainAfter./gainBefore;
                        else
                            innovation = ul4tEstimateCalibration(reference,received,calPilot.port, ...
                                nTx,noiseVar(i),ul.calibration.innovation_variance);
                            gainAfter = gainBefore.*innovation;
                        end
                        trackedGain(:,i,k) = gainAfter;
                        if feedback(c)
                            correction(:,i,k) = correctionBefore./innovation;
                        else
                            port = initialGain(:,c).*gainAfter;
                            gram = srsGram.*(conj(port)*port.');
                            [w,id] = selectPrecoder(gram,r,book{c},tpmi{c}, ...
                                ul.selection_noise_power,abs(gainAfter));
                        end
                    end
                    txCorrection = correction(:,i,k);
                    wWithCorrection = bsxfun(@times,initialGain(:,c).*txCorrection,w);
                    targetNorm = norm(w,'fro');
                    actualW = normalizePrecoder(bsxfun(@times,g,wWithCorrection),targetNorm);
                    predictedW = normalizePrecoder(bsxfun(@times,gainAfter,wWithCorrection),targetNorm);
                    G = reshape(Hpusch*actualW,N,F,r);
                    [K,dmrsCache] = cachedDmrsFilter(dmrsCache,i, ...
                        model,predictedW,X,XHX,noiseVar(i),ul.dmrs_filter_cache_enabled);
                    Y = dmrsObservation(G,dmrsTone,calP,calRot);
                    % Main throughput uses the common full PUSCH data budget.
                    [bits,mcs,eff,nmse,sinr,fixedEff] = evaluatePusch(Y,z,K,basis,G,noiseVar(i),linkTable);
                    result.good_bits(mc,u,i,c) = bits;
                    matchedBits(mc,u,i,c) = ul4tLinkMapping('evaluate',sinr,calLinkTable);
                    result.mcs(mc,u,i,c) = mcs;
                    result.eff_sinr_dB(mc,u,i,c) = eff;
                    result.fixed_beta_eff_sinr_dB(mc,u,i,c) = fixedEff;
                    result.dmrs_nmse(mc,u,i,c) = nmse;
                    result.calibration.estimated_gain(:,mc,u,i,k) = gainAfter;
                    result.calibration.innovation(:,mc,u,i,k) = innovation;
                    result.calibration.tx_correction(:,mc,u,i,k) = txCorrection;
                    result.calibration.precoder(:,:,mc,u,i,k) = w;
                    result.calibration.tpmi(mc,u,i,k) = id;
                    result.calibration.tx_power(mc,u,i,k) = norm(actualW,'fro')^2;
                    result.calibration.rx_energy(mc,u,i,k) = norm(G(:))^2;
                    relativeTruth = g/g(1);
                    relativeEstimate = gainAfter/gainAfter(1);
                    estimationResidual = relativeTruth./relativeEstimate;
                    result.calibration.estimation_phase_error_deg(:,mc,u,i,k) = ...
                        angle(estimationResidual)*180/pi;
                    result.calibration.estimation_amplitude_error_dB(:,mc,u,i,k) = ...
                        20*log10(abs(estimationResidual));
                end
                result.pusch_tx_power(mc,u,c) = result.calibration.tx_power(mc,u,1,k);
                result.pusch_effective_energy(mc,u,c) = result.calibration.rx_energy(mc,u,1,k);
                continue
            end
            effectiveNorm = referenceNorm;
            if ul.coherence_enabled && ul.cases(c).coherence_error
                actualW = bsxfun(@times,g,nominalW);
                if normalizeTxPower
                    actualW = actualW*(sqrt(txPower)/norm(actualW,'fro'));
                end
                txPower = norm(actualW,'fro')^2;
                G = reshape(Hpusch*actualW,N,F,r);
                effectiveNorm = norm(G(:));
                if ul.normalize_pusch_energy && effectiveNorm>0
                    % One scale for the full effective channel, shared by pilots and data.
                    G = G.*(referenceNorm/effectiveNorm);
                    effectiveNorm = norm(G(:));
                end
            end
            result.pusch_tx_power(mc,u,c) = txPower;
            result.pusch_effective_energy(mc,u,c) = effectiveNorm^2;
            Y = dmrsObservation(G,dmrsTone,P,rot);
            for i = 1:ns
                [bits,mcs,eff,nmse,sinr,fixedEff] = evaluatePusch(Y,z,filter{i,c},basis,G,noiseVar(i),linkTable);
                if hasCalibration
                    matchedBits(mc,u,i,c) = ul4tLinkMapping('evaluate',sinr,calLinkTable);
                end
                result.dmrs_nmse(mc,u,i,c) = nmse;
                result.good_bits(mc,u,i,c) = bits;
                result.mcs(mc,u,i,c) = mcs;
                result.eff_sinr_dB(mc,u,i,c) = eff;
                result.fixed_beta_eff_sinr_dB(mc,u,i,c) = fixedEff;
            end
        end
    end
    result.dmrs_filter_cache.hits(mc) = dmrsCache.hits;
    result.dmrs_filter_cache.misses(mc) = dmrsCache.misses;
    if mod(mc,cfg.progress_every)==0 || mc==cfg.MC
        fprintf('%dT uplink: %d/%d MC complete.\n',nTx,mc,cfg.MC);
    end
end
result.cfg = cfg;
result.ue_array = ue;
result.schedule = schedule;
result.rx_noise_power = noiseVar;
result.inverse_noise_power = 1./noiseVar;
result.inverse_noise_power_dB = 10*log10(result.inverse_noise_power);
result.max_doppler_Hz = model.max_doppler_Hz;
result.labels = labels;
result.sem = linkTable; % historical MAT field for the link mapping table
result.receiver_definition = struct('method','MMSE-IRC', ...
    'covariance','White receive noise (sigma^2 I)', ...
    'external_interference',false);
dmrsSegments = cfg.ul.dmrs_segments_rb;
if isempty(dmrsSegments), dmrsSegments = cfg.ofdm.num_rb; end
result.resource_grid_definition = struct('num_rb',cfg.ofdm.num_rb, ...
    'num_tones',cfg.ofdm.num_tones, ...
    'active_bandwidth_Hz',cfg.ofdm.active_bandwidth_Hz, ...
    'dmrs_segments_rb',dmrsSegments, ...
    'is_extended_research_grid',cfg.ofdm.num_rb>275);
result.throughput_definition = struct( ...
    'main','Common PUSCH data RE budget; calibration prefix overhead is not deducted.', ...
    'main_data_re_per_rb',linkTable.data_re_per_rb, ...
    'frame_scope','Includes the actual scheduled uplink-slot duty factor.', ...
    'link_abstraction','EESM/SEM and AWGN BLER-0.1 thresholds; expected good bits are 0.9 times the selected TBS.');
result.thp_uplink_Mbps = reshape(mean(mean(result.good_bits,1),2),ns,nc)/(schedule.slot_ms*1000);
result.thp_frame_Mbps = schedule.ul_duty*result.thp_uplink_Mbps;
if hasCalibration
    result.throughput_definition.matched_diagnostic = ...
        'All cases use the reduced data RE budget after the calibration prefix; auxiliary equal-data-RE diagnostic.';
    result.throughput_definition.matched_data_re_per_rb = calLinkTable.data_re_per_rb;
    result.calibration.matched_thp_uplink_Mbps = ...
        reshape(mean(mean(matchedBits,1),2),ns,nc)/(schedule.slot_ms*1000);
    result.calibration.matched_thp_frame_Mbps = ...
        schedule.ul_duty*result.calibration.matched_thp_uplink_Mbps;
end
end

function [pilot,groups] = buildDMRS(cfg,slotInFrame,startSymbol)
if strcmp(cfg.cdl.backend,'5g-toolbox')
    [pilot,groups] = ul4tToolboxDMRS(cfg,slotInFrame,startSymbol);
else
    [pilot,groups] = ul4tBaseDMRS(cfg,slotInFrame,startSymbol);
end
end

function [pilots,rotations] = slotPilots(cfg,schedule,dmrsTone,startSymbol)
pilots = cell(1,schedule.slots_per_frame);
rotations = cell(size(pilots));
for n = unique(mod(schedule.ul_slot,schedule.slots_per_frame))
    pilot = buildDMRS(cfg,n,startSymbol);
    pilot = pilot(dmrsTone,:,:);
    common = pilot(:,:,1);
    for layer = 2:cfg.ul.rank
        missing = (abs(common)==0);
        next = pilot(:,:,layer);
        common(missing) = next(missing);
    end
    common = common./abs(common);
    pilots{n+1} = pilot;
    rotations{n+1} = common(:);
end
end

function [w,id] = selectPrecoder(gram,rank,book,tpmi,selectionNoise,rfAmplitude)
id = -1;
if isempty(book)
    if nargin>=6
        % Select equal-power eigenmodes in the physical transmit coordinates.
        inverse = 1./rfAmplitude(:);
        [~,~,V] = svd(gram.*(inverse*inverse.'));
        w = normalizePrecoder(bsxfun(@times,inverse,V(:,1:rank)),1);
        return
    end
    [~,~,V] = svd(gram);
    w = V(:,1:rank)/sqrt(rank);
    return
end
score = zeros(1,numel(tpmi));
for b = 1:numel(score)
    candidate = book(:,:,b);
    projected = candidate'*gram*candidate;
    if nargin>=6
        actual = bsxfun(@times,rfAmplitude,candidate);
        projected = projected*(norm(candidate,'fro')/norm(actual,'fro'))^2;
    end
    if rank==1
        score(b) = real(projected);
    else
        score(b) = sum(log1p(max(0,real(eig(projected)))/selectionNoise));
    end
end
[~,b] = max(score);
w = book(:,:,b);
id = tpmi(b);
end

function w = normalizePrecoder(w,targetNorm)
w = w*(targetNorm/norm(w,'fro'));
end

function [K,cache] = cachedDmrsFilter(cache,index,model,w,X,XHX,noiseVar,enabled)
% Toolbox's prior depends on layer powers, not the precoder's port phases.
% Compare within roundoff only; every miss retains the original prior/solve.
useCache = enabled && isfield(model,'pathGains');
if useCache
    power = sum(abs(w).^2,1);
    previousPower = cache.power{index};
    if ~isempty(previousPower)
        tolerance = 64*eps(max(1,max(abs([power previousPower]))));
        if all(abs(power-previousPower)<=tolerance)
            K = cache.filter{index};
            cache.hits = cache.hits+1;
            return
        end
    end
end
K = dmrsFilter(model,w,X,XHX,noiseVar);
cache.misses = cache.misses+1;
if useCache
    cache.power{index} = power;
    cache.filter{index} = K;
end
end

function K = dmrsFilter(model,w,X,XHX,noiseVar)
if isfield(model,'pathGains')
    spatialPower = sum(abs(w).^2,1);
else
    spatialPower = abs(model.txPhase*w).^2;
end
prior = bsxfun(@times,spatialPower,model.pRx);
rootPrior = sqrt(prior(:));
XP = bsxfun(@times,X,rootPrior.');
A = XHX.*(rootPrior*rootPrior.');
K = (bsxfun(@times,rootPrior,(A+noiseVar*eye(size(A)))\XP'))';
end

function Y = dmrsObservation(G,dmrsTone,pilot,rotation)
N = size(G,1);
r = size(G,3);
Y = sum(bsxfun(@times,reshape(G(:,dmrsTone,:),N,numel(dmrsTone),1,r),pilot),4);
Y = bsxfun(@times,reshape(Y,N,[]),conj(rotation.'));
end

function [bits,mcs,eff,nmse,sinr,fixedEff] = evaluatePusch(Y,z,K,basis,G,noiseVar,table)
[N,F,r] = size(G);
B = size(basis,2);
coefficients = (Y+sqrt(noiseVar)*z)*conj(K);
Ghat = complex(zeros(N,F,r));
for layer = 1:r
    cols = (layer-1)*B+(1:B);
    Ghat(:,:,layer) = coefficients(:,cols)*basis.';
end
nmse = norm(Ghat(:)-G(:))^2/max(realmin,norm(G(:))^2);
sinr = postSinr(Ghat,G,noiseVar);
[bits,mcs,eff] = ul4tLinkMapping('evaluate',sinr,table);
% Stable EESM with fixed linear beta=1, independent of selected MCS/TBS.
s = max(0,sinr(:));
floorSINR = min(s);
fixedLinear = floorSINR-log(mean(exp(-(s-floorSINR))));
fixedEff = 10*log10(max(fixedLinear,realmin));
end

function sinr = postSinr(Hhat,H,noiseVar)
sinr = ul4tPostSinr(Hhat,H,noiseVar);
end

function cfg = validateULConfig(cfg)
if ~isfield(cfg.ul,'normalize_pusch_energy'), cfg.ul.normalize_pusch_energy = false; end
if ~isfield(cfg.ul,'hold_tx_power'), cfg.ul.hold_tx_power = true; end
if ~isfield(cfg.ul,'channel_aging_enabled'), cfg.ul.channel_aging_enabled = true; end
if ~isfield(cfg.ul,'dmrs_filter_cache_enabled'), cfg.ul.dmrs_filter_cache_enabled = true; end
if isfield(cfg.ul,'srs_timing'), cfg.ul = rmfield(cfg.ul,'srs_timing'); end
if ~isfield(cfg.ul,'coherence_model'), cfg.ul.coherence_model = 'uniform'; end
if strcmp(cfg.ul.coherence_model,'drift') && ~isfield(cfg.ul,'coherence_drift')
    defaults = configUL4T();
    cfg.ul.coherence_drift = defaults.ul.coherence_drift;
end
if ~isfield(cfg.ul,'rx_noise_power')
    defaults = configUL4T();
    cfg.ul.rx_noise_power = defaults.ul.rx_noise_power;
end
if ~isfield(cfg.ul,'selection_noise_power'), cfg.ul.selection_noise_power = 1; end
if ~isfield(cfg.ul,'dmrs_segments_rb'), cfg.ul.dmrs_segments_rb = []; end
if ~isfield(cfg.cdl,'backend') || strcmp(cfg.cdl.backend,'base-cdl-c')
    cfg.cdl.backend = 'base-cdl';
end
ul = cfg.ul;
assert(cfg.K==1,'This entry point requires one UE.');
assert(cfg.ofdm.enabled && strcmp(cfg.ofdm.resource_mode,'rb') && ...
    strcmp(cfg.ofdm.frequency_mode,'full') && cfg.ofdm.subcarriers_per_rb==12 && ...
    isempty(cfg.ofdm.rb_allocation), ...
    'Use an ungrouped allocation of contiguous 12-tone RBs.');
assert(isempty(cfg.irc.external_covariance),'The SU receiver currently assumes white noise.');
assert(isscalar(ul.rank) && ismember(ul.rank,1:min(ul.num_tx,cfg.N)), ...
    'Set rank between 1 and the smaller of the UE and BS port counts.');
assert(ismember(ul.nid,0:65535) && ismember(ul.nscid,0:1) && ...
    isnumeric(ul.prb_start_crb) && isreal(ul.prb_start_crb) && isscalar(ul.prb_start_crb) && ...
    isfinite(ul.prb_start_crb) && ul.prb_start_crb>=0 && ...
    ul.prb_start_crb==fix(ul.prb_start_crb),'Invalid DMRS scrambling or starting CRB.');
assert(any(strcmp(cfg.cdl.backend,{'base-cdl','5g-toolbox'})), ...
    'Choose base-cdl or the optional 5g-toolbox backend.');
if strcmp(cfg.cdl.backend,'5g-toolbox')
    segments = ul.dmrs_segments_rb;
    if isempty(segments)
        assert(ul.prb_start_crb+cfg.ofdm.num_rb<=275,'UL4T:DMRSGrid', ...
            'A Toolbox grid extending beyond CRB 274 requires explicit ul.dmrs_segments_rb.');
    else
        assert(isnumeric(segments) && isreal(segments) && isvector(segments) && ...
            all(isfinite(segments)) && all(segments==fix(segments)) && ...
            all(segments>=1 & segments<=275) && sum(segments)==cfg.ofdm.num_rb, ...
            'UL4T:DMRSGrid','Toolbox DM-RS window sizes must be integers from 1 to 275 and sum to num_rb.');
        lastWindowStart = ul.prb_start_crb+cfg.ofdm.num_rb-segments(end);
        assert(lastWindowStart<=2199,'UL4T:DMRSGrid', ...
            'Each Toolbox DM-RS window must start at CRB 2199 or below.');
        cfg.ul.dmrs_segments_rb = double(segments(:).');
    end
else
    cfg.ul.dmrs_segments_rb = [];
end
assert(all(isfinite([ul.ue_speed_kmh ul.coherence_dB ul.coherence_deg])) && ...
    min([ul.ue_speed_kmh ul.coherence_dB ul.coherence_deg])>=0, ...
    'Invalid speed or coherence model.');
assert(strcmp(ul.coherence_model,'uniform') || strcmp(ul.coherence_model,'drift'), ...
    'Set coherence_model to uniform or drift.');
profiles = {'CDL-A','CDL-B','CDL-C','CDL-D','CDL-E'};
assert(any(strcmp(cfg.cdl.toolbox_profile,profiles)) && ...
    isscalar(cfg.cdl.toolbox_delay_spread_s) && isfinite(cfg.cdl.toolbox_delay_spread_s) && ...
    cfg.cdl.toolbox_delay_spread_s>=0,'Use a supported CDL profile and delay spread.');
cfg.ofdm.frequency_dependent_steering = false;
v = ul.rx_noise_power;
assert(isnumeric(v) && isreal(v) && isvector(v) && ~isempty(v) && ...
    all(v(:)>0 & isfinite(v(:)) & isfinite(1./v(:))),'rx_noise_power must contain positive finite noise powers.');
cfg.ul.rx_noise_power = sort(double(v(:).'),'descend');
assert(isscalar(ul.selection_noise_power) && isfinite(ul.selection_noise_power) && ...
    ul.selection_noise_power>0,'selection_noise_power must be positive.');
assert(any(strcmp(ul.thp_scope,{'frame','uplink'})),'Invalid throughput scope.');
assert(isscalar(ul.normalize_pusch_energy) && ismember(ul.normalize_pusch_energy,[0 1]) && ...
    isscalar(ul.hold_tx_power) && ismember(ul.hold_tx_power,[0 1]) && ...
    isscalar(ul.channel_aging_enabled) && ismember(ul.channel_aging_enabled,[0 1]) && ...
    isscalar(ul.dmrs_filter_cache_enabled) && ismember(ul.dmrs_filter_cache_enabled,[0 1]), ...
    'Energy normalization, channel aging and DMRS cache switches must be logical flags.');
assert(~isempty(ul.cases) && all(ismember({ul.cases.precoder},{'svd','nr'})), ...
    'Configure at least one svd or nr case.');
if ~isfield(ul.cases,'calibration')
    [cfg.ul.cases.calibration] = deal('none');
end
modes = {cfg.ul.cases.calibration};
assert(all(ismember(modes,{'none','feedback','reselect','oracle'})), ...
    'Set case calibration to none, feedback, reselect or oracle.');
active = ~strcmp(modes,'none');
if any(active)
    if ~isfield(cfg.ul,'calibration')
        defaults = configUL4T();
        cfg.ul.calibration = defaults.ul.calibration;
    end
    assert(ul.coherence_enabled && all([ul.cases(active).coherence_error]), ...
        'Calibrated cases require the RF error to be enabled.');
    assert(ul.hold_tx_power && ~ul.normalize_pusch_energy && ul.codebook_unit_power, ...
        'Calibration comparisons require fixed, equal TX power and unnormalized receive power.');
end
end


% ===== Consolidated helper: ul4tUEArray.m =====
function ue = ul4tUEArray(cfg)
%UL4TUEARRAY Resolve logical ports to ordered UE antenna positions.
ul = cfg.ul;
nTx = 4; % configurations saved before num_tx use four ports
if isfield(ul,'num_tx'), nTx = ul.num_tx; end
assert(isnumeric(nTx) && isscalar(nTx) && ismember(nTx,[2 4]), ...
    'Set num_tx to 2 or 4.');
nTx = double(nTx);
ue.layout = 'upa';
if isfield(ul,'ue_layout'), ue.layout = ul.ue_layout; end
wavelength = cfg.c0/cfg.fc_Hz;
switch ue.layout
    case 'handheld'
        % TR 38.901 v19.3.0, Figures 7.3-2 and 7.3-3, UT reference LCS.
        halfLength = 0.075;
        halfWidth = 0.035;
        ue.candidate_positions_m = [
            -halfLength -halfWidth 0
                      0 -halfWidth 0
             halfLength -halfWidth 0
             halfLength          0 0
             halfLength  halfWidth 0
                      0  halfWidth 0
            -halfLength  halfWidth 0
            -halfLength          0 0];
        defaultIndices = [1 7 3 5];
        indices = defaultIndices(1:nTx);
        if isfield(ul,'ue_antenna_indices') && ~isempty(ul.ue_antenna_indices)
            indices = ul.ue_antenna_indices;
        end
        assert(isnumeric(indices) && isreal(indices) && isvector(indices) && ...
            numel(indices)==nTx && all(ismember(indices(:),1:8)) && numel(unique(indices))==nTx, ...
            'ue_antenna_indices must be empty or contain %d distinct antenna IDs from 1 to 8.',nTx);
        ue.antenna_indices = double(indices(:).');
        ue.positions_m = ue.candidate_positions_m(ue.antenna_indices,:);
        ue.positions_lambda = ue.positions_m/wavelength;
        % The Toolbox indexes each grid column top to bottom. After -90 deg
        % downtilt, the row coordinate is -z of the original panel.
        gridPorts = [1 2 3 6 9 8 7 4];
        ue.toolbox_ports = gridPorts(ue.antenna_indices);
        ue.toolbox_spacing_lambda = [halfLength halfWidth]/wavelength;
    case 'upa'
        pos = ul.tx_positions_lambda;
        assert(isnumeric(pos) && isreal(pos) && ismatrix(pos) && size(pos,2)==3 && ...
            ismember(size(pos,1),[nTx 4]) && all(isfinite(pos(:))), ...
            'Set num_tx or four real, finite UE port positions in wavelengths.');
        ue.antenna_indices = 1:nTx;
        ue.positions_lambda = pos(ue.antenna_indices,:);
        ue.positions_m = ue.positions_lambda*wavelength;
        ue.candidate_positions_m = pos*wavelength;
        ue.toolbox_ports = ue.antenna_indices;
    otherwise
        error('Use handheld or upa for ue_layout.');
end
end


% ===== Consolidated helper: ul4tSchedule.m =====
function schedule = ul4tSchedule(cfg)
%UL4TSCHEDULE Map TDD blocks and SRS occasions to NR physical slots.
ul = cfg.ul;
mu = log2(cfg.ofdm.scs_Hz/15e3);
assert(ismember(mu,0:3),'SCS must be 15, 30, 60 or 120 kHz.');
slot_ms = 2^(-mu);
assert(all(isfinite([ul.block_ms ul.duration_ms ul.srs_period_ms])) && ...
    all([ul.block_ms ul.duration_ms ul.srs_period_ms]>0),'Invalid TDD/SRS times.');
slotsPerBlock = round(ul.block_ms/slot_ms);
periodSlots = round(ul.srs_period_ms/slot_ms);
numSlots = round(ul.duration_ms/slot_ms);
assert(all(abs([slotsPerBlock periodSlots numSlots]*slot_ms- ...
    [ul.block_ms ul.srs_period_ms ul.duration_ms])<1e-10), ...
    'TDD/SRS times must be multiples of an NR slot.');
pattern = upper(ul.tdd_pattern);
assert(ischar(pattern) && all(ismember(pattern,'DSU')) && ...
    any(pattern=='S') && any(pattern=='U'),'TDD pattern requires S and U.');
slots = 0:numSlots-1;
tag = pattern(1+mod(floor(slots/slotsPerBlock),numel(pattern)));
specialStart = slots(tag=='S' & mod(slots,slotsPerBlock)==0);
assert(~isempty(specialStart),'Duration must include at least one S block.');
% Update on the first S block eligible since the last actual SRS.
srsSlot = specialStart(1);
for t = specialStart(2:end)
    if t-srsSlot(end)>=periodSlots
        srsSlot(end+1) = t; %#ok<AGROW>
    end
end
uplinkSlot = slots(tag=='U' & slots>=srsSlot(1));
assert(~isempty(uplinkSlot),'Duration must contain a U slot after the first SRS.');
% All uplink slots refer to their most recent accepted SRS.
srsId = sum(bsxfun(@ge,uplinkSlot(:),srsSlot),2).';
schedule.slot_ms = slot_ms;
schedule.slots_per_frame = 10/slot_ms;
schedule.srs_timing = 'special-only';
schedule.srs_ms = srsSlot*slot_ms;
schedule.srs_block_type = tag(srsSlot+1);
schedule.ul_slot = uplinkSlot;
schedule.ul_ms = uplinkSlot*slot_ms;
schedule.srs_id = srsId;
schedule.srs_age_ms = schedule.ul_ms-schedule.srs_ms(srsId);
schedule.ul_duty = numel(uplinkSlot)/numSlots;
end


% ===== Consolidated helper: ul4tCoherenceDrift.m =====
function [phaseDeg,amplitudeDb] = ul4tCoherenceDrift(nPorts,time_ms,p)
%UL4TCOHERENCEDRIFT Sample the attached RF model at the supplied time points.
time_s = time_ms(:).'/1000;
dt = diff(time_s);
assert(~isempty(time_s) && time_s(1)==0 && all(isfinite(time_s)) && all(dt>0), ...
    'RF drift sample times must start at zero and increase strictly.');
scale = {'sigmaPhi0Deg','sigmaA0dB','sigmaDeltaFHz','qPhiDeg2PerSec','sigmaADynDb','tauA'};
rho = {'rhoFreqCommon','rhoPhiCommon','rhoACommon'};
fields = [scale rho];
for i = 1:numel(fields)
    field = fields{i};
    assert(isfield(p,field),'Missing coherence_drift.%s.',field);
    v = p.(field);
    assert(isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v>=0, ...
        'coherence_drift.%s must be a finite nonnegative scalar.',field);
end
assert(p.tauA>0 && max([p.rhoFreqCommon p.rhoPhiCommon p.rhoACommon])<=1, ...
    'Set coherence_drift.tauA > 0 (seconds) and common correlations in [0,1].');

phi0 = p.sigmaPhi0Deg*randn(nPorts,1);
frequencyHz = p.sigmaDeltaFHz*correlatedNormal(nPorts,1,p.rhoFreqCommon);
phaseNoise = correlatedNormal(nPorts,numel(dt),p.rhoPhiCommon);
phaseStep = bsxfun(@times,phaseNoise,sqrt(p.qPhiDeg2PerSec*dt));
phaseDeg = bsxfun(@plus,360*frequencyHz*time_s,phi0) + ...
    [zeros(nPorts,1) cumsum(phaseStep,2)];

amplitude0 = p.sigmaA0dB*randn(nPorts,1);
innovations = p.sigmaADynDb*correlatedNormal(nPorts,numel(time_s),p.rhoACommon);
ou = zeros(nPorts,numel(time_s));
ou(:,1) = innovations(:,1);
decay = exp(-dt/p.tauA);
innovationScale = sqrt(-expm1(-2*dt/p.tauA));
for i = 2:numel(time_s)
    ou(:,i) = decay(i-1)*ou(:,i-1) + innovationScale(i-1)*innovations(:,i);
end
amplitudeDb = bsxfun(@plus,amplitude0,ou);
end

function z = correlatedNormal(nPorts,nSamples,rho)
common = randn(1,nSamples);
z = sqrt(rho)*repmat(common,nPorts,1) + sqrt(1-rho)*randn(nPorts,nSamples);
end


% ===== Consolidated helper: ul4tCalibrationPilots.m =====
function pilot = ul4tCalibrationPilots(cfg)
%UL4TCALIBRATIONPILOTS Orthogonal one-port REs before the PUSCH allocation.
p = cfg.ul.calibration;
nTx = cfg.ul.num_tx;
F = cfg.ofdm.num_tones;
assert(any(strcmp(p.layout,{'frequency','time'})), ...
    'Use frequency or time for the calibration RS layout.');
assert(isnumeric(p.num_re) && isreal(p.num_re) && isscalar(p.num_re) && ...
    isfinite(p.num_re) && p.num_re==fix(p.num_re) && p.num_re>=nTx, ...
    'Calibration RE count must be an integer and cover every TX port.');
prefix = p.prefix_symbols;
if isempty(prefix)
    prefix = 1;
    if strcmp(p.layout,'time'), prefix = min(3,p.num_re); end
end
assert(isscalar(prefix) && ismember(prefix,1:9), ...
    'Reserve 1 to 9 symbols before the PUSCH.');
assert(p.num_re<=prefix*F,'Calibration RE count exceeds the reserved symbols.');
if strcmp(p.layout,'frequency')
    assert(prefix==1,'The frequency layout uses one prefix symbol.');
else
    assert(prefix>=2 && p.num_re>=prefix,'The time layout needs at least two populated symbols.');
end
assert(isnumeric(p.pilot_power) && isreal(p.pilot_power) && isscalar(p.pilot_power) && ...
    isfinite(p.pilot_power) && p.pilot_power>0,'Set positive finite calibration pilot_power.');
assert(isnumeric(p.innovation_variance) && isreal(p.innovation_variance) && ...
    isscalar(p.innovation_variance) && isfinite(p.innovation_variance) && ...
    p.innovation_variance>0,'Set positive finite calibration innovation_variance.');

% Split REs as evenly as possible across prefix symbols, then spread tones.
pilot.tone = zeros(1,p.num_re);
pilot.symbol = zeros(1,p.num_re); % zero-based OFDM symbol
edges = floor((0:prefix)*p.num_re/prefix);
for s = 1:prefix
    indices = edges(s)+1:edges(s+1);
    count = numel(indices);
    pilot.tone(indices) = 1+floor(((0:count-1)+0.5)*F/count);
    pilot.symbol(indices) = s-1;
end
pilot.port = 1+mod(0:p.num_re-1,nTx);
pilot.amplitude = sqrt(p.pilot_power); % known unit-phase pilot
pilot.prefix_symbols = prefix;
pilot.num_re = p.num_re;
pilot.allocated_re = prefix*F;
end


% ===== Consolidated helper: ul4tCalibrationPathBasis.m =====
function pathPilot = ul4tCalibrationPathBasis(model,Hsrs,calIndices,ports,N,F)
%UL4TCALIBRATIONPATHBASIS Decompose ideal SRS into received path components.
% Current PUSCH path gains and Doppler coefficients are never used here.
B = size(model.delayBasis,2);
numRE = numel(ports);
pathPilot = complex(zeros(N*numRE,B));
% Both CDL backends use one delay per cluster in the channel synthesis.
inverseBasis = pinv(model.delayBasis,1e-5*norm(model.delayBasis,'fro'));
for p = 1:size(Hsrs,2)
    mask = ports==p;
    columns = reshape(find(mask),1,[]);
    rows = reshape(bsxfun(@plus,(1:N).',N*(columns-1)),[],1);
    tones = 1+floor(mod(calIndices(1,mask)-1,N*F)/N);
    srsCoefficients = reshape(Hsrs(:,p),N,F)*inverseBasis.';
    for b = 1:B
        contribution = srsCoefficients(:,b)*model.delayBasis(tones,b).';
        pathPilot(rows,b) = contribution(:);
    end
end
end


% ===== Consolidated helper: ul4tTrackCalibration.m =====
function gain = ul4tTrackCalibration(state,received,pathPilot,pilotWeight,ports,noiseVar,priorVariance)
%UL4TTRACKCALIBRATION Joint per-path propagation and relative TX-gain fit.
% A common scalar is unidentifiable. Anchor it to the previous TX state.
N = size(received,1);
nTx = numel(state);
B = size(pathPilot,2);
weight = reshape(repmat(pilotWeight,N,1),[],1);
paths = bsxfun(@times,pathPilot,weight);
rxPort = reshape(repmat(ports,N,1),[],1);
y = received(:);
gain = state;
pathState = ones(B,1);
for iteration = 1:3
    design = bsxfun(@times,paths,gain(rxPort(:)));
    % Ridge centered at the unchanged SRS propagation state.
    pathState = ones(B,1)+(design'*design+noiseVar*eye(B))\(design'*(y-design*ones(B,1)));
    reference = reshape(paths*pathState,N,[]);
    predicted = bsxfun(@times,reference,reshape(state(ports),1,[]));
    innovation = ul4tEstimateCalibration(predicted,received,ports,nTx,noiseVar,priorVariance);
    gain = state.*innovation;
    common = (state'*gain)/(state'*state);
    if isfinite(common) && abs(common)>1e-3
        gain = gain/common;
    else
        gain = state;
    end
end
end


% ===== Consolidated helper: ul4tEstimateCalibration.m =====
function innovation = ul4tEstimateCalibration(reference,received,ports,nTx,noiseVar,priorVariance)
%UL4TESTIMATECALIBRATION Complex LMMSE update centered at the tracked state.
innovation = complex(ones(nTx,1));
for p = 1:nTx
    h = reference(:,ports==p);
    y = received(:,ports==p);
    energy = sum(abs(h(:)).^2);
    update = 1+sum(conj(h(:)).*(y(:)-h(:)))/(energy+noiseVar/priorVariance);
    % Retain the last reliable state if the inverse would be singular.
    if isfinite(update) && abs(update)>1e-3
        innovation(p) = update;
    end
end
end


% ===== Consolidated helper: ul4tResponse.m =====
function H = ul4tResponse(model,time_ms,W)
%UL4TRESPONSE H(rx,frequency,layer). Identity W exposes the selected UE ports.
sample = 1;
if model.max_doppler_Hz>0
    sample = 1+round(time_ms/model.sample_step_ms);
    assert(sample>=1 && sample<=size(model.pathGains,1), ...
        'Requested channel time lies outside the generated sample interval.');
    assert(abs(model.sample_times_ms(sample)-time_ms)<1e-6, ...
        'Requested channel time is off the snapshot grid; explicit interpolation is required.');
end
N = size(model.pathGains,4);
P = size(model.pathGains,2);
nTx = size(model.pathGains,3);
F = size(model.delayBasis,1);
L = size(W,2);
paths = permute(double(model.pathGains(sample,:,:,:)),[4 3 2 1]);
if ~isequal(W,eye(nTx))
    paths = reshape(permute(paths,[1 3 2]),N*P,nTx)*W;
    paths = permute(reshape(paths,N,P,L),[1 3 2]);
end
H = permute(reshape(reshape(paths,N*L,P)*model.weightedDelayBasis.',N,L,F),[1 3 2]);
end


% ===== Consolidated helper: nrUL4TCodebook.m =====
function [W,tpmi,category] = nrUL4TCodebook(rankValue,sets,unitPower)
%NRUL4TCODEBOOK TS 38.211 v17.6.0 Tables 6.3.1.5-3/-5/-6/-7.
% The unnormalized option below retains the literal table scale factors.
if nargin<2 || isempty(sets), sets = {'non','partial','full'}; end
if nargin<3, unitPower = false; end
if ischar(sets), sets = {sets}; end
assert(iscell(sets) && all(ismember(sets,{'non','partial','full'})), ...
    'Use non, partial, full, or a cell array of these categories.');
e = eye(4);
q = [1 1i -1 -1i];
switch rankValue
    case 1
        W = zeros(4,1,28);
        W(:,1,1:4) = reshape(e/2,4,1,4);
        s = [1 -1 1i -1i];
        for n = 1:4
            W(:,1,4+n) = [1;0;s(n);0]/2;
            W(:,1,8+n) = [0;1;0;s(n)]/2;
        end
        for a = 1:4
            for b = 1:4
                last = q(a)*q(b);
                W(:,1,12+(a-1)*4+b) = [1;q(a);q(b);last]/2;
            end
        end
        boundaries = [4 12 28];
    case 2
        W = zeros(4,2,22);
        pair = nchoosek(1:4,2);
        for n = 1:6, W(:,:,n) = e(:,pair(n,:))/2; end
        a = [1 1 -1i -1i -1 -1 1i 1i];
        b = [-1i 1i 1 -1 -1i 1i 1 -1];
        for n = 1:8, W(:,:,6+n) = [1 0;0 1;a(n) 0;0 b(n)]/2; end
        for n = 1:4
            W(:,:,13+2*n) = [1 1;q(n) q(n);1 -1;q(n) -q(n)]/sqrt(8);
            W(:,:,14+2*n) = [1 1;q(n) q(n);1i -1i;1i*q(n) -1i*q(n)]/sqrt(8);
        end
        boundaries = [6 14 22];
    case 3
        W = zeros(4,3,7);
        W(:,:,1) = e(:,1:3)/2;
        W(:,:,2) = [1 0 0;0 1 0;1 0 0;0 0 1]/2;
        W(:,:,3) = [1 0 0;0 1 0;-1 0 0;0 0 1]/2;
        for a = 0:1
            for b = 0:1
                signA = (-1)^a;
                phaseB = 1i^b;
                W(:,:,4+2*a+b) = [1 1 1;signA*[1 -1 1]; ...
                    phaseB*[1 1 -1];signA*phaseB*[1 -1 -1]]/sqrt(12);
            end
        end
        boundaries = [1 3 7];
    case 4
        W = zeros(4,4,5);
        W(:,:,1) = e/2;
        W(:,:,2) = [1 1 0 0;0 0 1 1;1 -1 0 0;0 0 1 -1]/sqrt(8);
        W(:,:,3) = diag([1 1 1i 1i])*W(:,:,2);
        W(:,:,4) = [1 1 1 1;1 -1 1 -1;1 1 -1 -1;1 -1 -1 1]/4;
        W(:,:,5) = diag([1 1 1i 1i])*W(:,:,4);
        boundaries = [1 3 5];
    otherwise
        error('rank must be one of 1, 2, 3 or 4.');
end
category = [repmat({'non'},1,boundaries(1)), ...
    repmat({'partial'},1,boundaries(2)-boundaries(1)), ...
    repmat({'full'},1,boundaries(3)-boundaries(2))];
keep = ismember(category,sets);
assert(any(keep),'Empty codebook subset.');
tpmi = find(keep)-1;
W = W(:,:,keep);
category = category(keep);
if unitPower
    for n = 1:size(W,3), W(:,:,n) = W(:,:,n)/norm(W(:,:,n),'fro'); end
end
end


% ===== Consolidated helper: nrUL2TCodebook.m =====
function [W,tpmi,category] = nrUL2TCodebook(rankValue,sets,unitPower)
%NRUL2TCODEBOOK TS 38.211 v17.6.0 Tables 6.3.1.5-1 and 6.3.1.5-4.
% The default retains the table scale factors.
if nargin<2 || isempty(sets), sets = {'non','full'}; end
if nargin<3, unitPower = false; end
if ischar(sets), sets = {sets}; end
assert(iscell(sets) && all(ismember(sets,{'non','full'})), ...
    'Use non, full, or a cell array of these categories for two UE ports.');
switch rankValue
    case 1
        W = reshape([1 0 1 1 1 1;0 1 1 -1 1i -1i]/sqrt(2),2,1,6);
        category = {'non','non','full','full','full','full'};
    case 2
        W = cat(3,eye(2)/sqrt(2),[1 1;1 -1]/2,[1 1;1i -1i]/2);
        category = {'non','full','full'};
    otherwise
        error('For two UE ports, rank must be 1 or 2.');
end
keep = ismember(category,sets);
tpmi = find(keep)-1;
W = W(:,:,keep);
category = category(keep);
if unitPower
    for n = 1:size(W,3), W(:,:,n) = W(:,:,n)/norm(W(:,:,n),'fro'); end
end
end
