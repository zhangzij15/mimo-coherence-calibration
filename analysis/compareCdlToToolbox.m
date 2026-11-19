function [summary,details] = compareCdlToToolbox(seedCount,speedKmh,outputDir)
%COMPARECDLTOTOOLBOX Compare the project CDL backends at matched inputs.
% Run with the 5G Toolbox installed. The channel generator itself needs only MATLAB.
arguments
    seedCount (1,1) double {mustBeInteger,mustBePositive} = 80
    speedKmh (1,1) double {mustBePositive} = 30
    outputDir (1,1) string = fullfile(fileparts(mfilename('fullpath')),"cdl_toolbox_comparison")
end
assert(~isempty(which('nrCDLChannel')),'This comparison requires 5G Toolbox.');
savedRng = rng;
restoreRng = onCleanup(@()rng(savedRng));
if ~isfolder(outputDir), mkdir(outputDir); end

cfg = configUL4T('base');
cfg.users.mode = 'specified';
cfg.users.positions_m = [0 60 1.5];
cfg.ul.ue_speed_kmh = speedKmh;
cfg.seed = 26011201;
cfg = prepareFixedArrayConfig(cfg);
pattern = fixedArrayChannel('pattern',cfg);
ue = struct('layout','upa','positions_lambda',.22*[0 0 0;0 0 1;1 0 0;1 0 1], ...
    'antenna_indices',1:4,'toolbox_ports',1:4);
% Observe Doppler over 40 ms and report the lags relevant to SRS aging.
schedule = struct('srs_ms',0,'ul_ms',(0:1:40).','slot_ms',1);
lagMs = [1 2 4 8 12 16 20];
profiles = ["CDL-A";"CDL-B";"CDL-C";"CDL-D";"CDL-E"];
% Sample spatial covariance from 16 ports within the first nine rows and columns.
[row,col] = ndgrid([1 2 5 9],[1 2 5 9]);
bsPorts = sub2ind([cfg.Nv cfg.Nh],row(:),col(:));

summary = table();
details = struct([]);
for k = 1:numel(profiles)
    cfg.cdl.toolbox_profile = char(profiles(k));
    exactMatch = true;
    sampleErrorEnergy = 0;
    sampleRefEnergy = 0;
    for mc = 1:seedCount
        rng(cfg.seed+mc-1,'twister');
        [native,nativeScene] = ul4tNativeCDL(cfg,pattern,schedule,mc,ue);
        rng(cfg.seed+mc-1,'twister');
        [ref,refScene] = ul4tToolboxChannel(cfg,pattern,schedule,mc,ue);
        assert(isequal(nativeScene.positions_m,refScene.positions_m));
        assert(isequal(size(native.pathGains),size(ref.pathGains)));
        assert(max(abs(native.sample_times_ms-ref.sample_times_ms))<1e-7);
        assert(native.max_doppler_Hz==ref.max_doppler_Hz);
        lagSamples = round(lagMs/native.sample_step_ms);
        if mc==1
            nTap = size(native.pathGains,2);
            a = emptyStats(nTap,numel(bsPorts),numel(lagSamples));
            b = a;
            seedPower = zeros(seedCount,2);
        end
        exactMatch = exactMatch && isequal(native.pathGains,ref.pathGains);
        difference = native.pathGains-ref.pathGains;
        sampleErrorEnergy = sampleErrorEnergy + sum(abs(difference).^2,'all','double');
        sampleRefEnergy = sampleRefEnergy + sum(abs(ref.pathGains).^2,'all','double');
        nativeStats = measure(native.pathGains,native.diagnostics.cluster_pattern_power,bsPorts,lagSamples);
        refStats = measure(ref.pathGains,ref.diagnostics.cluster_pattern_power,bsPorts,lagSamples);
        seedPower(mc,:) = [sum(nativeStats.tapPower),sum(refStats.tapPower)];
        a = addStats(a,nativeStats);
        b = addStats(b,refStats);
    end

    delay = native.diagnostics.path_delays_s(:);
    nominal = native.diagnostics.cluster_powers_linear(:);
    profileRef = nrCDLChannel('DelayProfile','Custom','InitialDelayProfile',profiles(k), ...
        'DelaySpread',cfg.cdl.toolbox_delay_spread_s);
    refNominal = 10.^(profileRef.AveragePathGains(:)/10);
    refNominal = refNominal/sum(refNominal);
    assert(native.diagnostics.has_los==ref.diagnostics.has_los);
    assert(native.diagnostics.first_cluster_k_factor_dB==ref.diagnostics.first_cluster_k_factor_dB);
    pdp = [a.tapPower/sum(a.tapPower),b.tapPower/sum(b.tapPower)];
    ds = sqrt(sum(pdp.*(delay-sum(pdp.*delay,1)).^2,1));
    timeA = normalizedTime(a.time);
    timeB = normalizedTime(b.time);
    txDelta = matrixError(a.tx,b.tx);
    rxDelta = matrixError(a.rx,b.rx);
    [~,lag4] = min(abs(lagMs-4));
    [~,lag20] = min(abs(lagMs-20));
    summary = [summary;table(profiles(k),nTap, ...
        max(abs(delay-ref.diagnostics.path_delays_s(:)))*1e9,max(abs(nominal-refNominal)), ...
        exactMatch,sqrt(sampleErrorEnergy/sampleRefEnergy), ...
        10*log10(sum(a.tapPower)/sum(b.tapPower)),50*sum(abs(diff(pdp,1,2))), ...
        max(abs(10*log10(pdp(:,1)./pdp(:,2)))),ds(1)*1e9,ds(2)*1e9, ...
        100*txDelta,100*rxDelta,abs(timeA(lag4)-timeB(lag4)), ...
        abs(timeA(lag20)-timeB(lag20)),max(abs(timeA-timeB)), ...
        'VariableNames',{'Profile','Taps','MaxDelayDiff_ns','MaxNominalPowerDiff', ...
        'IdenticalAtSameSeed','SameSeedRelL2','MeanPowerDiff_dB','PdpTotalVariation_pct', ...
        'PdpMaxTapDiff_dB','BaseRmsDelay_ns','ToolboxRmsDelay_ns','UeCovDiff_pct', ...
        'BsCovDiff_pct','TimeCorr4msDiff','TimeCorr20msDiff','MaxTimeCorrDiff'})]; %#ok<AGROW>
    details(k).profile = profiles(k);
    details(k).delay_s = delay;
    details(k).pdp = [nominal,pdp];
    details(k).ue_covariance = cat(3,normalizeCov(a.tx),normalizeCov(b.tx));
    details(k).bs_subarray_covariance = cat(3,normalizeCov(a.rx),normalizeCov(b.rx));
    details(k).lag_ms = lagMs;
    details(k).time_correlation = [timeA,timeB];
    details(k).seed_power = seedPower;
    fprintf('%s finished, %d seeds, %g km/h, same-seed L2 %.3f, PDP TV %.3f%%, Tx covariance %.3f%%, time corr %.4f\n', ...
        profiles(k),seedCount,speedKmh,summary.SameSeedRelL2(end), ...
        summary.PdpTotalVariation_pct(end),summary.UeCovDiff_pct(end),summary.MaxTimeCorrDiff(end));
end
metadata = struct('seed_count',seedCount,'speed_kmh',speedKmh,'first_seed',cfg.seed, ...
    'matlab_version',version,'toolbox_version',ver('5g'),'schedule',schedule, ...
    'config',cfg,'ue',ue,'bs_covariance_ports',bsPorts);
stem = sprintf('cdl_comparison_%gkmh_%dseeds',speedKmh,seedCount);
writetable(summary,fullfile(outputDir,[stem '.csv']));
save(fullfile(outputDir,[stem '.mat']),'summary','details','metadata');
disp(summary);
end

function stats = emptyStats(nTap,nBs,nLag)
stats = struct('tapPower',zeros(nTap,1),'tx',zeros(4),'rx',zeros(nBs), ...
    'time',zeros(nLag,3));
end

function stats = measure(gains,patternPower,bsPorts,lags)
[nTime,nTap,nTx,~] = size(gains);
stats.tapPower = reshape(sum(abs(gains).^2,[1 3 4],'double'),[],1)/numel(gains)*nTap;
% The physical-array statistics include the project's external BS pattern.
weighted = gains.*reshape(single(sqrt(patternPower)),1,nTap);
txData = reshape(permute(weighted,[1 2 4 3]),[],nTx);
rxData = reshape(weighted(:,:,:,bsPorts),[],numel(bsPorts));
stats.tx = double(txData'*txData);
stats.rx = double(rxData'*rxData);
lagData = reshape(weighted,nTime,[]);
powerPerTime = sum(abs(lagData).^2,2,'double');
stats.time = zeros(numel(lags),3);
for j = 1:numel(lags)
    d = lags(j);
    stats.time(j,:) = [sum(lagData(1+d:end,:).*conj(lagData(1:end-d,:)),'all','double'), ...
        sum(powerPerTime(1:end-d)),sum(powerPerTime(1+d:end))];
end
end

function total = addStats(total,increment)
for f = fieldnames(total).'
    name = f{1};
    total.(name) = total.(name)+increment.(name);
end
end

function R = normalizeCov(G)
R = G./sqrt(real(diag(G))*real(diag(G)).');
end

function e = matrixError(a,b)
a = normalizeCov(a);
b = normalizeCov(b);
e = norm(a-b,'fro')/norm(b,'fro');
end

function r = normalizedTime(x)
r = x(:,1)./sqrt(real(x(:,2)).*real(x(:,3)));
end
