function report = runCoherence256Search(opts)
%RUNCOHERENCE256SEARCH Screen 4T256R cases then confirm on new MC seeds.
%   report = runCoherence256Search(opts) accepts a scalar options struct.
%   dryRun returns configurations without starting simulation or writing files.
%   Set stage to 'screen', 'confirm', or 'both'. Keep the screen options and
%   outputDir when resuming. A distinct outputDir is needed for concurrent runs.
if nargin==0, opts=struct(); end
root=fileparts(mfilename('fullpath'));
opts=normalizeOptions(opts,root);
scenarios=buildScenarios();
scenarios=selectNames(scenarios,opts.scenarioNames);
thresholds=struct('snr_dB',-20,'loss_dB',.5,'gain_dB',.4, ...
    'phase20_deg',40,'amp20_dB',4);
[sourceHash,sourceManifest]=fingerprint(root,opts.dryRun);
configs=cell(1,numel(scenarios));
for i=1:numel(scenarios)
    configs{i}=makeConfig(scenarios(i),opts,opts.screenDurationMs, ...
        opts.screenSNRdB,opts.screenSeed,opts.screenMC);
end
report=struct('schema','coherence256_search_v1','status','planned','options',opts, ...
    'source_sha256',sourceHash,'source_manifest',sourceManifest,'thresholds',thresholds, ...
    'plan',struct('scenarios',scenarios,'screen_configs',{configs}), ...
    'case_order',{{'Ideal','Uncal','Feedback','Reselect','Oracle'}}, ...
    'age_groups',{{'all','le20','gt20'}}, ...
    'selection_basis','Automatic selection ranks NR only. SVD rows are a separate sensitivity study.', ...
    'receiver_definition','Single-user linear MMSE-IRC with white noise covariance sigma2 I.', ...
    'throughput_definition','Main frame-average throughput without deducted calibration overhead.', ...
    'bandwidth_definition','288 RB research grid, 103.68 MHz at 30 kHz, contiguous base MATLAB NR DM-RS.', ...
    'rf_definition',['Four per-port exact 20 ms maxima and six port-pair maxima across all sampled lags up to 20 ms. ' ...
    'Both per-port and pairwise changes are also audited over all sampled lags up to 20 ms. ' ...
    'SRS-relative changes are grouped by actual age. Phase slew may exceed the limit at longer ages.']);
if opts.dryRun, return; end

if ~isfolder(opts.outputDir), mkdir(opts.outputDir); end
fid=fopen(fullfile(opts.outputDir,'run.log'),'a','n','UTF-8');
assert(fid>0,'UL256:Log','Cannot open run log.');
closeLog=onCleanup(@()fclose(fid)); %#ok<NASGU>
logLine(fid,'START %s, stage=%s, source=%s\n',timestamp(),opts.stage,sourceHash);
logLine(fid,'Channel profile: %s.\n',opts.channelProfile);
logLine(fid,'288 physical RB research grid with contiguous base MATLAB DM-RS.\n');
logLine(fid,'Screen -20 dB overall targets: loss >= %.2f dB, gain >= %.2f dB, both paired CI95 lower > 0.\n', ...
    thresholds.loss_dB,thresholds.gain_dB);
logLine(fid,'Frame THP excludes deduction of calibration overhead. SNR per-antenna reference uses white noise.\n');
writeTableSafe(struct2table(scenarios(:),'AsArray',true),fullfile(opts.outputDir,'scenario_manifest.csv'),fid);
report.status='running';
saveAtomic(fullfile(opts.outputDir,'plan.mat'),'plan',report);

screen=initializeStage('screen',scenarios,configs,opts,sourceHash,fid);
if ~strcmp(opts.stage,'confirm')
    screen=runStage(screen,opts,thresholds,fid);
else
    assert(all(ismember({screen.status},{'complete','failed'})) && any(strcmp({screen.status},'complete')), ...
        'UL256:ScreenRequired','Confirm requires matching screen checkpoints for all scenarios.');
end
report.screen=screen;
report.selection=selectCandidate(screen,opts,thresholds);
writeTableSafe(struct2table(report.selection(:),'AsArray',true),fullfile(opts.outputDir,'selection.csv'),fid);
report.confirm=[];
if ~strcmp(opts.stage,'screen')
    if isempty(report.selection)
        logLine(fid,'No completed non-control screening candidate. Confirmation is unavailable.\n');
    else
        [confirmScenarios,lockedMethods]=makeConfirmation(report.selection,scenarios);
        confirmConfigs=cell(1,numel(confirmScenarios));
        for i=1:numel(confirmScenarios)
            confirmConfigs{i}=makeConfig(confirmScenarios(i),opts,opts.finalDurationMs, ...
                opts.finalSNRdB,opts.finalSeed,opts.finalMC);
        end
        confirm=initializeStage('confirm',confirmScenarios,confirmConfigs,opts,sourceHash,fid);
        for i=1:numel(confirm), confirm(i).locked_method=lockedMethods{i}; end
        writeTableSafe(struct2table(confirmScenarios(:),'AsArray',true),fullfile(opts.outputDir,'confirmation_manifest.csv'),fid);
        confirm=runStage(confirm,opts,thresholds,fid);
        report.confirm=confirm;
    end
end
allStages=report.screen;
if ~isempty(report.confirm), allStages=[allStages,report.confirm]; end
report.status='complete';
if any(strcmp({allStages.status},'failed')) || (~strcmp(opts.stage,'screen') && isempty(report.confirm))
    report.status='complete_with_failures';
end
report.completed_at=timestamp();
saveAtomic(fullfile(opts.outputDir,'report.mat'),'report',report);
logLine(fid,'Campaign %s at %s.\n',upper(report.status),report.completed_at);
end

function o = normalizeOptions(in,root)
assert(isstruct(in)&&isscalar(in),'UL256:Options','Provide a scalar options struct.');
o=struct('dryRun',false,'stage','both','screenMC',4,'finalMC',12, ...
    'screenSNRdB',-20,'finalSNRdB',[-22 -20 -18], ...
    'screenDurationMs',50,'finalDurationMs',100,'screenSeed',202611181,'finalSeed',203611181, ...
    'outputDir','','scenarioNames',{{}},'confirmNames',{{}}, ...
    'calibrationRE',128,'resume',true,'channelProfile','CDL-C');
names=fieldnames(in);
assert(all(ismember(names,fieldnames(o))),'UL256:Options','Unknown options field.');
for i=1:numel(names), o.(names{i})=in.(names{i}); end
validateattributes(o.dryRun,{'logical','numeric'},{'scalar','binary'});
validateattributes(o.resume,{'logical','numeric'},{'scalar','binary'});
o.stage=validatestring(o.stage,{'both','screen','confirm'});
o.channelProfile=validatestring(o.channelProfile,{'CDL-A','CDL-B','CDL-C','CDL-D','CDL-E'});
if isempty(o.outputDir)
    o.outputDir=fullfile(root,['results_coherence256_' lower(o.channelProfile)]);
end
for name={'screenMC','finalMC','screenSeed','finalSeed','calibrationRE'}
    validateattributes(o.(name{1}),{'numeric'},{'scalar','integer','finite','positive'},mfilename,name{1});
end
assert(o.calibrationRE>=4 && o.calibrationRE<=3456,'UL256:Options','calibrationRE must be 4–3456.');
assert(o.finalSeed>o.screenSeed+o.screenMC-1 && o.finalSeed+o.finalMC-1<2^32, ...
    'UL256:Options','Choose distinct in-range screen and confirmation seed blocks.');
for name={'screenDurationMs','finalDurationMs'}
    validateattributes(o.(name{1}),{'numeric'},{'scalar','finite','>=',10.5},mfilename,name{1});
    assert(o.(name{1})*2==round(o.(name{1})*2),'UL256:Options','Durations must be multiples of 0.5 ms.');
end
for name={'screenSNRdB','finalSNRdB'}
    validateattributes(o.(name{1}),{'numeric'},{'vector','real','finite','nonempty'},mfilename,name{1});
    o.(name{1})=sort(reshape(double(o.(name{1})),1,[]));
    assert(numel(unique(o.(name{1})))==numel(o.(name{1})),'UL256:Options','Duplicate SNR points.');
end
if ~strcmp(o.stage,'screen')
    assert(any(o.screenSNRdB==-20) && any(o.finalSNRdB==-20),'UL256:Options', ...
        'Screen and confirmation both need the predeclared -20 dB endpoint.');
end
for name={'scenarioNames','confirmNames'}
    o.(name{1})=cellstr(string(o.(name{1})));
    assert(numel(unique(o.(name{1})))==numel(o.(name{1})),'UL256:Options','Duplicate scenario names.');
end
assert(numel(o.confirmNames)<=2,'UL256:Options','Select at most two confirmation candidates.');
o.outputDir=char(java.io.File(char(o.outputDir)).getCanonicalPath());
end

function s = buildScenarios()
baseline=struct('name','PA_LOAD_COMPACT_R1','layout','upa','spacing_lambda',.35,'precoder','nr', ...
    'speed_kmh',.3,'rank',1,'mode','load','grouping','balanced','phase20_deg',38.2, ...
    'amp20_dB',3.8,'tau_ms',1.5,'dv_lambda',1.5,'dh_lambda',.5,'is_control',false);
s=repmat(baseline,1,13);
s(2).name='PA_LOAD_UPA05_R1'; s(2).spacing_lambda=.5;
s(3).name='PA_LOAD_HANDHELD_R1'; s(3).layout='handheld'; s(3).spacing_lambda=NaN;
s(4).name='SLEW_COMPACT_R1'; s(4).mode='slew';
s(5).name='SLEW_SPREAD_R1'; s(5).mode='slew'; s(5).grouping='spread';
s(6).name='BOUNDED_COMPACT_R1'; s(6).mode='bounded';
s(7).name='PA_LOAD_COMPACT_3KMH_R1'; s(7).speed_kmh=3;
s(8).name='PA_LOAD_COMPACT_HALF_WAVE_RX_R1'; s(8).dv_lambda=.5;
s(9).name='PA_LOAD_COMPACT_R2'; s(9).rank=2;
s(10).name='ZERO_RF_COMPACT_R1'; s(10).phase20_deg=0; s(10).amp20_dB=0; s(10).is_control=true;
s(11).name='SLEW_COMPACT_SVD_R1'; s(11).mode='slew'; s(11).precoder='svd';
s(12).name='PA_LOAD_COMPACT_SVD_R1'; s(12).precoder='svd';
s(13).name='PA_LOAD_FAST_SVD_R1'; s(13).precoder='svd'; s(13).tau_ms=.5;
end

function out = selectNames(s,names)
if isempty(names), out=s; return; end
[found,idx]=ismember(names,{s.name});
assert(all(found),'UL256:Scenario','Unrecognized scenario name.');
out=reshape(s(idx),1,[]);
end

function cfg = makeConfig(s,o,duration,snr,seed,nMC)
cfg=configUL4T('base');
cfg.Nh=16; cfg.Nv=16; cfg.dh_lambda=s.dh_lambda; cfg.dv_lambda=s.dv_lambda;
cfg.K=1; cfg.fc_Hz=6.7e9;
cfg.ofdm.scs_Hz=30e3; cfg.ofdm.num_rb=288; cfg.ofdm.nfft=5120;
cfg.receiver='mmse-irc'; cfg.irc.external_covariance=[];
cfg.cdl.toolbox_profile=o.channelProfile; cfg.cdl.toolbox_delay_spread_s=100e-9;
cfg.MC=nMC; cfg.seed=seed; cfg.progress_every=1; cfg.output_dir=o.outputDir;
cfg.save_results=false; cfg.output.save_per_tone=false; cfg.output.save_realizations=false;
cfg.ul.num_tx=4; cfg.ul.rank=s.rank; cfg.ul.ue_layout=s.layout;
cfg.ul.ue_speed_kmh=s.speed_kmh; cfg.ul.channel_aging_enabled=true;
if strcmp(s.layout,'upa')
    cfg.ul.tx_positions_lambda=s.spacing_lambda*[0 0 0;0 0 1;1 0 0;1 0 1];
end
cfg.ul.tdd_pattern='DSUDD'; cfg.ul.block_ms=5; cfg.ul.srs_period_ms=40;
cfg.ul.duration_ms=duration; cfg.ul.rx_noise_power=10.^(-snr/10);
cfg.ul.dmrs_segments_rb=[]; cfg.ul.calibration.num_re=o.calibrationRE;
cfg.ul.coherence_enabled=true; cfg.ul.coherence_model='drift';
cfg.ul.coherence_drift.bounded=struct('mode',s.mode,'phase20_deg',s.phase20_deg, ...
    'amp20_dB',s.amp20_dB,'grouping',s.grouping,'tau_ms',s.tau_ms);
cfg.ul.cases=struct('precoder',repmat({s.precoder},1,5),'coherence_error',{false,true,true,true,true}, ...
    'sets',repmat({{'full'}},1,5),'calibration',{'none','none','feedback','reselect','oracle'});
end

function [h,man] = fingerprint(root,dryRun)
names={ 'awgn_bler10_t1.csv','configFixedArray.m','configUL4T.m','fixedArrayChannel.m', ...
    'linkScene.m','prepareFixedArrayConfig.m','runCoherence256Search.m','simulateUL4T.m', ...
    'ul4tBoundedRfDrift.m','ul4tLinkMapping.m','ul4tPostSinr.m','ul4tRunnerCacheKey.m', ...
    'generateCdlChannel.m','ul4tNativeCDL.m','ul4tBaseDMRS.m','ul4tBaseTBS.m', ...
    'ul4tToolboxChannel.m','ul4tToolboxDMRS.m' };
man=struct('name',names,'sha256',repmat({''},size(names)));
for i=1:numel(man)
    path=fullfile(root,names{i});
    if ~isfile(path) && dryRun, man(i).sha256='MISSING'; continue; end
    fid=fopen(path,'rb');
    assert(fid>0,'UL256:Source','Cannot hash production input %s.',names{i});
    cleanup=onCleanup(@()fclose(fid));
    b=fread(fid,Inf,'*uint8');
    clear cleanup
    digest=java.security.MessageDigest.getInstance('SHA-256');
    digest.update(b);
    man(i).sha256=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
h=ul4tRunnerCacheKey(man,'coherence256-production-sha256');
end

function cp = initializeStage(stage,scenarios,configs,o,sourceHash,fid)
base=struct('schema','coherence256_checkpoint_v1','key','','source_sha256',sourceHash, ...
    'scenario',struct(),'stage',stage,'locked_method','','cfg',struct(),'checkpoint_file','', ...
    'completed_mc',0,'requested_mc',0,'current_mc',0,'status','pending','entries',{{}}, ...
    'error_identifier','','error_message','');
cp=repmat(base,1,numel(scenarios));
for i=1:numel(scenarios)
    cfg=configs{i};
    expected=base;
    expected.scenario=scenarios(i); expected.cfg=cfg; expected.requested_mc=cfg.MC;
    expected.key=ul4tRunnerCacheKey(struct('stage',stage,'scenario',scenarios(i),'cfg',cfg),sourceHash);
    expected.checkpoint_file=fullfile(o.outputDir,sprintf('%s_%s_%s.mat', ...
        stage,scenarios(i).name,expected.key(1:12)));
    cp(i)=expected;
    if ~(o.resume||strcmp(o.stage,'confirm')) || ~isfile(expected.checkpoint_file), continue; end
    try
        file=load(expected.checkpoint_file,'checkpoint');
        c=file.checkpoint;
        assert(validCheckpoint(c,expected),'UL256:Resume','Checkpoint validation failed.');
        cp(i)=c;
        logLine(fid,'RESUME %s %s %d/%d.\n',stage,c.scenario.name,c.completed_mc,c.requested_mc);
    catch ME
        logLine(fid,'Ignoring %s because %s\n',expected.checkpoint_file,ME.message);
    end
end
end

function valid = validCheckpoint(c,e)
valid=isequal(sort(fieldnames(c)),sort(fieldnames(e))) && strcmp(c.schema,e.schema) && ...
    strcmp(c.key,e.key) && strcmp(c.source_sha256,e.source_sha256) && ...
    isequaln(c.cfg,e.cfg) && isequaln(c.scenario,e.scenario) && strcmp(c.stage,e.stage) && ...
    strcmp(c.checkpoint_file,e.checkpoint_file) && c.requested_mc==e.requested_mc && ...
    ismember(c.completed_mc,0:e.requested_mc) && iscell(c.entries) && ...
    numel(c.entries)==c.completed_mc && any(strcmp(c.status,{'pending','running','complete','failed'}));
if ~valid, return; end
needSize=[3,numel(e.cfg.ul.rx_noise_power),5];
for i=1:c.completed_mc
    a=c.entries{i};
    valid=isstruct(a) && isscalar(a) && a.mc_index==i && a.seed==e.cfg.seed+i-1 && ...
        isequal(size(a.mean_sinr_dB),needSize) && all(isfinite(a.mean_sinr_dB(1,:,:)),'all');
    if ~valid, return; end
end
end

function stage = runStage(stage,o,thresholds,fid)
exportStage(stage,o,thresholds,fid);
for i=1:numel(stage)
    c=stage(i);
    if strcmp(c.status,'complete'), continue; end
    c.error_identifier=''; c.error_message='';
    try
        for j=c.completed_mc+1:c.requested_mc
            c.current_mc=j; c.status='running';
            saveAtomic(c.checkpoint_file,'checkpoint',c);
            stage(i)=c;
            exportStatus(stage,o,fid);
            cfg=c.cfg; cfg.MC=1; cfg.seed=c.cfg.seed+j-1;
            logLine(fid,'%s %s MC %d/%d START seed=%d.\n',upper(c.stage),c.scenario.name,j,c.requested_mc,cfg.seed);
            timer=tic;
            r=simulateUL4T(cfg);
            a=summarizeMC(r,c.scenario,cfg.seed,j,thresholds);
            a.elapsed_s=toc(timer);
            c.entries{j}=a; c.completed_mc=j;
            saveAtomic(c.checkpoint_file,'checkpoint',c);
            stage(i)=c;
            exportStage(stage,o,thresholds,fid);
            logLine(fid,'%s %s MC %d/%d DONE %.2fs.\n',upper(c.stage),c.scenario.name,j,c.requested_mc,a.elapsed_s);
        end
        c.status='complete'; c.current_mc=0;
        saveAtomic(c.checkpoint_file,'checkpoint',c);
    catch ME
        c.status='failed'; c.error_identifier=ME.identifier; c.error_message=ME.message;
        saveAtomic(c.checkpoint_file,'checkpoint',c);
        logLine(fid,'FAILED %s %s MC %d: %s\n',c.stage,c.scenario.name,c.current_mc, ...
            getReport(ME,'extended','hyperlinks','off'));
    end
    stage(i)=c;
    exportStage(stage,o,thresholds,fid);
end
end

function a = summarizeMC(r,s,seed,index,t)
assert(r.fixed_beta_linear==1 && size(r.fixed_beta_eff_sinr_dB,4)==5, ...
    'UL256:Cases','Expected fixed beta 1 and five cases.');
assert(r.cfg.N==256 && r.cfg.ofdm.num_tones==3456 && r.cfg.ul.num_tx==4, ...
    'UL256:Grid','Expected 4T256R and 3456 occupied tones.');
ns=numel(r.rx_noise_power); slots=numel(r.schedule.ul_ms);
sinr=reshape(r.fixed_beta_eff_sinr_dB,slots,ns,5);
bits=reshape(r.good_bits,slots,ns,5);
assert(all(isfinite(sinr(:))) && all(isfinite(bits(:))),'UL256:Finite','Nonfinite link metric.');
groups={true(slots,1),r.schedule.srs_age_ms(:)<=20,r.schedule.srs_age_ms(:)>20};
gain=NaN(3,ns,5); slotTHP=NaN(3,ns,5); counts=zeros(1,3);
for i=1:3
    take=groups{i}; counts(i)=nnz(take);
    if counts(i)==0, continue; end
    gain(i,:,:)=mean(sinr(take,:,:),1);
    slotTHP(i,:,:)=mean(bits(take,:,:),1)/(1000*r.schedule.slot_ms);
end
d=r.coherence_drift;
phase=reshape(d.phase_deg,4,[]); amp=reshape(d.amplitude_dB,4,[]);
lag=round(20/r.schedule.slot_ms);
if size(phase,2)>lag
    absPhase=abs(wrapDegrees(phase(:,1+lag:end)-phase(:,1:end-lag)));
    absAmp=abs(amp(:,1+lag:end)-amp(:,1:end-lag));
    phase20=max(absPhase,[],2).'; amp20=max(absAmp,[],2).';
else
    phase20=NaN(1,4); amp20=NaN(1,4);
end
portPairs=nchoosek(1:4,2);
pairPhase=phase(portPairs(:,1),:)-phase(portPairs(:,2),:);
pairAmp=amp(portPairs(:,1),:)-amp(portPairs(:,2),:);
[~,~,over20]=slidingDrift(phase,amp,lag,t);
[pairPhase20,pairAmp20,pairOver20]=slidingDrift(pairPhase,pairAmp,lag,t);
assert(over20==0 && pairOver20==0,'UL256:RFEnvelope', ...
    'Per-port or port-pair drift violates the envelope at a lag up to 20 ms.');
idx=d.pusch_sample_index; ref=d.srs_sample_index(r.schedule.srs_id);
phaseSRS=abs(wrapDegrees(phase(:,idx)-phase(:,ref)));
ampSRS=abs(amp(:,idx)-amp(:,ref));
pairPhaseSRS=abs(wrapDegrees(pairPhase(:,idx)-pairPhase(:,ref)));
pairAmpSRS=abs(pairAmp(:,idx)-pairAmp(:,ref));
srsPhaseMax=NaN(3,4); srsAmpMax=NaN(3,4);
pairSrsPhaseMax=NaN(3,6); pairSrsAmpMax=NaN(3,6);
for i=1:3
    if counts(i)>0
        srsPhaseMax(i,:)=max(phaseSRS(:,groups{i}),[],2).';
        srsAmpMax(i,:)=max(ampSRS(:,groups{i}),[],2).';
        pairSrsPhaseMax(i,:)=max(pairPhaseSRS(:,groups{i}),[],2).';
        pairSrsAmpMax(i,:)=max(pairAmpSRS(:,groups{i}),[],2).';
    end
end
short=groups{2};
if any(short)
    assert(all(phaseSRS(:,short)<=t.phase20_deg+1e-8,'all') && all(ampSRS(:,short)<=t.amp20_dB+1e-8,'all') && ...
        all(pairPhaseSRS(:,short)<=t.phase20_deg+1e-8,'all') && all(pairAmpSRS(:,short)<=t.amp20_dB+1e-8,'all'), ...
        'UL256:ShortAgeEnvelope','SRS-relative RF drift exceeds the limit at age <= 20 ms.');
end
if ~strcmp(s.mode,'slew')
    assert(all(pairPhaseSRS<=t.phase20_deg+1e-8,'all') && all(pairAmpSRS<=t.amp20_dB+1e-8,'all'), ...
        'UL256:GlobalEnvelope','Load or bounded mode exceeds the pairwise RF envelope at an actual SRS age.');
end
energy=reshape(r.pusch_effective_energy,slots,5);
oracleEnergyError=max(abs(energy(:,5)-energy(:,1))./max(energy(:,1),realmin));
powerError=max(abs(r.pusch_tx_power(:)-1));
assert(powerError<1e-10 && oracleEnergyError<1e-9,'UL256:Invariant','TX power or Oracle effective energy mismatch.');
if s.is_control
    assert(all(abs(sinr(:,:,1)-sinr(:,:,2))<1e-10,'all') && all(bits(:,:,1)==bits(:,:,2),'all'), ...
        'UL256:ZeroRF','Ideal and Uncal must agree for zero drift.');
end
a=struct('mc_index',index,'seed',seed,'snr_dB',reshape(10*log10(1./r.rx_noise_power),1,[]), ...
    'mean_sinr_dB',gain,'ul_slot_goodput_Mbps',slotTHP, ...
    'frame_thp_Mbps',reshape(r.thp_frame_Mbps,ns,5),'slots_per_group',counts, ...
    'phase20_abs_max_deg',phase20,'amp20_abs_max_dB',amp20,'rf20_violation_count',over20, ...
    'srs_phase_abs_max_deg',srsPhaseMax,'srs_amp_abs_max_dB',srsAmpMax, ...
    'port_pairs',portPairs,'pair_phase20_abs_max_deg',pairPhase20,'pair_amp20_abs_max_dB',pairAmp20, ...
    'pair_rf20_violation_count',pairOver20,'pair_srs_phase_abs_max_deg',pairSrsPhaseMax, ...
    'pair_srs_amp_abs_max_dB',pairSrsAmpMax, ...
    'phaseSRS_abs_deg',phaseSRS,'ampSRS_abs_dB',ampSRS,'schedule',r.schedule, ...
    'tx_power_error',powerError,'oracle_relative_energy_error',oracleEnergyError, ...
    'tb_bits',bits,'fixed_beta_eff_sinr_dB',sinr,'channel_diagnostics',r.channel_diagnostics, ...
    'wideband_tx_corr',NaN);
if isfield(r,'srs_propagation_tx_correlation')
    a.wideband_tx_corr=mean(r.srs_propagation_tx_correlation,'all');
end
end

function y = wrapDegrees(x)
y=mod(x+180,360)-180;
end

function [maxPhase,maxAmp,violations] = slidingDrift(phase,amp,maxLag,t)
% Changes between every pair of sampled instants separated by <= 20 ms.
maxPhase=zeros(1,size(phase,1)); maxAmp=maxPhase; violations=0;
for lag=1:min(maxLag,size(phase,2)-1)
    p=abs(wrapDegrees(phase(:,lag+1:end)-phase(:,1:end-lag)));
    a=abs(amp(:,lag+1:end)-amp(:,1:end-lag));
    maxPhase=max(maxPhase,max(p,[],2).'); maxAmp=max(maxAmp,max(a,[],2).');
    violations=violations+nnz(p>t.phase20_deg+1e-8 | a>t.amp20_dB+1e-8);
end
end

function selected = selectCandidate(screen,o,t)
base=struct('scenario','','precoder','','method','','met_screen_targets',false,'exploratory',true, ...
    'loss_dB',NaN,'gain_dB',NaN,'loss_ci95_low_dB',NaN,'gain_ci95_low_dB',NaN, ...
    'shortfall',Inf,'decision','');
eligible=repmat(base,0,1);
for i=1:numel(screen)
    c=screen(i);
    if c.scenario.is_control || ~strcmp(c.status,'complete'), continue; end
    a=aggregate(c,t);
    row=find(strcmp({a.age_group},'all') & abs([a.snr_dB]-t.snr_dB)<1e-8,1);
    if isempty(row), continue; end
    p=a(row);
    methods={'Feedback','Reselect'};
    gains=[p.feedback_mean_db p.reselect_mean_db];
    lower=[p.feedback_ci_low_db p.reselect_ci_low_db];
    ok=p.loss_mean_db>=t.loss_dB & p.loss_ci_low_db>0 & gains>=t.gain_dB & lower>0;
    [~,k]=max(gains+1000*ok);
    shortage=hypot(max(0,1-p.loss_mean_db/t.loss_dB),max(0,1-gains(k)/t.gain_dB));
    q=base; q.scenario=c.scenario.name; q.precoder=c.scenario.precoder; q.method=methods{k};
    q.met_screen_targets=ok(k); q.exploratory=~ok(k); q.loss_dB=p.loss_mean_db; q.gain_dB=gains(k);
    q.loss_ci95_low_dB=p.loss_ci_low_db; q.gain_ci95_low_dB=lower(k); q.shortfall=shortage;
    q.decision='Screen targets met';
    if ~ok(k), q.decision='Exploratory candidate; screening targets not both met'; end
    eligible(end+1,1)=q; %#ok<AGROW>
end
selected=eligible([]);
if isempty(eligible)
    assert(isempty(o.confirmNames),'UL256:Selection','No confirmed screen results to lock a method.');
    return;
end
if ~isempty(o.confirmNames)
    [found,idx]=ismember(o.confirmNames,{eligible.scenario});
    assert(all(found),'UL256:Selection','confirmNames must have completed non-control screens at -20 dB.');
    selected=eligible(idx);
    return;
end
eligible=eligible(strcmp({eligible.precoder},'nr'));
if isempty(eligible), return; end
met=[eligible.met_screen_targets];
if any(met)
    indices=find(met);
    [~,ii]=max([eligible(indices).gain_dB]);
    selected=eligible(indices(ii));
else
    keys=[[eligible.shortfall].',-[eligible.gain_dB].'];
    [~,ii]=sortrows(keys,[1,2]);
    selected=eligible(ii(1));
end
end

function [s,methods] = makeConfirmation(selected,screenScenarios)
s=selectNames(screenScenarios,{selected.scenario});
methods={selected.method};
for i=1:numel(selected)
    matched=s(i);
    matched.name=['ZERO_RF_MATCHED_' matched.name];
    matched.phase20_deg=0; matched.amp20_dB=0; matched.is_control=true;
    s(end+1)=matched; %#ok<AGROW>
    methods{end+1}=methods{i}; %#ok<AGROW>
end
end

function rows = aggregate(c,t)
gNames={'all','le20','gt20'};
rows=repmat(summaryRow(),0,1);
if c.completed_mc==0, return; end
n=c.completed_mc;
for g=1:3
    if c.entries{1}.slots_per_group(g)==0, continue; end
    for i=1:numel(c.entries{1}.snr_dB)
        x=zeros(n,5); frame=zeros(n,5); slot=frame; ph=zeros(n,4); am=ph; ph20=ph; am20=ph;
        pairPh=zeros(n,6); pairAm=pairPh; pairPh20=pairPh; pairAm20=pairPh;
        violation=zeros(n,1); pairViolation=violation; pow=zeros(n,1); oracle=pow; elapsed=pow; corr=NaN(n,1);
        for j=1:n
            a=c.entries{j};
            x(j,:)=reshape(a.mean_sinr_dB(g,i,:),1,5);
            frame(j,:)=a.frame_thp_Mbps(i,:);
            slot(j,:)=reshape(a.ul_slot_goodput_Mbps(g,i,:),1,5);
            ph(j,:)=a.srs_phase_abs_max_deg(g,:); am(j,:)=a.srs_amp_abs_max_dB(g,:);
            ph20(j,:)=a.phase20_abs_max_deg; am20(j,:)=a.amp20_abs_max_dB;
            pairPh(j,:)=a.pair_srs_phase_abs_max_deg(g,:); pairAm(j,:)=a.pair_srs_amp_abs_max_dB(g,:);
            pairPh20(j,:)=a.pair_phase20_abs_max_deg; pairAm20(j,:)=a.pair_amp20_abs_max_dB;
            pairViolation(j)=a.pair_rf20_violation_count;
            violation(j)=a.rf20_violation_count; pow(j)=a.tx_power_error;
            oracle(j)=a.oracle_relative_energy_error; elapsed(j)=a.elapsed_s; corr(j)=a.wideband_tx_corr;
        end
        y=summaryRow(); y.stage=c.stage; y.scenario=c.scenario.name; y.status=c.status;
        y.age_group=gNames{g}; y.snr_dB=c.entries{1}.snr_dB(i); y.completed_mc=n;
        y.precoder=c.scenario.precoder;
        y.slots_per_mc=c.entries{1}.slots_per_group(g); y.rank=c.scenario.rank;
        y.locked_method=c.locked_method; y.is_control=c.scenario.is_control;
        [y.loss_mean_db,y.loss_ci_low_db,y.loss_ci_high_db]=pairedStats(x(:,1)-x(:,2));
        [y.feedback_mean_db,y.feedback_ci_low_db,y.feedback_ci_high_db]=pairedStats(x(:,3)-x(:,2));
        [y.reselect_mean_db,y.reselect_ci_low_db,y.reselect_ci_high_db]=pairedStats(x(:,4)-x(:,2));
        [y.oracle_ideal_mean_db,y.oracle_ideal_ci_low_db,y.oracle_ideal_ci_high_db]=pairedStats(x(:,5)-x(:,1));
        if g==1
            y.ideal_frame_thp_mbps=mean(frame(:,1)); y.uncal_frame_thp_mbps=mean(frame(:,2));
            y.feedback_frame_thp_mbps=mean(frame(:,3)); y.reselect_frame_thp_mbps=mean(frame(:,4));
            y.oracle_frame_thp_mbps=mean(frame(:,5));
            [y.feedback_frame_thp_gain_mbps,y.feedback_frame_thp_ci_low_mbps,y.feedback_frame_thp_ci_high_mbps]=pairedStats(frame(:,3)-frame(:,2));
            [y.reselect_frame_thp_gain_mbps,y.reselect_frame_thp_ci_low_mbps,y.reselect_frame_thp_ci_high_mbps]=pairedStats(frame(:,4)-frame(:,2));
        end
        y.feedback_ul_slot_goodput_gain_mbps=mean(slot(:,3)-slot(:,2));
        y.reselect_ul_slot_goodput_gain_mbps=mean(slot(:,4)-slot(:,2));
        y.phase20_abs_max_deg=max(ph20(:)); y.amp20_abs_max_db=max(am20(:));
        y.srs_phase_abs_max_deg=max(ph(:)); y.srs_amp_abs_max_db=max(am(:));
        y.pair_phase20_abs_max_deg=max(pairPh20(:)); y.pair_amp20_abs_max_db=max(pairAm20(:));
        y.pair_srs_phase_abs_max_deg=max(pairPh(:)); y.pair_srs_amp_abs_max_db=max(pairAm(:));
        y.pair_rf20_violation_count=sum(pairViolation);
        y.rf20_violation_count=sum(violation); y.tx_power_max_error=max(pow);
        y.oracle_max_relative_energy_error=max(oracle); y.mean_run_seconds=mean(elapsed);
        y.wideband_tx_corr=mean(corr,'omitnan');
        y.is_predeclared_endpoint=(g==1 && abs(y.snr_dB-t.snr_dB)<1e-8);
        if y.is_predeclared_endpoint && strcmp(y.status,'complete') && ~y.is_control
            idx=[1,2];
            if strcmp(c.locked_method,'Feedback'), idx=1; end
            if strcmp(c.locked_method,'Reselect'), idx=2; end
            means=[y.feedback_mean_db y.reselect_mean_db];
            lows=[y.feedback_ci_low_db y.reselect_ci_low_db];
            y.met_targets=y.loss_mean_db>=t.loss_dB && y.loss_ci_low_db>0 && ...
                any(means(idx)>=t.gain_dB & lows(idx)>0);
        end
        rows(end+1,1)=y; %#ok<AGROW>
    end
end
end

function r = summaryRow()
r=struct('stage','','scenario','','status','','precoder','','age_group','','snr_dB',NaN,'rank',0, ...
    'completed_mc',0,'slots_per_mc',0,'is_control',false,'locked_method','', ...
    'loss_mean_db',NaN,'loss_ci_low_db',NaN,'loss_ci_high_db',NaN, ...
    'feedback_mean_db',NaN,'feedback_ci_low_db',NaN,'feedback_ci_high_db',NaN, ...
    'reselect_mean_db',NaN,'reselect_ci_low_db',NaN,'reselect_ci_high_db',NaN, ...
    'oracle_ideal_mean_db',NaN,'oracle_ideal_ci_low_db',NaN,'oracle_ideal_ci_high_db',NaN, ...
    'ideal_frame_thp_mbps',NaN,'uncal_frame_thp_mbps',NaN,'feedback_frame_thp_mbps',NaN, ...
    'reselect_frame_thp_mbps',NaN,'oracle_frame_thp_mbps',NaN, ...
    'feedback_frame_thp_gain_mbps',NaN,'feedback_frame_thp_ci_low_mbps',NaN,'feedback_frame_thp_ci_high_mbps',NaN, ...
    'reselect_frame_thp_gain_mbps',NaN,'reselect_frame_thp_ci_low_mbps',NaN,'reselect_frame_thp_ci_high_mbps',NaN, ...
    'feedback_ul_slot_goodput_gain_mbps',NaN,'reselect_ul_slot_goodput_gain_mbps',NaN, ...
    'phase20_abs_max_deg',NaN,'amp20_abs_max_db',NaN,'srs_phase_abs_max_deg',NaN,'srs_amp_abs_max_db',NaN, ...
    'pair_phase20_abs_max_deg',NaN,'pair_amp20_abs_max_db',NaN,'pair_srs_phase_abs_max_deg',NaN, ...
    'pair_srs_amp_abs_max_db',NaN,'pair_rf20_violation_count',NaN, ...
    'rf20_violation_count',NaN,'tx_power_max_error',NaN,'oracle_max_relative_energy_error',NaN, ...
    'mean_run_seconds',NaN,'wideband_tx_corr',NaN,'is_predeclared_endpoint',false,'met_targets',false);
end

function [mu,lo,hi] = pairedStats(x)
n=numel(x); mu=mean(x); lo=NaN; hi=NaN;
if n<2, return; end
v=n-1; q=sqrt(v*(1/betaincinv(.05,v/2,.5)-1));
half=q*std(x,0)/sqrt(n);
lo=mu-half; hi=mu+half;
end

function exportStatus(stage,o,fid)
status=struct('scenario',{},'stage',{},'status',{},'completed_mc',{}, ...
    'requested_mc',{},'current_mc',{},'locked_method',{},'checkpoint_file',{},'error_message',{});
for i=1:numel(stage)
    c=stage(i);
    status(i)=struct('scenario',c.scenario.name,'stage',c.stage,'status',c.status, ...
        'completed_mc',c.completed_mc,'requested_mc',c.requested_mc,'current_mc',c.current_mc, ...
        'locked_method',c.locked_method,'checkpoint_file',c.checkpoint_file,'error_message',c.error_message);
end
writeTableSafe(struct2table(status(:),'AsArray',true),fullfile(o.outputDir,[stage(1).stage '_status.csv']),fid);
end

function exportStage(stage,o,t,fid)
exportStatus(stage,o,fid);
rows=repmat(summaryRow(),0,1);
mc=cell(0,16); rf=cell(0,10); pairRF=cell(0,12); groups={'all','le20','gt20'};
for i=1:numel(stage)
    c=stage(i);
    rows=[rows;aggregate(c,t)]; %#ok<AGROW>
    for j=1:c.completed_mc
        a=c.entries{j};
        for g=1:3
            if a.slots_per_group(g)==0, continue; end
            for p=1:4
                rf(end+1,:)={c.stage,c.scenario.name,j,a.seed,groups{g},p, ...
                    a.phase20_abs_max_deg(p),a.amp20_abs_max_dB(p), ...
                    a.srs_phase_abs_max_deg(g,p),a.srs_amp_abs_max_dB(g,p)}; %#ok<AGROW>
            end
            for p=1:6
                pairRF(end+1,:)={c.stage,c.scenario.name,j,a.seed,groups{g},a.port_pairs(p,1),a.port_pairs(p,2), ...
                    a.pair_phase20_abs_max_deg(p),a.pair_amp20_abs_max_dB(p), ...
                    a.pair_srs_phase_abs_max_deg(g,p),a.pair_srs_amp_abs_max_dB(g,p), ...
                    a.pair_rf20_violation_count}; %#ok<AGROW>
            end
            for k=1:numel(a.snr_dB)
                x=reshape(a.mean_sinr_dB(g,k,:),1,5);
                thp=NaN(1,5);
                if g==1, thp=a.frame_thp_Mbps(k,:); end
                mc(end+1,:)={c.stage,c.scenario.name,c.scenario.precoder,j,a.seed,groups{g},a.snr_dB(k), ...
                    x(1)-x(2),x(3)-x(2),x(4)-x(2),x(5)-x(1), ...
                    thp(2),thp(3),thp(4),a.wideband_tx_corr,a.elapsed_s}; %#ok<AGROW>
            end
        end
    end
end
prefix=stage(1).stage;
writeTableSafe(struct2table(rows,'AsArray',true),fullfile(o.outputDir,[prefix '_summary.csv']),fid);
tab=cell2table(mc,'VariableNames',{'stage','scenario','precoder','mc_index','seed','age_group','snr_dB', ...
    'loss_db','feedback_db','reselect_db','oracle_ideal_db','uncal_frame_thp_mbps', ...
    'feedback_frame_thp_mbps','reselect_frame_thp_mbps','wideband_tx_corr','elapsed_s'});
writeTableSafe(tab,fullfile(o.outputDir,[prefix '_mc.csv']),fid);
tab=cell2table(rf,'VariableNames',{'stage','scenario','mc_index','seed','age_group','port', ...
    'phase20_abs_max_deg','amp20_abs_max_db','srs_phase_abs_max_deg','srs_amp_abs_max_db'});
writeTableSafe(tab,fullfile(o.outputDir,[prefix '_rf.csv']),fid);
tab=cell2table(pairRF,'VariableNames',{'stage','scenario','mc_index','seed','age_group','port_i','port_j', ...
    'pair_phase20_abs_max_deg','pair_amp20_abs_max_db','pair_srs_phase_abs_max_deg', ...
    'pair_srs_amp_abs_max_db','pair_rf20_violation_count_all_pairs'});
writeTableSafe(tab,fullfile(o.outputDir,[prefix '_pair_rf.csv']),fid);
end

function saveAtomic(path,name,value)
tmp=[tempname(fileparts(path)) '.mat'];
payload=struct(); payload.(name)=value;
save(tmp,'-struct','payload','-v7');
movefile(tmp,path,'f');
end

function writeTableSafe(t,path,fid)
tmp=[tempname(fileparts(path)) '.csv'];
cleanup=onCleanup(@()deleteTemporary(tmp)); %#ok<NASGU>
try
    writetable(t,tmp,'Encoding','UTF-8');
    movefile(tmp,path,'f');
catch ME
    logLine(fid,'EXPORT WARNING %s: %s\n',path,ME.message);
end
end

function deleteTemporary(p)
if isfile(p), delete(p); end
end

function ts=timestamp()
ts=char(datetime('now','TimeZone','Asia/Shanghai','Format','yyyy-MM-dd HH:mm:ss'));
end

function logLine(fid,varargin)
fprintf(1,varargin{:});
fprintf(fid,varargin{:});
end
