function report = runNrStrict256(opts)
%RUNNRSTRICT256 Screen NR at actual SRS age <=20 ms, then independently confirm.
% opts.stage is 'screen', 'confirm', or 'both'. opts.dryRun writes nothing.
% Reuse the same outputDir and screen options for confirmation and resume.
% Base MATLAB is the default. Set opts.backend='5g-toolbox' for a reference run.
% Set opts.channelProfile to one of 'CDL-A' through 'CDL-E'.
if nargin==0, opts=struct(); end
root=fileparts(mfilename('fullpath'));
o=options(opts,root);
scenario=nrStrictScenarios();
if ~isempty(o.scenarioNames)
    [found,index]=ismember(o.scenarioNames,{scenario.name});
    assert(all(found),'ULNrStrict:Scenario','Unknown scenario name.');
    scenario=scenario(index);
end
scenario=reshape(scenario,1,[]);
assert(~isempty(scenario) && any([scenario.selectable]),'ULNrStrict:Scenario','Select at least one selectable scenario.');
base=runCoherence256Search(struct('dryRun',true, ...
    'scenarioNames',{{'PA_LOAD_COMPACT_R1'}},'channelProfile',o.channelProfile));
template=base.plan.screen_configs{1};
manifest=base.source_manifest;
for name={'runNrStrict256.m','nrStrictScenarios.m'}
    manifest(end+1)=struct('name',name{1},'sha256',fileHash(fullfile(root,name{1}))); %#ok<AGROW>
end
source=ul4tRunnerCacheKey(manifest,'nr-strict-256-production-sha256');
threshold=struct('snr_dB',-20,'loss_dB',.5,'gain_dB',.4, ...
    'phase_deg',40,'amplitude_dB',4,'age_ms',20,'duration_ms',50);
cfg=cell(1,numel(scenario));
for i=1:numel(scenario)
    cfg{i}=makeConfig(template,scenario(i),o.screenMC,o.screenSeed,o.screenSNRdB,o,threshold);
end
plan=struct('scenarios',scenario,'screen_configs',{cfg},'threshold',threshold, ...
    'tx_rx',[4 256],'allocated_RB',288,'center_frequency_Hz',6.7e9,'cdl_backend',o.backend, ...
    'channel_profile',o.channelProfile, ...
    'case_order',{{'Ideal','Uncal','Feedback','Reselect','Oracle'}}, ...
    'selection_rule',['Maximize min(loss/0.50, Feedback/0.40, Reselect/0.40) at -20 dB. ' ...
        'Only selectable NR scenarios enter selection. Meeting targets requires all three paired CI95 lower bounds above zero.'], ...
    'throughput_definition',['Good bits summed over the actual SRS-age subset and divided by the full 50 ms timeline. ' ...
        'Main throughput does not deduct calibration overhead.'], ...
    'rf_definition',['Report same-port and six port-pair changes. ' ...
        'Enforce the 40 degree and 4 dB global pairwise spans over the realized trajectories.']);
report=struct('schema','nr_strict_256_v1','status','planned','options',o, ...
    'source_sha256',source,'source_manifest',manifest,'plan',plan);
if o.dryRun, return; end

if ~isfolder(o.outputDir), mkdir(o.outputDir); end
fid=fopen(fullfile(o.outputDir,'run.log'),'a','n','UTF-8');
assert(fid>0,'ULNrStrict:Log','Cannot open campaign log.');
closeLog=onCleanup(@()fclose(fid));
logLine(fid,'START %s stage=%s source=%s\n',stamp(),o.stage,source);
logLine(fid,'Channel backend: %s.\n',o.backend);
logLine(fid,'NR only, 4T256R, 288 RB research grid, 6.7 GHz, %s, 50 ms, actual SRS age <=20 ms.\n',o.channelProfile);
logLine(fid,'Selection at -20 dB requires loss>=0.50 dB and both gains>=0.40 dB for the declared targets.\n');
logLine(fid,'Score: min(loss/0.50, Feedback/0.40, Reselect/0.40). CI unit is paired MC.\n');
safeTable(manifestTable(scenario),fullfile(o.outputDir,'scenario_manifest.csv'),fid);
report.status='running';
atomicSave(fullfile(o.outputDir,'plan.mat'),'plan',report);
screen=loadStage('screen',scenario,cfg,source,o,fid);
if strcmp(o.stage,'confirm')
    assert(all(ismember({screen.status},{'complete','failed'})) && any(strcmp({screen.status},'complete')), ...
        'ULNrStrict:ScreenRequired','Confirm requires completed matching screens and a saved selection.');
else
    screen=executeStage(screen,o,threshold,fid);
end
report.screen=screen;
screenKey=ul4tRunnerCacheKey({screen.key},source);
lockPath=fullfile(o.outputDir,['selection_' screenKey(1:12) '.mat']);
if isfile(lockPath) && (o.resume || strcmp(o.stage,'confirm'))
    locked=load(lockPath,'selection'); selection=locked.selection;
    assert(strcmp(selection.screen_key,screenKey) && any(strcmp(selection.scenario,{screen.scenario_name})), ...
        'ULNrStrict:Selection','Saved selection does not match the screens.');
elseif strcmp(o.stage,'confirm')
    error('ULNrStrict:Selection','Missing saved selection. Run the screen stage first.');
else
    selection=chooseCandidate(screen,threshold);
    if ~isempty(selection)
        selection.screen_key=screenKey;
        atomicSave(lockPath,'selection',selection);
    end
end
report.selection=selection;
if ~isempty(selection)
    safeTable(struct2table(selection,'AsArray',true),fullfile(o.outputDir,'selection.csv'),fid);
    logLine(fid,'LOCKED %s score=%.4f targets=%d.\n',selection.scenario,selection.score,selection.met_targets);
else
    logLine(fid,'No selectable scenario completed.\n');
end

report.confirm=[];
if ~strcmp(o.stage,'screen') && ~isempty(selection)
    chosen=screen(strcmp({screen.scenario_name},selection.scenario));
    assert(isscalar(chosen) && strcmp(chosen.status,'complete') && chosen.scenario.selectable, ...
        'ULNrStrict:Selection','Locked scenario must be a completed selectable candidate.');
    candidate=chosen.scenario;
    zero=candidate;
    zero.name=['ZERO_RF_MATCHED_' zero.name];
    zero.phase20_deg=0; zero.amp20_dB=0; zero.selectable=false;
    confirmCfg={makeConfig(template,candidate,o.finalMC,o.finalSeed,o.finalSNRdB,o,threshold), ...
        makeConfig(template,zero,o.finalMC,o.finalSeed,-20,o,threshold)};
    confirm=loadStage('confirm',[candidate zero],confirmCfg,source,o,fid);
    safeTable(manifestTable([candidate zero]),fullfile(o.outputDir,'confirmation_manifest.csv'),fid);
    report.confirm=executeStage(confirm,o,threshold,fid);
end
report.status='complete';
allStages=screen;
if ~isempty(report.confirm), allStages=[allStages,report.confirm]; end
if any(strcmp({allStages.status},'failed')) || (~strcmp(o.stage,'screen') && isempty(report.confirm))
    report.status='complete_with_failures';
end
report.completed_at=stamp();
atomicSave(fullfile(o.outputDir,'report.mat'),'report',report);
logLine(fid,'Campaign %s at %s.\n',upper(report.status),report.completed_at);
end

function o=options(in,root)
assert(isstruct(in)&&isscalar(in),'ULNrStrict:Options','Expected a scalar options struct.');
o=struct('stage','both','dryRun',false,'resume',true,'screenMC',8,'finalMC',24, ...
    'screenSeed',204611181,'finalSeed',205611181,'screenSNRdB',-20,'finalSNRdB',-24:2:-12, ...
    'scenarioNames',{{}},'backend','base-cdl','channelProfile','CDL-C','outputDir','');
names=fieldnames(in);
assert(all(ismember(names,fieldnames(o))),'ULNrStrict:Options','Unknown options field.');
for i=1:numel(names), o.(names{i})=in.(names{i}); end
o.stage=validatestring(o.stage,{'both','screen','confirm'});
o.backend=validatestring(o.backend,{'base-cdl','base-cdl-c','5g-toolbox'});
if strcmp(o.backend,'base-cdl-c'), o.backend='base-cdl'; end
o.channelProfile=validatestring(o.channelProfile,{'CDL-A','CDL-B','CDL-C','CDL-D','CDL-E'});
if isempty(o.outputDir)
    o.outputDir=fullfile(root,['results_nr_strict256_' o.backend '_' lower(o.channelProfile)]);
end
for name={'dryRun','resume'}, validateattributes(o.(name{1}),{'logical','numeric'},{'scalar','binary'}); end
for name={'screenMC','finalMC','screenSeed','finalSeed'}
    validateattributes(o.(name{1}),{'numeric'},{'scalar','integer','finite','positive'});
end
assert(o.finalSeed>o.screenSeed+o.screenMC-1 && o.finalSeed+o.finalMC-1<2^32, ...
    'ULNrStrict:Options','Screen and final MC seeds must be distinct and in range.');
for name={'screenSNRdB','finalSNRdB'}
    validateattributes(o.(name{1}),{'numeric'},{'vector','real','finite','nonempty'});
    x=sort(reshape(double(o.(name{1})),1,[]));
    assert(numel(unique(x))==numel(x) && any(x==-20),'ULNrStrict:Options','Unique SNR points must include -20 dB.');
    o.(name{1})=x;
end
o.scenarioNames=cellstr(string(o.scenarioNames));
assert(numel(unique(o.scenarioNames))==numel(o.scenarioNames),'ULNrStrict:Options','Duplicate scenarios.');
o.outputDir=char(java.io.File(char(o.outputDir)).getCanonicalPath());
end

function cfg=makeConfig(cfg,s,nMC,seed,snr,o,t)
cfg.MC=nMC; cfg.seed=seed; cfg.output_dir=o.outputDir;
cfg.cdl.backend=o.backend;
if strcmp(o.backend,'5g-toolbox'), cfg.ul.dmrs_segments_rb=[144 144]; end
cfg.ul.duration_ms=t.duration_ms;
cfg.ul.rx_noise_power=10.^(-snr/10);
cfg.ul.rank=s.rank; cfg.ul.ue_speed_kmh=s.speed_kmh;
if isnan(s.spacing_lambda)
    cfg.ul.ue_layout='handheld';
else
    cfg.ul.ue_layout='upa';
    cfg.ul.tx_positions_lambda=s.spacing_lambda*[0 0 0;0 0 1;1 0 0;1 0 1];
end
p=cfg.ul.coherence_drift.bounded;
p.tau_ms=s.tau_ms; p.phase20_deg=s.phase20_deg; p.amp20_dB=s.amp20_dB;
if ~isempty(s.direction), p.direction=s.direction; end
if isfinite(s.initial_offset_deg)
    p.initial_offset_deg=s.initial_offset_deg;
    p.initial_jitter_deg=s.initial_jitter_deg;
end
cfg.ul.coherence_drift.bounded=p;
for i=1:numel(cfg.ul.cases), cfg.ul.cases(i).sets=s.sets; end
end

function h=fileHash(path)
fid=fopen(path,'rb');
assert(fid>0,'ULNrStrict:Source','Cannot open %s for hashing.',path);
cleanup=onCleanup(@()fclose(fid));
bytes=fread(fid,Inf,'*uint8');
digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(bytes);
h=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end

function m=manifestTable(s)
% CSV fields are scalar even when a preset has a vector of port directions.
m=table(string({s.name}).',[s.spacing_lambda].',[s.speed_kmh].',[s.rank].', ...
    [s.tau_ms].',[s.phase20_deg].',[s.amp20_dB].',[s.selectable].', ...
    'VariableNames',{'scenario','spacing_lambda','speed_kmh','rank','tau_ms','phase_bound_deg','amp_bound_db','selectable'});
m.direction=string(arrayfun(@(x)mat2str(x.direction),s,'UniformOutput',false)).';
m.initial_offset_deg=string(arrayfun(@(x)mat2str(x.initial_offset_deg),s,'UniformOutput',false)).';
m.initial_jitter_deg=string(arrayfun(@(x)mat2str(x.initial_jitter_deg),s,'UniformOutput',false)).';
m.sets=string(arrayfun(@(x)strjoin(x.sets,','),s,'UniformOutput',false)).';
end

function cp=loadStage(name,scenario,cfg,source,o,fid)
base=struct('schema','nr_strict_256_mc_v1','scenario_name','','scenario',struct(), ...
    'cfg',struct(),'key','','source_sha256',source,'stage',name,'checkpoint_file','', ...
    'status','pending','completed_mc',0,'requested_mc',0,'current_mc',0,'entries',{{}}, ...
    'error_identifier','','error_message','');
cp=repmat(base,1,numel(scenario));
for i=1:numel(scenario)
    e=base; e.scenario=scenario(i); e.scenario_name=scenario(i).name; e.cfg=cfg{i};
    e.requested_mc=e.cfg.MC;
    e.key=ul4tRunnerCacheKey(struct('cfg',e.cfg,'scenario',e.scenario,'stage',name),source);
    e.checkpoint_file=fullfile(o.outputDir,sprintf('%s_%s_%s.mat',name,e.scenario_name,e.key(1:12)));
    cp(i)=e;
    if ~(o.resume || strcmp(o.stage,'confirm')) || ~isfile(e.checkpoint_file), continue; end
    try
        old=load(e.checkpoint_file,'checkpoint'); c=old.checkpoint;
        assert(validCheckpoint(c,e),'ULNrStrict:Checkpoint','Configuration or MC data mismatch.');
        cp(i)=c;
        logLine(fid,'RESUME %s %s %d/%d.\n',name,c.scenario_name,c.completed_mc,c.requested_mc);
    catch ME
        logLine(fid,'Ignoring %s because %s\n',e.checkpoint_file,ME.message);
    end
end
end

function valid=validCheckpoint(c,e)
valid=isequal(sort(fieldnames(c)),sort(fieldnames(e))) && strcmp(c.schema,e.schema) && ...
    strcmp(c.key,e.key) && strcmp(c.source_sha256,e.source_sha256) && ...
    strcmp(c.stage,e.stage) && strcmp(c.scenario_name,e.scenario_name) && ...
    isequaln(c.cfg,e.cfg) && isequaln(c.scenario,e.scenario) && ...
    strcmp(c.checkpoint_file,e.checkpoint_file) && c.requested_mc==e.requested_mc && ...
    ismember(c.completed_mc,0:e.requested_mc) && numel(c.entries)==c.completed_mc && ...
    any(strcmp(c.status,{'pending','running','complete','failed'}));
if ~valid, return; end
if strcmp(c.status,'complete') && c.completed_mc~=c.requested_mc, valid=false; return; end
snr=10*log10(1./c.cfg.ul.rx_noise_power);
for i=1:c.completed_mc
    a=c.entries{i};
    valid=a.mc_index==i && a.seed==c.cfg.seed+i-1 && ...
        isequal(a.snr_dB,snr) && isequal(size(a.thp_Mbps),[numel(snr),5]) && ...
        isequal(size(a.mean_sinr_dB),[numel(snr),5]) && all(isfinite(a.mean_sinr_dB),'all') && ...
        all(isfinite(a.thp_Mbps),'all') && size(a.slot_bits,1)==a.selected_slots && size(a.slot_bits,3)==5;
    if ~valid, return; end
end
end

function cp=executeStage(cp,o,t,fid)
exportStage(cp,o,t,fid);
for i=1:numel(cp)
    c=cp(i);
    if strcmp(c.status,'complete'), continue; end
    c.error_identifier=''; c.error_message='';
    try
        for j=c.completed_mc+1:c.requested_mc
            c.current_mc=j; c.status='running';
            atomicSave(c.checkpoint_file,'checkpoint',c); cp(i)=c;
            exportStatus(cp,o,fid);
            cfg=c.cfg; cfg.MC=1; cfg.seed=c.cfg.seed+j-1;
            logLine(fid,'%s %s MC %d/%d START seed=%d.\n',upper(c.stage),c.scenario_name,j,c.requested_mc,cfg.seed);
            timer=tic;
            result=simulateUL4T(cfg);
            a=measureMC(result, j,cfg.seed,t);
            a.elapsed_s=toc(timer);
            c.entries{j}=a; c.completed_mc=j;
            atomicSave(c.checkpoint_file,'checkpoint',c); cp(i)=c;
            exportStage(cp,o,t,fid);
            logLine(fid,'%s %s MC %d/%d DONE %.2fs.\n',upper(c.stage),c.scenario_name,j,c.requested_mc,a.elapsed_s);
        end
        c.status='complete'; c.current_mc=0;
        atomicSave(c.checkpoint_file,'checkpoint',c);
    catch ME
        c.status='failed'; c.error_identifier=ME.identifier; c.error_message=ME.message;
        atomicSave(c.checkpoint_file,'checkpoint',c);
        logLine(fid,'FAILED %s %s MC %d: %s\n',c.stage,c.scenario_name,c.current_mc, ...
            getReport(ME,'extended','hyperlinks','off'));
    end
    cp(i)=c;
    exportStage(cp,o,t,fid);
end
end

function a=measureMC(r,idx,seed,t)
assert(r.fixed_beta_linear==1 && size(r.fixed_beta_eff_sinr_dB,4)==5, ...
    'ULNrStrict:Cases','Expected five cases scored with fixed beta 1.');
assert(r.cfg.N==256 && r.cfg.ofdm.num_tones==3456 && r.cfg.ul.num_tx==4 && ...
    r.cfg.ul.calibration.num_re==128 && r.cfg.ul.duration_ms==t.duration_ms, ...
    'ULNrStrict:Grid','Unexpected link configuration.');
ns=numel(r.rx_noise_power); total=numel(r.schedule.ul_ms);
selected=r.schedule.srs_age_ms(:)<=t.age_ms;
assert(any(selected),'ULNrStrict:Schedule','No PUSCH with SRS age up to 20 ms.');
sinr=reshape(r.fixed_beta_eff_sinr_dB,total,ns,5);
bits=reshape(r.good_bits,total,ns,5);
assert(all(isfinite(sinr),'all') && all(isfinite(bits),'all'),'ULNrStrict:Finite','Nonfinite link metric.');
shortSinr=sinr(selected,:,:); shortBits=bits(selected,:,:);
d=r.coherence_drift;
phase=reshape(d.phase_deg,4,[]); amp=reshape(d.amplitude_dB,4,[]);
pairs=nchoosek(1:4,2);
pairPhase=phase(pairs(:,1),:)-phase(pairs(:,2),:);
pairAmp=amp(pairs(:,1),:)-amp(pairs(:,2),:);
phaseSpan=max(phase,[],2)-min(phase,[],2); ampSpan=max(amp,[],2)-min(amp,[],2);
pairPhaseSpan=max(pairPhase,[],2)-min(pairPhase,[],2);
pairAmpSpan=max(pairAmp,[],2)-min(pairAmp,[],2);
checkLimits([phaseSpan;pairPhaseSpan],[ampSpan;pairAmpSpan],t,'global trajectory');
pu=d.pusch_sample_index; srs=d.srs_sample_index(r.schedule.srs_id);
srsPhase=abs(wrap(phase(:,pu)-phase(:,srs))); srsAmp=abs(amp(:,pu)-amp(:,srs));
pairSrsPhase=abs(wrap(pairPhase(:,pu)-pairPhase(:,srs))); pairSrsAmp=abs(pairAmp(:,pu)-pairAmp(:,srs));
checkLimits([srsPhase;pairSrsPhase],[srsAmp;pairSrsAmp],t,'all actual SRS ages');
powerError=max(abs([r.pusch_tx_power(:);r.calibration.tx_power(:)]-1));
energy=reshape(r.pusch_effective_energy,total,5);
oracleError=max(abs(energy(:,5)-energy(:,1))./max(energy(:,1),realmin));
assert(powerError<1e-10 && oracleError<1e-9,'ULNrStrict:Invariant','TX power or Oracle energy mismatch.');
if r.cfg.ul.coherence_drift.bounded.phase20_deg==0 && r.cfg.ul.coherence_drift.bounded.amp20_dB==0
    assert(all(abs(sinr(:,:,1)-sinr(:,:,2))<1e-10,'all') && all(bits(:,:,1)==bits(:,:,2),'all'), ...
        'ULNrStrict:ZeroRF','Ideal and Uncal differ for zero drift.');
end
baseline=reshape(r.precoder_tpmi(1,r.schedule.srs_id,2),total,1);
reselect=reshape(r.calibration.tpmi(1,:,:,2),total,ns);
switchRate=mean(reselect(selected,:)~=baseline(selected),1);
a=struct('mc_index',idx,'seed',seed,'snr_dB',reshape(10*log10(1./r.rx_noise_power),1,[]), ...
    'mean_sinr_dB',reshape(mean(shortSinr,1),ns,5), ...
    'thp_Mbps',reshape(sum(shortBits,1),ns,5)/(1000*t.duration_ms), ...
    'selected_slots',nnz(selected),'duty',nnz(selected)*r.schedule.slot_ms/t.duration_ms, ...
    'selected_ul_ms',r.schedule.ul_ms(selected),'selected_srs_age_ms',r.schedule.srs_age_ms(selected), ...
    'srs_ms',r.schedule.srs_ms,'slot_bits',shortBits,'slot_sinr_dB',shortSinr, ...
    'port_phase_span_deg',phaseSpan,'port_amp_span_dB',ampSpan,'pair_ports',pairs, ...
    'pair_phase_span_deg',pairPhaseSpan,'pair_amp_span_dB',pairAmpSpan, ...
    'pair_srs_phase_abs_deg',pairSrsPhase,'pair_srs_amp_abs_dB',pairSrsAmp, ...
    'port_srs_phase_abs_deg',srsPhase,'port_srs_amp_abs_dB',srsAmp, ...
    'oracle_relative_energy_error',oracleError,'tx_power_error',powerError, ...
    'reselect_tpmi_switch_rate',switchRate,'wideband_tx_corr',mean(r.srs_propagation_tx_correlation,'all'));
end

function checkLimits(phase,amp,t,label)
assert(all(phase<=t.phase_deg+1e-8,'all') && all(amp<=t.amplitude_dB+1e-8,'all'), ...
    'ULNrStrict:RFEnvelope','Phase or amplitude envelope exceeded in %s.',label);
end

function y=wrap(x)
y=mod(x+180,360)-180;
end

function selection=chooseCandidate(screen,t)
selection=[];
for i=1:numel(screen)
    c=screen(i);
    if ~c.scenario.selectable || ~strcmp(c.status,'complete'), continue; end
    s=summarize(c,t);
    r=s(abs([s.snr_dB]-t.snr_dB)<1e-8);
    if isempty(selection) || r.met_targets>selection.met_targets || ...
            (r.met_targets==selection.met_targets && r.score>selection.score)
        selection=struct('scenario',r.scenario,'score',r.score,'met_targets',r.met_targets, ...
            'exploratory',~r.met_targets,'loss_db',r.loss_db,'feedback_db',r.feedback_db, ...
            'reselect_db',r.reselect_db,'loss_ci_low_db',r.loss_ci_low_db, ...
            'feedback_ci_low_db',r.feedback_ci_low_db,'reselect_ci_low_db',r.reselect_ci_low_db);
    end
end
end

function rows=summarize(c,t)
rows=repmat(rowTemplate(),0,1);
if c.completed_mc==0, return; end
n=c.completed_mc; a=c.entries{1};
for k=1:numel(a.snr_dB)
    sinr=zeros(n,5); thp=sinr; power=zeros(n,1); oracle=power; corr=power; switchRate=power;
    phaseSpan=power; ampSpan=power; pairPhase=power; pairAmp=power;
    for m=1:n
        e=c.entries{m};
        assert(e.selected_slots==a.selected_slots && e.duty==a.duty,'ULNrStrict:Schedule','Schedule differs between MCs.');
        sinr(m,:)=e.mean_sinr_dB(k,:); thp(m,:)=e.thp_Mbps(k,:);
        power(m)=e.tx_power_error; oracle(m)=e.oracle_relative_energy_error;
        corr(m)=e.wideband_tx_corr; switchRate(m)=e.reselect_tpmi_switch_rate(k);
        phaseSpan(m)=max(e.port_phase_span_deg); ampSpan(m)=max(e.port_amp_span_dB);
        pairPhase(m)=max(e.pair_phase_span_deg); pairAmp(m)=max(e.pair_amp_span_dB);
    end
    r=rowTemplate(); r.scenario=c.scenario_name; r.stage=c.stage; r.status=c.status;
    r.snr_dB=a.snr_dB(k); r.mc=n; r.selected_slots_per_mc=a.selected_slots; r.duty=a.duty;
    r.selectable=c.scenario.selectable;
    [r.loss_db,r.loss_ci_low_db,r.loss_ci_high_db]=pairedCI(sinr(:,1)-sinr(:,2));
    [r.feedback_db,r.feedback_ci_low_db,r.feedback_ci_high_db]=pairedCI(sinr(:,3)-sinr(:,2));
    [r.reselect_db,r.reselect_ci_low_db,r.reselect_ci_high_db]=pairedCI(sinr(:,4)-sinr(:,2));
    [r.oracle_ideal_db,r.oracle_ideal_ci_low_db,r.oracle_ideal_ci_high_db]=pairedCI(sinr(:,5)-sinr(:,1));
    avg=mean(thp,1);
    r.ideal_thp_mbps=avg(1); r.uncal_thp_mbps=avg(2); r.feedback_thp_mbps=avg(3);
    r.reselect_thp_mbps=avg(4); r.oracle_thp_mbps=avg(5);
    [r.feedback_thp_gain_mbps,r.feedback_thp_ci_low_mbps,r.feedback_thp_ci_high_mbps]=pairedCI(thp(:,3)-thp(:,2));
    [r.reselect_thp_gain_mbps,r.reselect_thp_ci_low_mbps,r.reselect_thp_ci_high_mbps]=pairedCI(thp(:,4)-thp(:,2));
    r.wideband_tx_corr=mean(corr); r.reselect_tpmi_switch_rate=mean(switchRate);
    r.port_phase_span_deg=max(phaseSpan); r.port_amp_span_db=max(ampSpan);
    r.pair_phase_span_deg=max(pairPhase); r.pair_amp_span_db=max(pairAmp);
    r.oracle_max_relative_energy_error=max(oracle); r.tx_power_max_error=max(power);
    r.score=min([r.loss_db/t.loss_dB,r.feedback_db/t.gain_dB,r.reselect_db/t.gain_dB]);
    r.met_targets=strcmp(c.status,'complete') && r.score>=1 && ...
        all([r.loss_ci_low_db,r.feedback_ci_low_db,r.reselect_ci_low_db]>0);
    rows(end+1,1)=r; %#ok<AGROW>
end
end

function r=rowTemplate()
r=struct('scenario','','stage','','status','','snr_dB',NaN,'mc',0,'selected_slots_per_mc',0,'duty',NaN, ...
    'selectable',false,'score',NaN,'met_targets',false, ...
    'loss_db',NaN,'loss_ci_low_db',NaN,'loss_ci_high_db',NaN, ...
    'feedback_db',NaN,'feedback_ci_low_db',NaN,'feedback_ci_high_db',NaN, ...
    'reselect_db',NaN,'reselect_ci_low_db',NaN,'reselect_ci_high_db',NaN, ...
    'oracle_ideal_db',NaN,'oracle_ideal_ci_low_db',NaN,'oracle_ideal_ci_high_db',NaN, ...
    'ideal_thp_mbps',NaN,'uncal_thp_mbps',NaN,'feedback_thp_mbps',NaN,'reselect_thp_mbps',NaN,'oracle_thp_mbps',NaN, ...
    'feedback_thp_gain_mbps',NaN,'feedback_thp_ci_low_mbps',NaN,'feedback_thp_ci_high_mbps',NaN, ...
    'reselect_thp_gain_mbps',NaN,'reselect_thp_ci_low_mbps',NaN,'reselect_thp_ci_high_mbps',NaN, ...
    'port_phase_span_deg',NaN,'port_amp_span_db',NaN,'pair_phase_span_deg',NaN,'pair_amp_span_db',NaN, ...
    'wideband_tx_corr',NaN,'reselect_tpmi_switch_rate',NaN, ...
    'oracle_max_relative_energy_error',NaN,'tx_power_max_error',NaN);
end

function [mu,lo,hi]=pairedCI(x)
mu=mean(x); lo=NaN; hi=NaN; n=numel(x);
if n<2, return; end
v=n-1; q=sqrt(v*(1/betaincinv(.05,v/2,.5)-1));
h=q*std(x,0)/sqrt(n); lo=mu-h; hi=mu+h;
end

function exportStatus(cp,o,fid)
T=table(string({cp.scenario_name}).',string({cp.stage}).',string({cp.status}).', ...
    [cp.completed_mc].',[cp.requested_mc].',[cp.current_mc].',string({cp.checkpoint_file}).', ...
    string({cp.error_identifier}).',string({cp.error_message}).', ...
    'VariableNames',{'scenario','stage','status','completed_mc','requested_mc','current_mc', ...
    'checkpoint_file','error_identifier','error_message'});
safeTable(T,fullfile(o.outputDir,[cp(1).stage '_status.csv']),fid);
end

function exportStage(cp,o,t,fid)
exportStatus(cp,o,fid);
rows=repmat(rowTemplate(),0,1); mc=cell(0,21);
for i=1:numel(cp)
    c=cp(i); rows=[rows;summarize(c,t)]; %#ok<AGROW>
    for j=1:c.completed_mc
        a=c.entries{j};
        for k=1:numel(a.snr_dB)
            g=a.mean_sinr_dB(k,:); h=a.thp_Mbps(k,:);
            mc(end+1,:)={c.scenario_name,c.stage,j,a.seed,a.snr_dB(k),a.selected_slots,a.duty, ...
                g(1)-g(2),g(3)-g(2),g(4)-g(2),g(5)-g(1),h(1),h(2),h(3),h(4),h(5), ...
                max(a.pair_phase_span_deg),max(a.pair_amp_span_dB),a.reselect_tpmi_switch_rate(k), ...
                a.wideband_tx_corr,a.elapsed_s}; %#ok<AGROW>
        end
    end
end
stage=cp(1).stage;
T=struct2table(rows,'AsArray',true);
safeTable(T,fullfile(o.outputDir,[stage '_summary.csv']),fid);
Tmc=cell2table(mc,'VariableNames',{'scenario','stage','mc_index','seed','snr_dB','selected_slots_per_mc','duty', ...
    'loss_db','feedback_db','reselect_db','oracle_ideal_db','ideal_thp_mbps','uncal_thp_mbps', ...
    'feedback_thp_mbps','reselect_thp_mbps','oracle_thp_mbps','pair_phase_span_deg','pair_amp_span_db', ...
    'reselect_tpmi_switch_rate','wideband_tx_corr','elapsed_s'});
safeTable(Tmc,fullfile(o.outputDir,[stage '_mc.csv']),fid);
if strcmp(stage,'confirm')
    vars={'snr_dB','ideal_thp_mbps','uncal_thp_mbps','feedback_thp_mbps','reselect_thp_mbps', ...
        'mc','selected_slots_per_mc','duty','feedback_thp_gain_mbps','feedback_thp_ci_low_mbps', ...
        'feedback_thp_ci_high_mbps','reselect_thp_gain_mbps','reselect_thp_ci_low_mbps','reselect_thp_ci_high_mbps'};
    candidate=T(strcmp(T.scenario,cp(1).scenario_name),vars);
    safeTable(sortrows(candidate,'snr_dB'),fullfile(o.outputDir,'confirm_thp.csv'),fid);
end
end

function atomicSave(path,key,value)
tmp=[tempname(fileparts(path)) '.mat'];
payload=struct(); payload.(key)=value;
save(tmp,'-struct','payload','-v7');
movefile(tmp,path,'f');
end

function safeTable(T,path,fid)
tmp=[tempname(fileparts(path)) '.csv'];
cleanup=onCleanup(@()deleteTemporary(tmp));
try
    writetable(T,tmp,'Encoding','UTF-8');
    movefile(tmp,path,'f');
catch ME
    logLine(fid,'EXPORT WARNING %s: %s\n',path,ME.message);
end
end

function deleteTemporary(path)
if isfile(path), delete(path); end
end

function ts=stamp()
ts=char(datetime('now','TimeZone','Asia/Shanghai','Format','yyyy-MM-dd HH:mm:ss'));
end

function logLine(fid,varargin)
fprintf(1,varargin{:});
fprintf(fid,varargin{:});
end
