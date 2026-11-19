function tests = test_ul_native_cdl
%TEST_UL_NATIVE_CDL Standalone CDL-C geometry, scaling, and time evolution.
tests = functiontests(localfunctions);
end

function testProfileScalingAndSnapshotAlignment(testCase)
[cfg,pattern,schedule,ue] = fixture();
model = ul4tNativeCDLC(cfg,pattern,schedule,1,ue);
d = model.diagnostics;
verifyEqual(testCase,d.backend,'base-cdl-c');
verifyEqual(testCase,d.rays_per_cluster,20);
verifyEqual(testCase,size(model.pathGains),[6 24 4 16]);
verifyEqual(testCase,d.path_delays_s([1 2 24]),[0;0.2099;8.6523]*100e-9, ...
    'AbsTol',1e-15);
verifyEqual(testCase,sum(d.cluster_powers_linear),1,'AbsTol',1e-14);
verifyEqual(testCase,model.pRx,d.cluster_powers_linear.*d.cluster_pattern_power, ...
    'AbsTol',1e-14);
verifyEqual(testCase,model.sample_times_ms,(0:.5:2.5).','AbsTol',1e-12);
verifyEqual(testCase,d.response_sample_indices,[2;3;5;6]);
verifyEqual(testCase,model.delayBasis(:,24), ...
    exp(-1i*2*pi*cfg.ofdm.freq_offsets_Hz(:)*8.6523e-7),'AbsTol',1e-12);
verifyEqual(testCase,model.weightedDelayBasis(:,2), ...
    model.delayBasis(:,2)*sqrt(d.cluster_pattern_power(2)),'AbsTol',1e-12);
verifyTrue(testCase,all(isfinite(model.pathGains(:))));
end

function testSeededRealizationAndMobility(testCase)
[cfg,pattern,schedule,ue] = fixture();
rng(101,'twister');
first = ul4tNativeCDLC(cfg,pattern,schedule,1,ue);
rng(101,'twister');
repeat = ul4tNativeCDLC(cfg,pattern,schedule,1,ue);
rng(101,'twister');
other = ul4tNativeCDLC(cfg,pattern,schedule,2,ue);
verifyEqual(testCase,first.pathGains,repeat.pathGains,'AbsTol',single(1e-7));
verifyGreaterThan(testCase,norm(first.pathGains(:)-other.pathGains(:)),1e-3);
verifyGreaterThan(testCase,norm(reshape(first.pathGains(2,:,:,:)- ...
    first.pathGains(end,:,:,:),[],1)),1e-3);
verifyLessThanOrEqual(testCase,max(abs(first.diagnostics.ray_doppler_Hz(:))), ...
    first.max_doppler_Hz+1e-12);
end

function testStaticCoincidentPorts(testCase)
[cfg,pattern,schedule,ue] = fixture();
cfg.ul.ue_speed_kmh = 0;
ue.positions_lambda = zeros(size(ue.positions_lambda));
model = ul4tNativeCDLC(cfg,pattern,schedule,1,ue);
verifyEqual(testCase,size(model.pathGains),[1 24 4 16]);
verifyEqual(testCase,model.sample_times_ms,0);
verifyEqual(testCase,model.diagnostics.response_sample_indices,ones(4,1));
verifyEqual(testCase,model.pathGains(:,:,1,:),model.pathGains(:,:,4,:), ...
    'AbsTol',single(1e-7));
end

function [cfg,pattern,schedule,ue] = fixture()
cfg = configFixedArray();
cfg.Nh = 4;
cfg.Nv = 4;
cfg.pattern.source = 'virtual';
cfg.pattern.file = '';
cfg.pathloss.enabled = false;
cfg.cdl.toolbox_profile = 'CDL-C';
cfg.cdl.toolbox_delay_spread_s = 100e-9;
cfg.ul.ue_speed_kmh = 0.3;
cfg.ul.channel_aging_enabled = true;
cfg = prepareFixedArrayConfig(cfg);
pattern = fixedArrayChannel('pattern',cfg);
schedule = struct('slot_ms',.5,'srs_ms',[.5 2],'ul_ms',[1 2.5]);
ue = struct('layout','upa','positions_lambda',.22*[0 0 0;0 0 1;1 0 0;1 0 1], ...
    'antenna_indices',1:4);
end
