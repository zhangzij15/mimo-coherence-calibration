function [model,scene] = ul4tNativeCDL(cfg,pattern,schedule,mc,ue)
%UL4TNATIVECDL Adapt the portable CDL generator to the calibration simulator.
scene = linkScene('users',cfg);
heading = 360*rand-180;
speed = cfg.ul.ue_speed_kmh/3.6;
maxDoppler = cfg.fc_Hz/cfg.c0*speed;
seed = mod(double(cfg.seed)+mc-1,2^32);

% Project broadside is +y. The CDL array panel has +x broadside.
projectToPanel = [0 1 0;-1 0 0;0 0 1];
mapCenters = @(u,p) projectToPanel*fixedArrayChannel('orient-clusters', ...
    cfg,scene,1,struct('U',u,'p',p));
[row,column] = ndgrid(0:cfg.Nv-1,0:cfg.Nh-1);
bsPositions = [zeros(numel(row),1), ...
    (column(:)-(cfg.Nh-1)/2)*cfg.dh_lambda,((cfg.Nv-1)/2-row(:))*cfg.dv_lambda];
uePositions = ue.positions_lambda;
if strcmp(ue.layout,'upa')
    uePositions = uePositions-(min(uePositions,[],1)+max(uePositions,[],1))/2;
    ueUp = [0;0;-1];
else
    assert(strcmp(ue.layout,'handheld'),'UL256:NativeCDLArray','Unknown UE layout.');
    ueUp = [-1;0;0];
end

samplesPerSlot = max(1,ceil(2*maxDoppler*schedule.slot_ms/1000));
model.sample_step_ms = schedule.slot_ms/samplesPerSlot;
requestedTimes = unique([schedule.srs_ms(:);schedule.ul_ms(:)]);
responseTimes = schedule.srs_ms(:);
if cfg.ul.channel_aging_enabled
    responseTimes = requestedTimes;
end
if maxDoppler==0
    model.sample_times_ms = 0;
    responseIndices = ones(size(responseTimes));
else
    responseIndices = 1+round(responseTimes/model.sample_step_ms);
    model.sample_times_ms = (0:max(responseIndices)-1).'*model.sample_step_ms;
    assert(all(abs(model.sample_times_ms(responseIndices)-responseTimes)<1e-7), ...
        'UL256:NativeCDLTime','SRS/PUSCH times are not on the CDL snapshot grid.');
end

channel = generateCdlChannel(cfg.cdl.toolbox_profile, ...
    CarrierFrequencyHz=cfg.fc_Hz,LightSpeedMps=cfg.c0, ...
    DelaySpreadSec=cfg.cdl.toolbox_delay_spread_s, ...
    SampleTimesSec=model.sample_times_ms/1000, ...
    BaseStationPositionsLambda=bsPositions,UserPositionsLambda=uePositions, ...
    UserUpDirection=ueUp,SpeedMps=speed,TravelAnglesDeg=[heading 0], ...
    Seed=seed,BsClusterTransform=mapCenters);
tau = channel.pathDelays_s;
assert(max(tau)-min(tau)<=cfg.ofdm.cp_duration_s+1e-15,'UL256:NativeCDLDelay', ...
    'The CDL excess delay must fit in the configured cyclic prefix.');
clusterGain = fixedArrayChannel('pattern-power', ...
    projectToPanel.'*channel.bsClusterDirections,pattern);
pathPowerScale = clusterGain(:)*scene.pathloss_power_gain(1);
model.pRx = channel.pathPowersLinear.*pathPowerScale;
assert(sum(model.pRx)>0,'UL256:NativeCDLPattern','All CDL clusters fall in the zero-gain region.');
model.max_doppler_Hz = channel.maxDoppler_Hz;
model.pathGains = channel.pathGains;
model.delayBasis = exp(-1i*2*pi*cfg.ofdm.freq_offsets_Hz(:)*tau.');
model.weightedDelayBasis = model.delayBasis.*sqrt(pathPowerScale).';

d = struct();
d.backend = 'base-cdl';
d.profile = channel.profile;
d.seed = seed;
d.mc = mc;
d.delay_spread_s = cfg.cdl.toolbox_delay_spread_s;
d.path_delays_s = tau;
d.cluster_powers_linear = channel.pathPowersLinear;
d.cluster_pattern_power = clusterGain(:);
d.rays_per_cluster = size(channel.rayDoppler_Hz,2);
d.ray_doppler_Hz = channel.rayDoppler_Hz;
d.has_los = channel.hasLOS;
d.first_cluster_k_factor_dB = channel.firstClusterKFactor_dB;
d.los_power_fraction = channel.losPowerFraction;
d.los_doppler_Hz = channel.losDoppler_Hz;
d.path_gains_size = size(model.pathGains);
d.selected_path_gains_size = size(model.pathGains);
d.sample_times_ms = model.sample_times_ms;
d.requested_times_ms = requestedTimes;
d.response_times_ms = responseTimes;
d.response_sample_indices = responseIndices;
d.sample_step_ms = model.sample_step_ms;
d.snapshot_sample_rate_Hz = 1000/model.sample_step_ms;
d.max_doppler_Hz = model.max_doppler_Hz;
d.ue_layout = ue.layout;
d.selected_physical_ports = ue.antenna_indices;
d.ue_selected_positions_lambda = uePositions;
d.bs_array_size = [cfg.Nv cfg.Nh 1 1 1];
d.bs_spacing_lambda = [cfg.dv_lambda cfg.dh_lambda 0 0];
d.pattern_source = cfg.pattern.source;
d.virtual_hpbw_deg = cfg.pattern.virtual_hpbw_deg;
d.pattern_spherical_integral = pattern.integral;
model.diagnostics = d;
end
