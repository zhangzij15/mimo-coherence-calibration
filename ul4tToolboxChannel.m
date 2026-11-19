function [model,scene] = ul4tToolboxChannel(cfg,pattern,schedule,mc,ue)
%UL4TTOOLBOXCHANNEL 5G Toolbox CDL coefficients at the NR slot sampling rate.
assert(nargin>=5 && ~isempty(ue), ...
    'Pass the UE-array description explicitly as the fifth input.');
% Start with the standard DL cluster definitions and reverse the link to UL.
channel = nrCDLChannel('DelayProfile','Custom', ...
    'InitialDelayProfile',cfg.cdl.toolbox_profile, ...
    'DelaySpread',cfg.cdl.toolbox_delay_spread_s);
tau = channel.PathDelays(:);
p = 10.^(channel.AveragePathGains(:)/10);
p = p/sum(p);
assert(max(tau)-min(tau)<=cfg.ofdm.cp_duration_s+1e-15, ...
    'Toolbox CDL excess delay exceeds the configured CP. Reduce toolbox_delay_spread_s.');

az = channel.AnglesAoD(:).';
zo = channel.AnglesZoD(:).';
profile.U = [sind(zo).*cosd(az);sind(zo).*sind(az);cosd(zo)];
profile.p = p;
scene = linkScene('users',cfg);
u = fixedArrayChannel('orient-clusters',cfg,scene,1,profile);
gain = fixedArrayChannel('pattern-power',u,pattern);
pathPowerScale = gain(:)*scene.pathloss_power_gain(1);
model.pRx = p.*pathPowerScale;
assert(sum(model.pRx)>0,'All Toolbox CDL clusters fall in a zero-gain pattern region.');

% Toolbox array broadside is +x, while the project's BS broadside is +y.
channel.AnglesAoD = atan2d(-u(1,:),u(2,:));
channel.AnglesZoD = acosd(max(-1,min(1,u(3,:))));
channel.TransmitAntennaArray = makePanel(channel.TransmitAntennaArray, ...
    cfg.Nv,cfg.Nh,cfg.dv_lambda,cfg.dh_lambda);
if strcmp(ue.layout,'handheld')
    % The eight handset positions form the perimeter of a 3x3 XY grid.
    % Rotate the native YZ panel so +z maps to -x and the broadside to +z.
    spacing = ue.toolbox_spacing_lambda;
    channel.ReceiveAntennaArray = makePanel(channel.ReceiveAntennaArray,3,3,spacing(1),spacing(2));
    channel.ReceiveArrayOrientation = [0;-90;0];
else
    channel.ReceiveAntennaArray = makeRegularUEPanel(channel.ReceiveAntennaArray,ue.positions_lambda);
    % Native panel rows run along -z and columns along +y. This proper
    % rotation maps those port increments to the configured +z and +x.
    channel.ReceiveArrayOrientation = [90;0;180];
end
nGridPorts = prod(channel.ReceiveAntennaArray.Size);
channel.CarrierFrequency = cfg.fc_Hz;
model.max_doppler_Hz = cfg.fc_Hz/cfg.c0*(cfg.ul.ue_speed_kmh/3.6);
channel.MaximumDopplerShift = model.max_doppler_Hz;
channel.UTDirectionOfTravel = [360*rand-180;90];
channel.NormalizePathGains = true;
channel.NormalizeChannelOutputs = false; % retain the BS array gain
channel.RandomStream = 'mt19937ar with seed';
channel.Seed = mod(double(cfg.seed)+mc-1,2^32);
channel.ChannelFiltering = false;
channel.OutputDataType = 'single';
channel.SampleDensity = Inf;
samplesPerSlot = max(1,ceil(2*model.max_doppler_Hz*schedule.slot_ms/1000));
model.sample_step_ms = schedule.slot_ms/samplesPerSlot;
channel.SampleRate = 1000/model.sample_step_ms;
lastTime_ms = max(schedule.srs_ms);
if cfg.ul.channel_aging_enabled
    lastTime_ms = max(lastTime_ms,max(schedule.ul_ms));
end
numSamples = 1+round(lastTime_ms/model.sample_step_ms);
if model.max_doppler_Hz==0, numSamples = 1; end
channel.NumTimeSamples = numSamples;
swapTransmitAndReceive(channel);
[pathGains,sampleTimes] = channel();
if channel.HasLOSCluster && size(pathGains,2)==numel(tau)+1
    % The custom profile stores aggregate first-cluster power. R2024b
    % returns its specular and diffuse coefficients separately at one delay.
    pathInfo = info(channel);
    types = pathInfo.ClusterTypes;
    observedDelays = pathInfo.PathDelays(:);
    assert(numel(observedDelays)==numel(tau)+1 && ...
        strcmp(types{1},'LOS') && strcmp(types{2},'NLOS') && ...
        observedDelays(1)==observedDelays(2) && ...
        all(abs(observedDelays(2:end)-tau)<1e-15), ...
        'UL256:ToolboxCDLLOS','Unexpected Toolbox LOS path ordering.');
    pathGains = cat(2,pathGains(:,1,:,:)+pathGains(:,2,:,:),pathGains(:,3:end,:,:));
end
assert(size(pathGains,2)==numel(tau) && ...
    size(pathGains,3)==nGridPorts && size(pathGains,4)==cfg.N, ...
    'Toolbox CDL output dimensions must be time x cluster x UE grid x BS Rx.');
assert(numel(sampleTimes)>=numSamples && size(pathGains,1)>=numSamples, ...
    'Toolbox CDL did not cover the requested simulation times.');
expectedTimes = (0:numSamples-1)'/channel.SampleRate;
assert(max(abs(double(sampleTimes(1:numSamples))-expectedTimes))<1e-10, ...
    'Toolbox CDL path gains are not aligned to the requested snapshot times.');

model.pathGains = pathGains(:,:,ue.toolbox_ports,:);
assert(all(isfinite(model.pathGains(:))) && any(abs(model.pathGains(:))>0), ...
    'Selected Toolbox CDL path gains must be finite and nonzero.');
model.sample_times_ms = double(sampleTimes(:))*1000;
requestedTimes_ms = unique([schedule.srs_ms(:);schedule.ul_ms(:)]);
alignmentTolerance_ms = 1e-7;
if model.max_doppler_Hz>0
    requestedIndices = 1+round(requestedTimes_ms/model.sample_step_ms);
    assert(all(abs((requestedIndices-1)*model.sample_step_ms-requestedTimes_ms) ...
        <=alignmentTolerance_ms), ...
        'SRS/PUSCH times must align exactly with the Toolbox snapshot grid.');
    % Without propagation aging, PUSCH reuses its SRS channel snapshot.
    responseTimes_ms = schedule.srs_ms(:);
    if cfg.ul.channel_aging_enabled, responseTimes_ms = requestedTimes_ms; end
    responseIndices = 1+round(responseTimes_ms/model.sample_step_ms);
    assert(all(responseIndices>=1 & responseIndices<=size(model.pathGains,1)), ...
        'Toolbox snapshots do not cover every requested channel response time.');
    assert(all(abs(model.sample_times_ms(responseIndices)-responseTimes_ms) ...
        <=alignmentTolerance_ms), ...
        'Actual Toolbox sample times do not match the requested response times.');
else
    % A zero-Doppler channel is constant; all target times use snapshot one.
    responseTimes_ms = requestedTimes_ms;
    responseIndices = ones(size(responseTimes_ms));
    assert(abs(model.sample_times_ms(1))<=alignmentTolerance_ms, ...
        'A static Toolbox channel must provide its initial snapshot at time zero.');
end
model.delayBasis = exp(-1i*2*pi*cfg.ofdm.freq_offsets_Hz(:)*tau.');
model.weightedDelayBasis = bsxfun(@times,model.delayBasis,sqrt(pathPowerScale).');

% Keep metadata small: no channel coefficients are duplicated in diagnostics.
model.diagnostics.backend = '5g-toolbox';
model.diagnostics.profile = cfg.cdl.toolbox_profile;
model.diagnostics.has_los = channel.HasLOSCluster;
model.diagnostics.first_cluster_k_factor_dB = -Inf;
if channel.HasLOSCluster
    model.diagnostics.first_cluster_k_factor_dB = channel.KFactorFirstCluster;
end
model.diagnostics.delay_spread_s = cfg.cdl.toolbox_delay_spread_s;
model.diagnostics.seed = channel.Seed;
model.diagnostics.mc = mc;
model.diagnostics.path_gains_size = size(pathGains);
model.diagnostics.selected_path_gains_size = size(model.pathGains);
model.diagnostics.sample_times_ms = model.sample_times_ms;
model.diagnostics.requested_times_ms = requestedTimes_ms;
model.diagnostics.response_times_ms = responseTimes_ms;
model.diagnostics.response_sample_indices = responseIndices;
model.diagnostics.sample_step_ms = model.sample_step_ms;
model.diagnostics.snapshot_sample_rate_Hz = channel.SampleRate;
model.diagnostics.max_doppler_Hz = model.max_doppler_Hz;
model.diagnostics.selected_toolbox_ports = ue.toolbox_ports;
model.diagnostics.selected_physical_ports = ue.antenna_indices;
model.diagnostics.ue_layout = ue.layout;
% These are the actual UL array roles after swapTransmitAndReceive.
model.diagnostics.bs_array_size = channel.ReceiveAntennaArray.Size;
model.diagnostics.ue_grid_size = channel.TransmitAntennaArray.Size;
model.diagnostics.bs_spacing_lambda = channel.ReceiveAntennaArray.ElementSpacing;
model.diagnostics.ue_grid_spacing_lambda = channel.TransmitAntennaArray.ElementSpacing;
model.diagnostics.bs_orientation_deg = channel.ReceiveArrayOrientation;
model.diagnostics.ue_orientation_deg = channel.TransmitArrayOrientation;
% Same direction transform used above: project BS +y becomes Toolbox +x.
model.diagnostics.bs_direction_project_to_toolbox = [0 1 0;-1 0 0;0 0 1];
% Derive coordinates from the applied public array properties. Tests compare
% these values against the Toolbox's own coordinate construction separately.
gridPositions = panelPositions(channel.TransmitAntennaArray,channel.TransmitArrayOrientation);
selectedPositions = gridPositions(ue.toolbox_ports,:);
actualRelative = bsxfun(@minus,selectedPositions,selectedPositions(1,:));
configuredRelative = bsxfun(@minus,ue.positions_lambda,ue.positions_lambda(1,:));
model.diagnostics.ue_grid_positions_lambda = gridPositions;
model.diagnostics.ue_selected_positions_lambda = selectedPositions;
model.diagnostics.ue_selected_relative_positions_lambda = actualRelative;
model.diagnostics.ue_configured_relative_positions_lambda = configuredRelative;
model.diagnostics.ue_geometry_error_lambda = norm(actualRelative-configuredRelative,'fro');
model.diagnostics.ue_coordinate_frame = ...
    'Toolbox GCS, centered grid, selected coordinates in logical TX-port order';
assert(model.diagnostics.ue_geometry_error_lambda<1e-10, ...
    'Configured UE positions do not match the applied Toolbox panel and port order.');
model.diagnostics.normalize_path_gains = channel.NormalizePathGains;
model.diagnostics.normalize_channel_outputs = channel.NormalizeChannelOutputs;
model.diagnostics.transmit_receive_swapped = channel.TransmitAndReceiveSwapped;
model.diagnostics.pattern_source = cfg.pattern.source;
model.diagnostics.virtual_hpbw_deg = cfg.pattern.virtual_hpbw_deg;
model.diagnostics.path_delays_s = tau;
model.diagnostics.cluster_pattern_power = gain(:);
model.diagnostics.pattern_spherical_integral = pattern.integral;
end

function array = makeRegularUEPanel(array,positions)
pos = bsxfun(@minus,positions,positions(1,:));
dv = pos(2,3);
if size(pos,1)==2
    dh = pos(2,1);
    if dv>0 && norm(pos(2,1:2))<1e-10
        array = makePanel(array,2,1,dv,.5);
    elseif dh>0 && norm(pos(2,2:3))<1e-10
        array = makePanel(array,1,2,.5,dh);
    else
        error('Toolbox mode needs a vertical or horizontal pair for a two-port regular UE array.');
    end
else
    dh = pos(3,1);
    assert(dv>0 && dh>0 && norm(pos-[0 0 0;0 0 dv;dh 0 0;dh 0 dv],'fro')<1e-10, ...
        'Toolbox mode requires a rectangular 2x2 UE panel, with vertical-first port order.');
    array = makePanel(array,2,2,dv,dh);
end
end

function array = makePanel(array,nVertical,nHorizontal,dVertical,dHorizontal)
array.Size = [nVertical nHorizontal 1 1 1];
array.ElementSpacing = [dVertical dHorizontal 0 0];
array.PolarizationAngles = 0;
array.Element = 'isotropic';
end

function positions = panelPositions(array,orientation)
% One-polarization, one-panel geometry in wavelengths, using TR 38.901 ZYX
% bearing/downtilt/slant. No Toolbox internal API is used in production.
assert(isequal(array.Size(3:5),[1 1 1]), ...
    'Coordinate diagnostics require a single-polarization, single-panel array.');
[row,column] = ndgrid(0:array.Size(1)-1,0:array.Size(2)-1);
native = [zeros(numel(row),1), ...
    (column(:)-(array.Size(2)-1)/2)*array.ElementSpacing(2), ...
    ((array.Size(1)-1)/2-row(:))*array.ElementSpacing(1)];
c = cosd(orientation(:));
s = sind(orientation(:));
Rz = [c(1) -s(1) 0;s(1) c(1) 0;0 0 1];
Ry = [c(2) 0 s(2);0 1 0;-s(2) 0 c(2)];
Rx = [1 0 0;0 c(3) -s(3);0 s(3) c(3)];
positions = (Rz*Ry*Rx*native.').';
end
