function tests = test_ul_toolbox_geometry
%TEST_UL_TOOLBOX_GEOMETRY Real Toolbox coordinates, port order and snapshots.
% Run: results = runtests('test_ul_toolbox_geometry.m'); assertSuccess(results)
% Toolbox internals are used only here to inspect the real panel coordinates.
tests = functiontests(localfunctions);
end

function setup(testCase)
assumeNotEmpty(testCase,which('nrCDLChannel'),'Optional 5G Toolbox reference.');
end

function testUpaFourPortCoordinates(testCase)
% A missing/wrong UE orientation reverses vertical ports and swaps X/Y.
[cfg,pattern,schedule] = fixture();
ue = regularUE([0 0 0;0 0 .7;.4 0 0;.4 0 .7]);
model = ul4tToolboxChannel(cfg,pattern,schedule,1,ue);
actual = toolboxUECoordinates(model.diagnostics);
verifyEqual(testCase,relativeToFirst(actual), ...
    [0 0 0;0 0 .7;.4 0 0;.4 0 .7],'AbsTol',1e-12);
verifyEqual(testCase,model.diagnostics.ue_grid_size,[2 2 1 1 1]);
verifyEqual(testCase,size(model.pathGains,3),4);
verifyEqual(testCase,size(model.pathGains,4),64);
end

function testUpaTwoPortVerticalCoordinates(testCase)
[cfg,pattern,schedule] = fixture();
ue = regularUE([0 0 0;0 0 .7]);
model = ul4tToolboxChannel(cfg,pattern,schedule,1,ue);
verifyEqual(testCase,relativeToFirst(toolboxUECoordinates(model.diagnostics)), ...
    [0 0 0;0 0 .7],'AbsTol',1e-12);
verifyEqual(testCase,model.diagnostics.ue_grid_size,[2 1 1 1 1]);
end

function testUpaTwoPortHorizontalCoordinates(testCase)
[cfg,pattern,schedule] = fixture();
ue = regularUE([0 0 0;.4 0 0]);
model = ul4tToolboxChannel(cfg,pattern,schedule,1,ue);
verifyEqual(testCase,relativeToFirst(toolboxUECoordinates(model.diagnostics)), ...
    [0 0 0;.4 0 0],'AbsTol',1e-12);
verifyEqual(testCase,model.diagnostics.ue_grid_size,[1 2 1 1 1]);
end

function testHandheldPortOrderIsPreserved(testCase)
[cfg,pattern,schedule] = fixture();
lambda = cfg.c0/cfg.fc_Hz;
ue = handheldUE(lambda);
model = ul4tToolboxChannel(cfg,pattern,schedule,1,ue);
actual = toolboxUECoordinates(model.diagnostics)*lambda;
verifyEqual(testCase,actual, ...
    [-.075 -.035 0;-.075 .035 0;.075 -.035 0;.075 .035 0], ...
    'AbsTol',1e-12);
verifyEqual(testCase,model.diagnostics.selected_toolbox_ports,[1 7 3 9]);
verifyEqual(testCase,model.diagnostics.selected_physical_ports,[1 7 3 5]);
verifyEqual(testCase,model.diagnostics.ue_grid_size,[3 3 1 1 1]);
verifyEqual(testCase,model.diagnostics.ue_orientation_deg,[0;-90;0]);
end

function testCoordinateDiagnosticsMatchToolbox(testCase)
[cfg,pattern,schedule] = fixture();
for layout = {'upa','handheld'}
    if strcmp(layout{1},'upa')
        ue = regularUE([0 0 0;0 0 .7;.4 0 0;.4 0 .7]);
    else
        ue = handheldUE(cfg.c0/cfg.fc_Hz);
    end
    model = ul4tToolboxChannel(cfg,pattern,schedule,1,ue);
    d = model.diagnostics;
    actual = toolboxUECoordinates(d);
    verifyTrue(testCase,isfield(d,'ue_selected_positions_lambda'), ...
        'Channel diagnostics must preserve the actual ordered UE coordinates.');
    if ~isfield(d,'ue_selected_positions_lambda'), continue; end
    verifyEqual(testCase,d.ue_selected_positions_lambda,actual,'AbsTol',1e-12);
    verifyEqual(testCase,d.ue_selected_relative_positions_lambda, ...
        relativeToFirst(actual),'AbsTol',1e-12);
    verifyEqual(testCase,d.ue_configured_relative_positions_lambda, ...
        relativeToFirst(ue.positions_lambda),'AbsTol',1e-12);
    verifyLessThanOrEqual(testCase,d.ue_geometry_error_lambda,1e-12);
end
end

function testUpaPoseAndBSDirectionTransformAreProperRotations(testCase)
[cfg,pattern,schedule] = fixture();
ue = regularUE([0 0 0;0 0 .7;.4 0 0;.4 0 .7]);
model = ul4tToolboxChannel(cfg,pattern,schedule,1,ue);
d = model.diagnostics;
% Use the Toolbox's own Euler convention, independent of project helpers.
R = wireless.internal.channelmodels.rotationMatrix(d.ue_orientation_deg.');
verifyEqual(testCase,R*[0;0;-1],[0;0;1],'AbsTol',1e-12);
verifyEqual(testCase,R*[0;1;0],[1;0;0],'AbsTol',1e-12);
verifyEqual(testCase,R.'*R,eye(3),'AbsTol',1e-12);
verifyEqual(testCase,det(R),1,'AbsTol',1e-12);
verifyEqual(testCase,d.bs_orientation_deg,[0;0;0]);
verifyTrue(testCase,isfield(d,'bs_direction_project_to_toolbox'));
if isfield(d,'bs_direction_project_to_toolbox')
    B = d.bs_direction_project_to_toolbox;
    verifyEqual(testCase,B*[0;1;0],[1;0;0],'AbsTol',1e-12);
    verifyEqual(testCase,B*[1;0;0],[0;-1;0],'AbsTol',1e-12);
    verifyEqual(testCase,B.'*B,eye(3),'AbsTol',1e-12);
    verifyEqual(testCase,det(B),1,'AbsTol',1e-12);
end
end

function testNlosProfilesAndZeroDopplerSnapshots(testCase)
[cfg,pattern,schedule] = fixture();
cfg.ul.ue_speed_kmh = 0;
ue = regularUE([0 0 0;0 0 .7;.4 0 0;.4 0 .7]);
for profile = {'CDL-A','CDL-B','CDL-C'}
    cfg.cdl.toolbox_profile = profile{1};
    model = ul4tToolboxChannel(cfg,pattern,schedule,1,ue);
    verifyEqual(testCase,size(model.pathGains,1),1);
    verifyEqual(testCase,size(model.pathGains,2),numel(model.diagnostics.path_delays_s));
    verifyEqual(testCase,size(model.pathGains,3),4);
    verifyEqual(testCase,size(model.pathGains,4),64);
    verifyEqual(testCase,model.sample_times_ms,0,'AbsTol',1e-12);
    verifyTrue(testCase,all(model.diagnostics.response_sample_indices==1));
    verifyTrue(testCase,all(isfinite(model.pathGains(:))));
    verifyGreaterThan(testCase,norm(model.pathGains(:)),0);
    verifyTrue(testCase,all(isfinite(model.weightedDelayBasis(:))));
end
end

function testLosReferenceMergesTheCoDelayedTap(testCase)
[cfg,pattern,schedule] = fixture();
cfg.ul.ue_speed_kmh = 0;
ue = regularUE([0 0 0;0 0 .7;.4 0 0;.4 0 .7]);
profiles = {'CDL-D','CDL-E'};
for k = 1:numel(profiles)
    cfg.cdl.toolbox_profile = profiles{k};
    model = ul4tToolboxChannel(cfg,pattern,schedule,1,ue);
    verifySize(testCase,model.pathGains,[1 12+k 4 64]);
    verifyEqual(testCase,model.diagnostics.path_delays_s(1:2)==0,[true;false]);
    verifyTrue(testCase,all(isfinite(model.pathGains),'all'));
end
end

function [cfg,pattern,schedule] = fixture()
cfg = configUL4T('toolbox');
cfg = prepareFixedArrayConfig(cfg);
pattern = fixedArrayChannel('pattern',cfg);
schedule = struct('slot_ms',.5,'srs_ms',[.5 2],'ul_ms',[1 2.5]);
end

function ue = regularUE(positions)
ue = struct('layout','upa','positions_lambda',positions, ...
    'antenna_indices',1:size(positions,1),'toolbox_ports',1:size(positions,1));
end

function ue = handheldUE(lambda)
ue = struct('layout','handheld','antenna_indices',[1 7 3 5], ...
    'toolbox_ports',[1 7 3 9],'toolbox_spacing_lambda',[.075 .035]/lambda, ...
    'positions_lambda',[-.075 -.035 0;-.075 .035 0;.075 -.035 0;.075 .035 0]/lambda);
end

function coordinates = toolboxUECoordinates(d)
panel = struct('Size',d.ue_grid_size, ...
    'ElementSpacing',d.ue_grid_spacing_lambda,'PolarizationAngles',0);
array = wireless.internal.channelmodels.makeAntennaArray(panel,d.ue_orientation_deg);
grid = reshape(array.ElementPositions,3,[]).';
coordinates = grid(d.selected_toolbox_ports,:);
end

function relative = relativeToFirst(coordinates)
relative = bsxfun(@minus,coordinates,coordinates(1,:));
end
