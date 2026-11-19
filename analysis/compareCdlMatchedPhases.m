function summary = compareCdlMatchedPhases(outputDir)
%COMPARECDLMATCHEDPHASES Compare deterministic CDL path gains to 5G Toolbox.
%   Compare 4T256R uplink path gains at 6.7 GHz with explicitly shared
%   diffuse ray coupling and phases and the default LOS reference phase.
%   A control repeats CDL-A with the unrotated cluster angles from the table.
%   MismatchedTaps counts taps with a relative L2 error above 1e-6.
arguments
    outputDir (1,1) string = fullfile(fileparts(mfilename('fullpath')),"cdl_toolbox_comparison")
end
assert(~isempty(which('nrCDLChannel')),'This comparison requires 5G Toolbox.');

profiles = ["CDL-A";"CDL-B";"CDL-C";"CDL-D";"CDL-E";"CDL-A"];
geometries = [repmat("project",5,1);"unrotated"];
seed = 8173;
carrierFrequency = 6.7e9;
delaySpread = 100e-9;
speed = 8;
sampleRate = 2000;
times = (0:40).'/sampleRate;
travelAngles = [17 0];
ueSpacing = .22;
uePositions = ueSpacing*[-.5 0 -.5;-.5 0 .5;.5 0 -.5;.5 0 .5];
cfg = configUL4T('base');
cfg.users.mode = 'specified';
cfg.users.positions_m = [0 60 1.5];
scene = linkScene('users',cfg);
pattern = fixedArrayChannel('pattern',cfg);
projectToPanel = [0 1 0;-1 0 0;0 0 1];
projectCenters = @(u,p) projectToPanel*fixedArrayChannel('orient-clusters', ...
    cfg,scene,1,struct('U',u,'p',p));
[row,column] = ndgrid(0:cfg.Nv-1,0:cfg.Nh-1);
bsPositions = [zeros(numel(row),1), ...
    (column(:)-(cfg.Nh-1)/2)*cfg.dh_lambda,((cfg.Nv-1)/2-row(:))*cfg.dv_lambda];
summary = table();

for k = 1:numel(profiles)
    profile = profiles(k);
    geometry = geometries(k);
    mapCenters = [];
    if geometry=="project", mapCenters = projectCenters; end
    baseArgs = {'CarrierFrequencyHz',carrierFrequency,'DelaySpreadSec',delaySpread, ...
        'SpeedMps',speed,'SampleTimesSec',times,'UserPositionsLambda',uePositions, ...
        'BaseStationPositionsLambda',bsPositions,'UserUpDirection',[0;0;-1], ...
        'TravelAnglesDeg',travelAngles,'BsClusterTransform',mapCenters,'Seed',seed};
    base = generateCdlChannel(profile,baseArgs{:});
    ref = nrCDLChannel('DelayProfile','Custom','InitialDelayProfile',char(profile), ...
        'DelaySpread',delaySpread);
    [baseRows,baseOrder] = sortrows([base.pathDelays_s,base.pathPowersLinear]);
    refPower = 10.^(ref.AveragePathGains(:)/10);
    [refRows,refOrder] = sortrows([ref.PathDelays(:),refPower/sum(refPower)]);
    assert(max(abs(baseRows-refRows),[],'all')<1e-12,'Profile rows differ.');

    nTaps = numel(baseOrder);
    nRays = size(base.rayDoppler_Hz,2);
    [phases,coupling] = matchedRayInputs(seed,nTaps,nRays);
    baseAtRef = zeros(size(refOrder));
    baseAtRef(refOrder) = baseOrder;
    ref.InitialPhases = phases(baseAtRef,:,:);
    ref.RayCoupling = coupling(baseAtRef,:,:);
    bsCenters = base.bsClusterDirections(:,baseAtRef);
    ref.AnglesAoD = atan2d(bsCenters(2,:),bsCenters(1,:));
    ref.AnglesZoD = acosd(max(-1,min(1,bsCenters(3,:))));
    ref.TransmitAntennaArray.Size = [cfg.Nv cfg.Nh 1 1 1];
    ref.TransmitAntennaArray.ElementSpacing = [cfg.dv_lambda cfg.dh_lambda 0 0];
    ref.ReceiveAntennaArray.Size = [2 2 1 1 1];
    ref.ReceiveAntennaArray.ElementSpacing = [ueSpacing ueSpacing 0 0];
    ref.ReceiveArrayOrientation = [90;0;180];
    ref.TransmitAntennaArray.PolarizationAngles = 0;
    ref.ReceiveAntennaArray.PolarizationAngles = 0;
    ref.TransmitAntennaArray.Element = 'isotropic';
    ref.ReceiveAntennaArray.Element = 'isotropic';
    ref.CarrierFrequency = carrierFrequency;
    ref.MaximumDopplerShift = base.maxDoppler_Hz;
    ref.UTDirectionOfTravel = [travelAngles(1);90-travelAngles(2)];
    ref.ChannelFiltering = false;
    ref.ChannelResponseOutput = 'path-gains';
    ref.NormalizePathGains = true;
    ref.NormalizeChannelOutputs = false;
    ref.SampleDensity = Inf;
    ref.SampleRate = sampleRate;
    ref.NumTimeSamples = numel(times);
    swapTransmitAndReceive(ref);
    [refGains,refTimes] = ref();
    assert(max(abs(refTimes(:)-times))<1e-12,'CDL sample times differ.');
    assert(ref.HasLOSCluster==base.hasLOS,'CDL LOS flag differs.');
    assert(isequal(size(refGains),[numel(times) nTaps+double(base.hasLOS) 4 cfg.Nv*cfg.Nh]), ...
        'Unexpected R2024b reference output shape.');

    baseGains = double(base.pathGains);
    baseDiffuse = baseGains;
    refDiffuse = refGains;
    losError = NaN;
    losPhaseOffset = NaN;
    if base.hasLOS
        pathInfo = info(ref);
        assert(strcmp(pathInfo.ClusterTypes{1},'LOS') && strcmp(pathInfo.ClusterTypes{2},'NLOS') && ...
            abs(pathInfo.PathDelays(1)-pathInfo.PathDelays(2))<1e-15,'Unexpected LOS path ordering.');

        % A pi shift changes only the coherent LOS term. Compare the two
        % isolated LOS terms directly at the unadjusted reference phase.
        shifted = generateCdlChannel(profile,baseArgs{:},'LosPhaseRad',pi);
        baseLos = (baseGains(:,1,:,:)-double(shifted.pathGains(:,1,:,:)))/2;
        baseDiffuse(:,1,:,:) = baseGains(:,1,:,:)-baseLos;
        refLos = refGains(:,1,:,:);
        refDiffuse = refGains(:,2:end,:,:);
        losPhaseOffset = angle(refLos(1)/baseLos(1));
        losError = relativeError(baseLos,refLos);
        refGains = refDiffuse;
        refGains(:,1,:,:) = refDiffuse(:,1,:,:)+refLos;
    end

    diffuseError = relativeError(baseDiffuse(:,baseOrder,:,:),refDiffuse(:,refOrder,:,:));
    refGains = refGains(:,refOrder,:,:);
    baseGains = baseGains(:,baseOrder,:,:);
    totalError = relativeError(baseGains,refGains);
    maxAbsError = max(abs(baseGains-refGains),[],'all');
    tapErrorEnergy = sum(abs(baseGains-refGains).^2,[1 3 4]);
    refTapEnergy = sum(abs(refGains).^2,[1 3 4]);
    tapError = sqrt(tapErrorEnergy./refTapEnergy);
    differingTaps = nnz(tapError>1e-6);
    patternPower = fixedArrayChannel('pattern-power', ...
        projectToPanel.'*base.bsClusterDirections(:,baseOrder),pattern);
    weightedError = sqrt(sum(tapErrorEnergy.*patternPower(:).')/sum(refTapEnergy.*patternPower(:).'));
    [worstTapError,worstSortedTap] = max(tapError);
    worstTap = baseOrder(worstSortedTap);
    worstTapPower = 100*base.pathPowersLinear(worstTap);
    worstCenterZo = acosd(max(-1,min(1,base.bsClusterDirections(3,worstTap))));
    summary = [summary;table(profile,geometry,nTaps,diffuseError,losError, ...
        rad2deg(losPhaseOffset),totalError,weightedError,maxAbsError,differingTaps, ...
        worstTap,worstTapError,worstTapPower,worstCenterZo, ...
        'VariableNames',{'Profile','Geometry','Taps','DiffuseRelL2','LosRelL2', ...
        'ToolboxLosPhaseOffset_deg','TotalRelL2','PatternWeightedRelL2','MaxAbsCoefficientError', ...
        'MismatchedTaps', ...
        'WorstTap','WorstTapRelL2','WorstTapNominalPower_pct','WorstTapCenterZenith_deg'})]; %#ok<AGROW>
end

disp(summary);
if ~isfolder(outputDir), mkdir(outputDir); end
writetable(summary,fullfile(outputDir,'cdl_matched_phases.csv'));
end

function [phases,coupling] = matchedRayInputs(seed,nTaps,nRays)
stream = RandStream('mt19937ar','Seed',seed);
phases = zeros(nTaps,nRays,4);
coupling = zeros(nTaps,nRays,3);
for tap = 1:nTaps
    % Replay the three native independent permutations. Toolbox couples
    % the two zenith lists first, then attaches them to the ordered BS AoD.
    bsZenithPermutation = randperm(stream,nRays);
    ueAzimuthPermutation = randperm(stream,nRays);
    ueZenithPermutation = randperm(stream,nRays);
    phaseDraws = rand(stream,2,nRays);
    coupling(tap,:,1) = ueAzimuthPermutation;
    coupling(tap,bsZenithPermutation,2) = ueZenithPermutation;
    coupling(tap,:,3) = bsZenithPermutation;
    phases(tap,:,1) = mod(360*phaseDraws(1,:)+180,360)-180;
    phases(tap,:,3) = mod(360*phaseDraws(2,:)+180,360)-180;
end
end

function error = relativeError(actual,reference)
error = norm(actual(:)-reference(:))/norm(reference(:));
end
