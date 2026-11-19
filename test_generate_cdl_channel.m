function tests = test_generate_cdl_channel
%TEST_GENERATE_CDL_CHANNEL Exercise the portable base-MATLAB CDL generator.
tests = functiontests(localfunctions);
end

function testFiveProfilesHaveDistinctDelayAndPowerTables(testCase)
profiles = {'CDL-A','CDL-B','CDL-C','CDL-D','CDL-E'};
tapCounts = [23 23 24 13 14];
isLos = [false false false true true];
kFirstDb = [-Inf -Inf -Inf 13.3 22];
largestDelay = [9.6586 4.7834 8.6523 12.525 20.6419];
delaySpread = 80e-9;
for k = 1:numel(profiles)
    model = generateCdlChannel(profiles{k},'DelaySpreadSec',delaySpread,'Seed',17);
    verifySize(testCase,model.pathGains,[1 tapCounts(k)]);
    verifyNumElements(testCase,model.pathDelays_s,tapCounts(k));
    verifyNumElements(testCase,model.pathPowersLinear,tapCounts(k));
    verifyEqual(testCase,max(model.pathDelays_s),largestDelay(k)*delaySpread,'AbsTol',1e-14);
    verifyEqual(testCase,sum(model.pathPowersLinear),1,'AbsTol',1e-13);
    verifyEqual(testCase,model.hasLOS,isLos(k));
    if isLos(k)
        verifyEqual(testCase,model.firstClusterKFactor_dB,kFirstDb(k),'AbsTol',1e-12);
        verifyGreaterThan(testCase,model.losPowerFraction,0);
        verifyLessThan(testCase,model.losPowerFraction,model.pathPowersLinear(1));
        diffuseFirst = model.pathPowersLinear(1)-model.losPowerFraction;
        verifyEqual(testCase,10*log10(model.losPowerFraction/diffuseFirst), ...
            kFirstDb(k),'AbsTol',1e-10);
    else
        verifyEqual(testCase,model.firstClusterKFactor_dB,-Inf);
        verifyEqual(testCase,model.losDoppler_Hz,0);
        verifyEqual(testCase,model.losPowerFraction,0);
    end
    verifyTrue(testCase,all(isfinite(model.pathGains(:))));
end
end

function testLosPhaseAndDopplerAreSeparateFromDiffuseFading(testCase)
sampleTimes = [0;0.13e-3;1.1e-3];
common = {'SampleTimesSec',sampleTimes,'Seed',29,'SpeedMps',8, ...
    'TravelAnglesDeg',[170 0],'UserUpDirection',[0;0;1]};
for profile = {'CDL-D','CDL-E'}
    a = generateCdlChannel(profile{1},common{:},'LosPhaseRad',0);
    b = generateCdlChannel(profile{1},common{:},'LosPhaseRad',pi);
    % Flipping only the deterministic LOS phase cancels the common diffuse part.
    hLos = double(a.pathGains(:,1)-b.pathGains(:,1))/2;
    verifyEqual(testCase,a.pathGains(:,2:end),b.pathGains(:,2:end));
    verifyEqual(testCase,abs(hLos).^2,repmat(a.losPowerFraction,numel(sampleTimes),1), ...
        'AbsTol',2e-7);
    verifyGreaterThan(testCase,abs(a.losDoppler_Hz),1);
    verifyEqual(testCase,hLos/hLos(1),exp(1i*2*pi*sampleTimes*a.losDoppler_Hz), ...
        'AbsTol',2e-7);
end
end

function testSeededSnapshotsSupportReciprocalDirection(testCase)
sampleTimes = [0;.2e-3;1.05e-3];
common = {'SampleTimesSec',sampleTimes,'Seed',71,'SpeedMps',7, ...
    'TravelAnglesDeg',[35 0],'UserPositionsLambda',[0 0 0;0 0 .22], ...
    'BaseStationPositionsLambda',[0 0 0;0 .5 0;0 0 1.5]};
before = rng;
ul = generateCdlChannel('CDL-B',common{:},'LinkDirection','uplink');
repeated = generateCdlChannel('CDL-B',common{:},'LinkDirection','uplink');
dl = generateCdlChannel('CDL-B',common{:},'LinkDirection','downlink');
after = rng;
verifyEqual(testCase,after,before);
verifyEqual(testCase,ul.sampleTimes_s(:),sampleTimes);
verifySize(testCase,ul.pathGains,[3 23 2 3]);
verifySize(testCase,dl.pathGains,[3 23 3 2]);
verifyEqual(testCase,repeated.pathGains,ul.pathGains);
verifyEqual(testCase,dl.pathGains,permute(ul.pathGains,[1 2 4 3]));
verifyGreaterThan(testCase,norm(reshape(ul.pathGains(1,:,:,:)-ul.pathGains(end,:,:,:),[],1)),1e-4);
verifyLessThanOrEqual(testCase,max(abs(ul.rayDoppler_Hz(:))),ul.maxDoppler_Hz+1e-12);
end

function testZenithOffsetsFoldAtPole(testCase)
% Rotate the CDL-A centers near the south pole. The ray zenith folds at
% 180 degrees, while its azimuth stays fixed (TR 38.901 §7.5 Step 7).
zenith = 178;
seed = 23;
bsPositions = [0 0 0;0 1.5 0];
direction = [0;sind(zenith);cosd(zenith)];
model = generateCdlChannel('CDL-A','Seed',seed,'BaseStationPositionsLambda',bsPositions, ...
    'BsClusterTransform',@(u,p)repmat(direction,1,size(u,2)));

magnitudes = [.0447 .1413 .2492 .3715 .5129 .6797 .8844 1.1481 1.5195 2.1551];
offsets = reshape([magnitudes;-magnitudes],1,[]);
stream = RandStream('mt19937ar','Seed',seed);
zenithOffsets = offsets(randperm(stream,20));
randperm(stream,20); % UE azimuth coupling
randperm(stream,20); % UE zenith coupling
phases = rand(stream,2,20);
rayZenith = zenith+3*zenithOffsets;
verifyGreaterThan(testCase,max(rayZenith),180);
rayZenith = min(rayZenith,360-rayZenith);
rayAzimuth = 90+5*offsets;
y = sind(rayZenith).*sind(rayAzimuth);
expected = sqrt(model.pathPowersLinear(1)/20)* ...
    sum(exp(1i*2*pi*(phases(1,:)+bsPositions(2,2)*y)));
verifyEqual(testCase,double(model.pathGains(1,1,1,2)),expected,'AbsTol',2e-7);
end

function testProfileParametersMatchToolboxWhenAvailable(testCase)
assumeNotEmpty(testCase,which('nrCDLChannel'),'Optional 5G Toolbox reference.');
delaySpread = 80e-9;
for profile = {'CDL-A','CDL-B','CDL-C','CDL-D','CDL-E'}
    model = generateCdlChannel(profile{1},'DelaySpreadSec',delaySpread);
    ref = nrCDLChannel('DelayProfile','Custom','InitialDelayProfile',profile{1}, ...
        'DelaySpread',delaySpread);
    powers = 10.^(ref.AveragePathGains(:)/10);
    powers = powers/sum(powers);
    % Compare sorted delay/power pairs to allow representation order differences.
    rows = sortrows([model.pathDelays_s(:),model.pathPowersLinear(:)]);
    reference = sortrows([ref.PathDelays(:),powers]);
    verifyEqual(testCase,rows(:,1),reference(:,1),'AbsTol',1e-14);
    verifyEqual(testCase,rows(:,2),reference(:,2),'AbsTol',1e-12);
    verifyEqual(testCase,model.hasLOS,ref.HasLOSCluster);
    if ref.HasLOSCluster
        verifyEqual(testCase,model.firstClusterKFactor_dB,ref.KFactorFirstCluster,'AbsTol',1e-12);
    end
end
end
