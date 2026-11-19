function model = generateCdlChannel(profile,options)
%GENERATECDLCHANNEL Single-file Base MATLAB CDL-A through CDL-E path gains.
%   M = generateCdlChannel('CDL-D',SampleTimesSec=0:0.0005:0.01) returns
%   pathGains(time,tap,tx,rx), pathDelays_s, and nominal pathPowersLinear.
%   Array positions are in carrier wavelengths in the coordinate frames of
%   the profile's BS and UE angles. The BS elements are theta-polarized.
%   UserUpDirection sets the UE's local theta polarization axis.
%   LinkDirection defaults to 'uplink'. 'downlink' transposes the antenna
%   dimensions of the reciprocal realization.
%   LosPhaseRad sets the reference LOS phase at t=0 for CDL-D and CDL-E.
%   BsClusterTransform is an optional function Unew = f(U,p) for mapping
%   the 3-by-N BS center directions, with p the N nominal tap powers.
%
%   Profile and ray data: 3GPP TR 38.901 V18.0.0, Tables 7.7.1-1..5,
%   Table 7.5-3, and Clause 7.7.3. Default LOS powers use the table values.
arguments
    profile (1,1) string = "CDL-C"
    options.CarrierFrequencyHz (1,1) double {mustBeReal,mustBePositive,mustBeFinite} = 6.7e9
    options.DelaySpreadSec (1,1) double {mustBeReal,mustBeNonnegative,mustBeFinite} = 100e-9
    options.SampleTimesSec double {mustBeVector,mustBeNonempty,mustBeReal,mustBeNonnegative,mustBeFinite} = 0
    options.BaseStationPositionsLambda (:,3) double {mustBeNonempty,mustBeReal,mustBeFinite} = [0 0 0]
    options.UserPositionsLambda (:,3) double {mustBeNonempty,mustBeReal,mustBeFinite} = [0 0 0]
    options.SpeedMps (1,1) double {mustBeReal,mustBeNonnegative,mustBeFinite} = 0
    options.TravelAnglesDeg (1,2) double {mustBeReal,mustBeFinite} = [0 0]
    options.UserUpDirection (3,1) double {mustBeReal,mustBeFinite} = [0;0;1]
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 1
    options.LinkDirection (1,1) string = "uplink"
    options.LosPhaseRad (1,1) double {mustBeReal,mustBeFinite} = 0
    options.BsClusterTransform = []
    options.LightSpeedMps (1,1) double {mustBeReal,mustBePositive,mustBeFinite} = 299792458
end

profile = validatestring(upper(profile),{'CDL-A','CDL-B','CDL-C','CDL-D','CDL-E'});
direction = validatestring(options.LinkDirection,{'uplink','downlink'});
assert(options.Seed<2^32,'CDL:Seed','Seed must be between 0 and 2^32-1.');
upNorm = norm(options.UserUpDirection);
assert(upNorm>0,'CDL:Polarization','UserUpDirection must be nonzero.');
ueUp = options.UserUpDirection/upNorm;
assert(abs(options.TravelAnglesDeg(2))<=90,'CDL:TravelAngle', ...
    'The elevation of the UE velocity must be between -90 and 90 degrees.');

[rows,spreads,xprDb,hasLos] = cdlProfileData(profile);
rawPower = 10.^(rows(:,2)/10);
power = rawPower/sum(rawPower);
losPower = 0;
firstK = -Inf;
if hasLos
    % The first two table rows are co-located specular and diffuse power.
    % Generate them separately and return their coherent sum in one tap.
    losPower = power(1);
    firstK = rows(1,2)-rows(2,2);
    power = power(2:end);
    rows = rows(2:end,:);
end
diffusePower = power;
power(1) = power(1)+losPower;
tapCount = size(rows,1);
center = cdlUnitVectors(rows(:,3).',rows(:,5).');
if ~isempty(options.BsClusterTransform)
    assert(isa(options.BsClusterTransform,'function_handle'),'CDL:ClusterTransform', ...
        'BsClusterTransform must be a function handle.');
    center = options.BsClusterTransform(center,power);
    assert(isnumeric(center) && isreal(center) && isequal(size(center),[3 tapCount]) && ...
        all(isfinite(center),'all') && all(abs(vecnorm(center)-1)<1e-9), ...
        'CDL:ClusterTransform','BsClusterTransform must return one unit vector per tap.');
end
bsAzimuth = atan2d(center(2,:),center(1,:));
bsZenith = acosd(max(-1,min(1,center(3,:))));

times = options.SampleTimesSec(:);
bsPositions = options.BaseStationPositionsLambda;
uePositions = options.UserPositionsLambda;
ueCount = size(uePositions,1);
bsCount = size(bsPositions,1);
speedAngles = options.TravelAnglesDeg;
velocityDirection = cdlUnitVectors(speedAngles(1),90-speedAngles(2));
maxDoppler = options.CarrierFrequencyHz/options.LightSpeedMps*options.SpeedMps;
stream = RandStream('mt19937ar','Seed',options.Seed);
offsetMagnitudes = [.0447 .1413 .2492 .3715 .5129 .6797 .8844 1.1481 1.5195 2.1551];
rayOffsets = reshape([offsetMagnitudes;-offsetMagnitudes],1,[]);
rayCount = numel(rayOffsets);
crossPolarAmplitude = 10^(-xprDb/20);

pathGains = complex(zeros(numel(times),tapCount,ueCount,bsCount,'single'));
rayDoppler = zeros(tapCount,rayCount);
losDoppler = 0;
for k = 1:tapCount
    % CDL tables describe DL. The UE incidence angles become departure
    % angles for the reciprocal UL while still setting UE Doppler.
    bsRayAz = bsAzimuth(k)+spreads(1)*rayOffsets;
    bsRayZo = foldCdlZenith(bsZenith(k)+spreads(3)*rayOffsets(randperm(stream,rayCount)));
    ueRayAz = rows(k,4)+spreads(2)*rayOffsets(randperm(stream,rayCount));
    ueRayZo = foldCdlZenith(rows(k,6)+spreads(4)*rayOffsets(randperm(stream,rayCount)));
    bsDirection = cdlUnitVectors(bsRayAz,bsRayZo);
    ueDirection = cdlUnitVectors(ueRayAz,ueRayZo);
    rayDoppler(k,:) = maxDoppler*(velocityDirection.'*ueDirection);

    [fTheta,fPhi] = cdlUeField(ueRayAz,ueRayZo,ueDirection,ueUp);
    phases = exp(1i*2*pi*rand(stream,2,rayCount));
    rayGain = sqrt(diffusePower(k)/rayCount)*(fTheta.*phases(1,:) + ...
        crossPolarAmplitude*fPhi.*phases(2,:));
    tx = exp(1i*2*pi*(uePositions*ueDirection));
    rx = exp(1i*2*pi*(bsPositions*bsDirection));
    spatial = reshape(reshape(tx,ueCount,1,rayCount).*reshape(rx,1,bsCount,rayCount), ...
        ueCount*bsCount,rayCount);
    temporal = exp(1i*2*pi*(times*rayDoppler(k,:))).*rayGain;
    tap = temporal*spatial.';

    if k==1 && hasLos
        % TR 38.901 (7.5-29). The LOS polarization matrix is diag(1,-1).
        % A theta-polarized BS selects its first row in this model.
        ueLos = cdlUnitVectors(rows(1,4),rows(1,6));
        losDoppler = maxDoppler*(velocityDirection.'*ueLos);
        fLos = cdlUeField(rows(1,4),rows(1,6),ueLos,ueUp);
        txLos = exp(1i*2*pi*(uePositions*ueLos));
        rxLos = exp(1i*2*pi*(bsPositions*center(:,1)));
        spatialLos = reshape(txLos*rxLos.',1,[]);
        temporalLos = exp(1i*(options.LosPhaseRad+2*pi*(times*losDoppler)));
        tap = tap+sqrt(losPower)*fLos*(temporalLos*spatialLos);
    end
    pathGains(:,k,:,:) = reshape(tap,numel(times),1,ueCount,bsCount);
end
if strcmp(direction,'downlink')
    pathGains = permute(pathGains,[1 2 4 3]);
end

model = struct('profile',profile,'linkDirection',direction,'pathGains',pathGains, ...
    'pathDelays_s',rows(:,1)*options.DelaySpreadSec,'pathPowersLinear',power, ...
    'sampleTimes_s',times,'hasLOS',hasLos,'firstClusterKFactor_dB',firstK, ...
    'losPowerFraction',losPower,'losDoppler_Hz',losDoppler, ...
    'rayDoppler_Hz',rayDoppler,'maxDoppler_Hz',maxDoppler, ...
    'bsClusterDirections',center,'angleSpreads_deg',spreads,'xpr_dB',xprDb);
end

function zenith = foldCdlZenith(zenith)
% TR 38.901 §7.5 Step 7: wrap then fold into [0,180] with azimuth unchanged.
outside = zenith<0 | zenith>180;
wrapped = mod(zenith(outside),360);
zenith(outside) = min(wrapped,360-wrapped);
end

function u = cdlUnitVectors(azimuth,zenith)
u = [sind(zenith).*cosd(azimuth);sind(zenith).*sind(azimuth);cosd(zenith)];
end

function [fTheta,fPhi] = cdlUeField(azimuth,zenith,direction,localUp)
eTheta = [cosd(zenith).*cosd(azimuth);cosd(zenith).*sind(azimuth);-sind(zenith)];
ePhi = [-sind(azimuth);cosd(azimuth);zeros(size(azimuth))];
localTheta = direction.*(localUp.'*direction)-localUp;
normTheta = sqrt(sum(localTheta.^2,1));
pole = normTheta<1e-12;
localTheta(:,pole) = eTheta(:,pole);
normTheta(pole) = 1;
localTheta = localTheta./normTheta;
fTheta = sum(localTheta.*eTheta,1);
fPhi = sum(localTheta.*ePhi,1);
end


function [rows,spreads,xprDb,hasLos] = cdlProfileData(profile)
% CDL parameters from 3GPP TR 38.901 V18.0.0, Tables 7.7.1-1 to 7.7.1-5.
% Columns: normalized delay, power (dB), DL AoD, AoA, ZoD, ZoA (degrees).
% In CDL-D/E, the first two rows are the specular and diffuse components
% of the same zero-delay cluster. The caller combines their coefficients.
% Spreads are [cASD cASA cZSD cZSA], with XPR given separately in dB.
profile = upper(char(profile));
hasLos = any(strcmp(profile,{'CDL-D','CDL-E'}));
switch profile
    case 'CDL-A'
        rows = [
             0.0000  -13.4  -178.1    51.3    50.2   125.4
             0.3819      0    -4.2  -152.7    93.2    91.3
             0.4025   -2.2    -4.2  -152.7    93.2    91.3
             0.5868     -4    -4.2  -152.7    93.2    91.3
             0.4610     -6    90.2    76.6     122      94
             0.5375   -8.2    90.2    76.6     122      94
             0.6708   -9.9    90.2    76.6     122      94
             0.5750  -10.5   121.5    -1.8   150.2    47.1
             0.7618   -7.5   -81.7   -41.9    55.2      56
             1.5375  -15.9   158.4    94.2    26.4    30.1
             1.8978   -6.6     -83    51.9   126.4    58.8
             2.2242  -16.7   134.8  -115.9   171.6      26
             2.1718  -12.4    -153    26.6   151.4    49.2
             2.4942  -15.2    -172    76.6   157.2   143.1
             2.5119  -10.8  -129.9      -7    47.2   117.4
             3.0582  -11.3    -136     -23    40.4   122.7
             4.0810  -12.7   165.4   -47.2    43.3   123.2
             4.4579  -16.2   148.4   110.4   161.8    32.6
             4.5695  -18.3   132.7   144.5    10.8    27.2
             4.7966  -18.9  -118.6   155.3    16.7    15.2
             5.0066  -16.6  -154.1     102   171.7     146
             5.3043  -19.9   126.5  -151.8    22.7   150.7
             9.6586  -29.7   -56.2    55.2   144.9   156.1];
        spreads = [5 11 3 3];
        xprDb = 10;
    case 'CDL-B'
        rows = [
             0.0000      0     9.3  -173.3   105.8    78.9
             0.1072   -2.2     9.3  -173.3   105.8    78.9
             0.2155     -4     9.3  -173.3   105.8    78.9
             0.2095   -3.2   -34.1   125.5   115.3    63.3
             0.2870   -9.8   -65.4   -88.0   119.3    59.9
             0.2986   -1.2   -11.4   155.1   103.2    67.5
             0.3752   -3.4   -11.4   155.1   103.2    67.5
             0.5055   -5.2   -11.4   155.1   103.2    67.5
             0.3681   -7.6   -67.2   -89.8   118.2    82.6
             0.3697     -3    52.5   132.1   102.0    66.3
             0.5700   -8.9     -72   -83.6   100.4    61.6
             0.5283     -9    74.3    95.3    98.3    58.0
             1.1021   -4.8   -52.2   103.7   103.4    78.2
             1.2756   -5.7   -50.5   -87.8   102.5    82.0
             1.5474   -7.5    61.4   -92.5   101.4    62.4
             1.7842   -1.9    30.6  -139.1   103.0    78.0
             2.0169   -7.6   -72.5   -90.6   100.0    60.9
             2.8294  -12.2   -90.6    58.6   115.2    82.9
             3.0219   -9.8   -77.6   -79.0   100.5    60.8
             3.6187  -11.4   -82.6    65.8   119.6    57.3
             4.1067  -14.9  -103.6    52.7   118.7    59.9
             4.2790   -9.2    75.6    88.7   117.8    60.1
             4.7834  -11.3   -77.6   -60.4   115.7    62.3];
        spreads = [10 22 3 7];
        xprDb = 8;
    case 'CDL-C'
        rows = [
                  0   -4.4   -46.6    -101    97.2    87.6
             0.2099   -1.2   -22.8     120    98.6    72.1
             0.2219   -3.5   -22.8     120    98.6    72.1
             0.2329   -5.2   -22.8     120    98.6    72.1
             0.2176   -2.5   -40.7  -127.5   100.6    70.1
             0.6366      0     0.3   170.4    99.2    75.3
             0.6448   -2.2     0.3   170.4    99.2    75.3
             0.6560   -3.9     0.3   170.4    99.2    75.3
             0.6584   -7.4    73.1    55.4   105.2    67.4
             0.7935   -7.1   -64.5    66.5    95.3    63.8
             0.8213  -10.7    80.2   -48.1   106.1    71.4
             0.9336  -11.1   -97.1    46.9    93.5    60.5
             1.2285   -5.1   -55.3    68.1   103.7    90.6
             1.3083   -6.8   -64.3   -68.7   104.2    60.1
             2.1704   -8.7   -78.5    81.5    93.0    61.0
             2.7105  -13.2   102.7    30.7   104.2   100.7
             4.2589  -13.9    99.2   -16.4    94.9    62.3
             4.6003  -13.9    88.8     3.8    93.1    66.7
             5.4902  -15.8  -101.9   -13.7    92.2    52.9
             5.6077  -17.1    92.2     9.7   106.7    61.8
             6.3065    -16    93.3     5.6    93.0    51.9
             6.6374  -15.7   106.6     0.7    92.9    61.7
             7.0427  -21.6   119.5   -21.9   105.2      58
             8.6523  -22.8  -123.8    33.6   107.8      57];
        spreads = [2 15 3 7];
        xprDb = 7;
    case 'CDL-D'
        rows = [
                  0   -0.2       0    -180    98.5    81.5
                  0  -13.5       0    -180    98.5    81.5
              0.035  -18.8    89.2    89.2    85.5    86.9
              0.612    -21    89.2    89.2    85.5    86.9
              1.363  -22.8    89.2    89.2    85.5    86.9
              1.405  -17.9      13     163    97.5    79.4
              1.804  -20.1      13     163    97.5    79.4
              2.596  -21.9      13     163    97.5    79.4
              1.775  -22.9    34.6    -137    98.5    78.2
              4.042  -27.8   -64.5    74.5    88.4    73.6
              7.937  -23.6   -32.9   127.7    91.3    78.3
              9.424  -24.8    52.6  -119.6   103.8      87
              9.708  -30.0  -132.1    -9.1    80.3    70.6
             12.525  -27.7    77.2   -83.8    86.5    72.9];
        spreads = [5 8 3 3];
        xprDb = 11;
    case 'CDL-E'
        rows = [
              0.000  -0.03       0    -180    99.6    80.4
              0.000 -22.03       0    -180    99.6    80.4
             0.5133  -15.8    57.5    18.2   104.2    80.4
             0.5440  -18.1    57.5    18.2   104.2    80.4
             0.5630  -19.8    57.5    18.2   104.2    80.4
             0.5440  -22.9   -20.1   101.8    99.4    80.8
             0.7112  -22.4    16.2   112.9   100.8    86.3
             1.9092  -18.6     9.3  -155.5    98.8    82.7
             1.9293  -20.8     9.3  -155.5    98.8    82.7
             1.9589  -22.6     9.3  -155.5    98.8    82.7
             2.6426  -22.3      19  -143.3   100.8    82.9
             3.7136  -25.6    32.7   -94.7    96.4      88
             5.4524  -20.2     0.5     147    98.9      81
            12.0034  -29.8    55.9   -36.2    95.6    88.6
            20.6419  -29.2    57.6     -26   104.6    78.3];
        spreads = [5 11 3 7];
        xprDb = 8;
    otherwise
        error('UL256:NativeCDLProfile','Use CDL-A, CDL-B, CDL-C, CDL-D or CDL-E.');
end
end
