function varargout = fixedArrayChannel(action,varargin)
% 固定平面阵列信道。输入 cfg 由 prepareFixedArrayConfig 准备。
switch action
    case 'prepare'
        cfg = varargin{1};
        varargout{1} = loadAndNormalizePattern(cfg);
        [varargout{2},varargout{3}] = loadCDLProfile(cfg);
    case 'realization'
        varargout{1} = buildChannelRealization(varargin{:});
    case 'pattern'
        varargout{1} = loadAndNormalizePattern(varargin{:});
    case 'pattern-power'
        varargout{1} = patternPower(varargin{:});
    case 'orient-clusters'
        [varargout{1:nargout}] = orientClusters(varargin{:});
    case 'response'
        varargout{1} = buildFixedResponse(varargin{:});
    otherwise
        error('Unknown fixedArrayChannel action: %s',action);
end
end

function pattern = loadAndNormalizePattern(cfg)
%LOADANDNORMALIZEPATTERN Build the measured or virtual power grid.
phi = 0:359;
theta = 0:180;
[P,T] = meshgrid(phi,theta);
native = [sind(T(:)).*cosd(P(:)),sind(T(:)).*sind(P(:)),cosd(T(:))].';
if strcmp(cfg.pattern.source,'virtual')
    widths = cfg.pattern.virtual_hpbw_deg;
    assert(isnumeric(widths) && isreal(widths) && numel(widths)==2 && ...
        all(isfinite(widths) & widths>0 & widths<=180),'Invalid virtual pattern HPBW.');
    az = mod(90-P+180,360)-180;
    el = 90-T;
    rolloff_dB = 10*log10(2)*((2*az/widths(1)).^2+(2*el/widths(2)).^2);
    G = 10.^(-min(30,rolloff_dB)/10);
    configuredRotation = eye(3);
    actualRotation = configuredRotation;
    headers = {'Virtual BS element power pattern','Azimuth/elevation half-power widths'};
    sourceFile = '';
else
    assert(strcmp(cfg.pattern.source,'measured'),'Pattern source must be measured or virtual.');
    assert(exist(cfg.pattern.file,'file')==2,'Missing measured pattern TXT: %s',cfg.pattern.file);
    [Vgain,headers] = readPattern(cfg.pattern.file);
    assert(isequal(size(Vgain),[181 360]),'Pattern must be a Theta x Phi grid.');
    assert(~any(isnan(Vgain(:))) && ~any(Vgain(:)==Inf),'Invalid measured gain samples.');
    G = 10.^(Vgain/10);
    configuredRotation = cfg.pattern.native_to_body;
    actualRotation = configuredRotation;
    if strcmp(cfg.pattern.alignment,'peak')
        [~,peakIndex] = max(G(:));
        preliminaryPeak = configuredRotation*native(:,peakIndex);
        correction = directionFrame([0;1;0])*directionFrame(preliminaryPeak)';
        actualRotation = correction*configuredRotation;
    end
    sourceFile = cfg.pattern.file;
end
body = actualRotation*native;
if cfg.front_only
    G(reshape(body(2,:)<-1e-12,size(G))) = 0;
end
phi = 0:360;
G = [G,G(:,1)];
rawIntegral = trapz(theta*pi/180,trapz(phi*pi/180,G,2).*sind(theta(:)));
assert(isfinite(rawIntegral) && rawIntegral>0,'Pattern spherical integral must be positive.');
G = G*(4*pi/rawIntegral);
pattern = struct('phi',phi,'theta',theta,'G',G,'native_to_body',actualRotation,'front_only',cfg.front_only);
pattern.configured_native_to_body = configuredRotation;
pattern.alignment = cfg.pattern.alignment;
pattern.source = cfg.pattern.source;
pattern.source_file = sourceFile;
pattern.file_headers = headers;
pattern.is_test_pattern = ~isempty(strfind(upper([headers{1} ' ' headers{2}]),'TEST'));
pattern.rawIntegral = rawIntegral;
pattern.scale = 4*pi/rawIntegral;
pattern.integral = trapz(theta*pi/180,trapz(phi*pi/180,G,2).*sind(theta(:)));
assert(abs(pattern.integral-4*pi)<cfg.check.norm_tol*4*pi,'Pattern normalization failed.');
end

function [gain,headers] = readPattern(filename)
fid = fopen(filename,'r');
assert(fid~=-1,'Cannot open configured pattern TXT: %s',filename);
cleanup = onCleanup(@() fclose(fid));
headers = {fgetl(fid),fgetl(fid)};
assert(ischar(headers{1}) && ischar(headers{2}),'Pattern TXT must contain two header lines.');
data = fscanf(fid,'%f %f %f %f',[4 inf]).';
assert(size(data,1)==360*181,'Expected 360*181 Phi/Theta/gain/phase rows after two headers.');
[P,T] = ndgrid(0:359,0:180);
assert(all(abs(data(:,1)-P(:))<1e-8) && all(abs(data(:,2)-T(:))<1e-8),'Pattern order must be Phi-fast and Theta-slow.');
gain = reshape(data(:,3),360,181).';
end

function [profile,raw] = loadCDLProfile(cfg)
%LOADCDLPROFILE 唯一 Julia/CDL 适配入口；不调用 5G Toolbox。
if ~isempty(cfg.cdl.profile_mat)
    raw = readCDLProfileFile(cfg.cdl.profile_mat);
else
    generator = cfg.cdl.generator;
    assert(~isempty(which(generator)),'CDL generator is not on the MATLAB path: %s',generator);
    modelType = 'NR-CDL-C';
    if isfield(cfg.cdl,'model_type') && ~isempty(cfg.cdl.model_type)
        modelType = cfg.cdl.model_type;
    end
    raw = feval(generator,modelType,'CenterFrequency',cfg.fc_Hz,cfg.cdl.extra_args{:});
end
required = {'pow','taus'};
if strcmp(cfg.cdl.bs_angle_side,'arrival')
    af = 'AoA';
    ef = 'EoA';
else
    af = 'AoD';
    ef = 'EoD';
end
required = [required,{af,ef}];
for i = 1:numel(required)
    assert(isfield(raw,required{i}),['Missing CDL field: ' required{i}]);
    v = raw.(required{i});
    assert(isnumeric(v) && isreal(v) && isvector(v) && all(isfinite(v(:))),['Invalid CDL field: ' required{i}]);
end
p = double(raw.pow(:));
tau = double(raw.taus(:));
az = raw.(af);
ea = raw.(ef);
az = double(az(:));
ea = double(ea(:));
assert(numel(p)==numel(tau) && numel(p)==numel(az) && numel(p)==numel(ea),'CDL field lengths differ.');
assert(all(p>=0) && sum(p)>0 && all(tau>=0),'pow must be nonnegative linear powers, taus must be seconds.');
assert(all(abs(az)<=2*pi+1e-8),'AoA/AoD must be radians, not degrees.');
if strcmp(cfg.cdl.eoa_convention,'zenith')
    assert(all(ea>=-1e-8 & ea<=pi+1e-8),'Zenith angles must lie in [0,pi].');
    el = pi/2-ea;
else
    assert(all(abs(ea)<=pi/2+1e-8),'EoA/EoD is outside elevation range; check whether Julia returns zenith.');
    el = ea;
end
profile.p = p/sum(p);
profile.tau = tau;
profile.U = [cos(el).*cos(az),cos(el).*sin(az),sin(el)].';
profile.numPaths = numel(p);
profile.model = 'NR-CDL-C scalar cluster approximation';
profile.bs_angle_side = cfg.cdl.bs_angle_side;
profile.eoa_convention = cfg.cdl.eoa_convention;
if isfield(raw,'numClusters')
    assert(raw.numClusters==numel(p),'Expected one angle/power entry per cluster.');
end
% KF/XPR/AoD/EoD 保留在 raw 中。此标量 NLoS 首版不构造极化矩阵或子径。
end

function raw = readCDLProfileFile(profileFile)
% 局部作用域内无返回值 load，兼容不支持 s=load(...) 的 MLang。
% 外部 MAT 的变量不会覆盖主信道函数的 cfg/profile 等已有状态。
load(profileFile);
if exist('rawCDL','var')==1
    raw = rawCDL;
elseif exist('result','var')==1 && isstruct(result) && isfield(result,'rawCDL')
    raw = result.rawCDL;
else
    error('MAT must contain rawCDL or result.rawCDL.');
end
end

function [U,info] = rotateCDLAnglesToUser(U0,p,target,cfg)
%ROTATECDLANGLESTOUSER 单个三维刚性旋转保持所有簇之间的夹角。
p = p(:);
target = target(:)/norm(target);
assert(size(U0,1)==3 && size(U0,2)==numel(p),'Profile shape mismatch.');
switch cfg.cdl.center_mode
    case 'power_mean'
        reference = U0*p;
        assert(norm(reference)>1e-8,'CDL power-weighted direction is ambiguous; configure explicit reference.');
        reference = reference/norm(reference);
    case 'explicit'
        a = cfg.cdl.reference_az_el_deg;
        assert(numel(a)==2 && all(isfinite(a)) && abs(a(2))<=90,'Set reference [standard az, elevation].');
        reference = [cosd(a(2))*cosd(a(1));cosd(a(2))*sind(a(1));sind(a(2))];
    otherwise
        error('Unknown CDL center_mode.');
end
% 先用一个正交局部框架移走 profile 中心，再放到目标；不用 az/el 直接相加。
R = directionFrame(target)*directionFrame(reference)';
U = R*U0;
info.R = R;
info.reference = reference;
info.target = target;
end

function F = directionFrame(u)
% 右手框架的第一列指向中心，使用当前参考坐标的 +z 定义其余轴，极点改用 +x。
up = [0;0;1];
if abs(up'*u)>1-1e-10
    up = [1;0;0];
end
side = cross(up,u);
side = side/norm(side);
F = [u,side,cross(u,side)];
end

function G = patternPower(U,pattern)
% U 已在 BS 阵列局部坐标；将来波方向映射到方向图 TXT 的原生坐标。
native = pattern.native_to_body'*U;
phi = mod(atan2d(native(2,:),native(1,:)),360);
theta = acosd(max(-1,min(1,native(3,:))));
G = interp2(pattern.phi,pattern.theta,pattern.G,phi,theta,'linear');
assert(all(isfinite(G)),'Pattern query returned invalid gain.');
G = max(G,0);
if pattern.front_only
    G(U(2,:)<-1e-12) = 0;
end
end

function [u,uGlobal,rotation] = orientClusters(cfg,scene,k,profile)
Rbs = scene.bs.R_global_from_array;
if strcmp(cfg.users.mode,'angles')
    [u,rotation] = rotateCDLAnglesToUser(profile.U,profile.p,scene.direction_array(:,k),cfg);
    uGlobal = Rbs*u;
    rotation.reference_frame = 'array_local_legacy_angles';
    rotation.R_array = rotation.R;
    rotation.R_global = Rbs*rotation.R;
else
    [uGlobal,rotation] = rotateCDLAnglesToUser(profile.U,profile.p,scene.direction_global(:,k),cfg);
    u = Rbs'*uGlobal;
    rotation.reference_frame = 'global_environment';
    rotation.R_global = rotation.R;
    rotation.R_array = Rbs'*rotation.R;
    rotation.R = rotation.R_array;
end
rotation.target_global = scene.direction_global(:,k);
rotation.target_array = scene.direction_array(:,k);
rotation.target = rotation.target_array;
end

function positions = planarPositions(cfg)
lambda = cfg.c0/cfg.fc_Hz;
x = ((0:cfg.Nh-1)-(cfg.Nh-1)/2)*cfg.dh_lambda*lambda;
z = ((0:cfg.Nv-1)-(cfg.Nv-1)/2)*cfg.dv_lambda*lambda;
[X,Z] = ndgrid(x,z);
positions = [X(:),zeros(numel(X),1),Z(:)];
end

function realization = buildChannelRealization(cfg,profile)
%BUILDCHANNELREALIZATION 一次抽样，所有频点和 SNR 共用。
% alpha 是每条簇级路径的复系数；频率选择性在信道构造时由 tau 产生。
% 部署位置先经 BS 的全局姿态变换到阵列局部坐标。
% 假定完美上行定时对齐：每 UE 的公共距离 d/c 已补偿，不加入 CDL 超额时延/CP。
% BS 的全局平移对应的参考载频相位吸收入 alpha；不得重复加入局部空间相位。
scene = linkScene('users',cfg);
realization.userScene = scene;
realization.bs_frame = scene.bs;
realization.coordinate_frame = 'BS array local: +x horizontal, +y front, +z up';
realization.timing_note = 'Perfect uplink timing alignment removes each UE common d/c; tau contains CDL excess delay only.';
realization.users.az = scene.azimuth_array_deg(:);
realization.users.el = scene.elevation_array_deg(:);
realization.users.U = scene.direction_array;
assert(isequal(size(realization.users.U),[3 cfg.K]),'Scene must provide local BS-frame 3 x K directions.');
assert(numel(scene.pathloss_power_gain)==cfg.K && all(isfinite(scene.pathloss_power_gain(:))) && all(scene.pathloss_power_gain(:)>=0),'Invalid scene path-loss power gains.');
realization.paths = cell(cfg.K,1);
earliest = Inf;
latest = -Inf;
activePaths = profile.p(:)>0;
assert(any(activePaths),'At least one positive-power CDL path is required.');
positions = planarPositions(cfg);
for k = 1:cfg.K
    [u,uGlobal,rotation] = orientClusters(cfg,scene,k,profile);
    alpha = sqrt(profile.p(:)/2).*(randn(profile.numPaths,1)+1i*randn(profile.numPaths,1));
    % 路损为功率系数，只在这里以平方根乘入幅度一次；不再按用户/频点归一化。
    alpha = sqrt(scene.pathloss_power_gain(k))*alpha;
    realization.paths{k} = struct('U',u,'U_local',u,'U_global',uGlobal,'alpha',alpha,'tau',profile.tau(:),'rotation',rotation,'pathloss_power_gain',scene.pathloss_power_gain(k));
    % h_n(f) 中的总时延为 tau_l - r_n'*u_l/c0。
    % 检查所有 UE、阵元、非零功率路径的共同 FFT 窗口所需时延跨度。
    % 即使关闭频变 steering，仍保留物理阵列几何时延的保守 CP 约束。
    effectiveDelay = bsxfun(@minus,profile.tau(activePaths).',positions*u(:,activePaths)/cfg.c0);
    earliest = min(earliest,min(effectiveDelay(:)));
    latest = max(latest,max(effectiveDelay(:)));
end
realization.delay_min_s = earliest;
realization.delay_max_s = latest;
realization.delay_span_s = latest-earliest;
realization.array_positions_m = positions;
if cfg.ofdm.enabled && realization.delay_span_s>cfg.ofdm.cp_duration_s+1e-15
    error('fixedArrayChannel:CPTooShort','CP %.6g us is shorter than common FFT-window delay span %.6g us. Increase cfg.ofdm.cp_ratio (ratio mode) or cp_samples (samples mode), or reduce channel delay spread; this model does not simulate ISI/ICI.',cfg.ofdm.cp_duration_s*1e6,realization.delay_span_s*1e6);
end
% 不按 realization 或频点归一化 alpha/H；保留频率选择性深衰落。
end

function H = buildFixedResponse(realization,pattern,cfg)
% H(n,k,f) 为接收阵元 n、单天线 UE k、计算频点 f 的复信道。
F = cfg.ofdm.num_tones;
H = complex(zeros(cfg.N,cfg.K,F));
for k = 1:cfg.K
    path = realization.paths{k};
    projection_s = realization.array_positions_m*path.U/cfg.c0;
    amplitudes = sqrt(patternPower(path.U,pattern)).*path.alpha(:).';
    for f = 1:F
        nu = cfg.ofdm.freq_offsets_Hz(f);
        spatialFrequency = cfg.fc_Hz;
        if cfg.ofdm.frequency_dependent_steering
            spatialFrequency = cfg.fc_Hz+nu;
        end
        spatialPhase = exp(1i*2*pi*spatialFrequency*projection_s);
        delayPhase = exp(-1i*2*pi*nu*path.tau(:));
        H(:,k,f) = bsxfun(@times,spatialPhase,amplitudes)*delayPhase;
    end
end
assert(all(isfinite(H(:))),'Fixed array channel contains nonfinite coefficients.');
end
