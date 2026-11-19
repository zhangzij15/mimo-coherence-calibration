function varargout = linkScene(action,varargin)
%LINKSCENE 基站部署与 UE 三维位置的公共接口，长度单位均为米。
% scene = linkScene('prepare',cfg) -> 基站位置、阵列姿态和坐标说明。
% users = linkScene('users',cfg) -> 本次 MC 的 UE 位置、距离、方向和路径损耗。
% 这里只描述 BS 到 UE 的 nominal 几何方向，不替代 CDL 的簇内多径角度。
switch action
    case 'prepare'
        varargout{1} = prepareScene(varargin{1});
    case 'users'
        varargout{1} = generateSceneUsers(varargin{1});
    otherwise
        error('Unknown linkScene action: %s',action);
end
end

function scene = prepareScene(cfg)
% 阵列局部 +y 为正面，+x 水平，+z 向上；正 downtilt 令局部 +y 朝地下俯。
az = cfg.bs.azimuth_deg;
down = cfg.bs.downtilt_deg;
Rz = [cosd(az) sind(az) 0;-sind(az) cosd(az) 0;0 0 1];
Rx = [1 0 0;0 cosd(down) sind(down);0 -sind(down) cosd(down)];
scene.bs.position_m = [cfg.bs.position_xy_m(:).',cfg.bs.ground_height_m+cfg.bs.height_m];
scene.bs.R_global_from_array = Rz*Rx;
scene.bs.boresight_global = scene.bs.R_global_from_array*[0;1;0];
scene.bs.ground_height_m = cfg.bs.ground_height_m;
scene.bs.height_m = cfg.bs.height_m;
scene.bs.azimuth_deg = az;
scene.bs.downtilt_deg = down;
% 顶层别名便于只需要位姿矩阵的调用者使用；与 bs 内字段完全相同。
scene.position_m = scene.bs.position_m;
scene.R_global_from_array = scene.bs.R_global_from_array;
scene.metadata.global_axes = '+x right, +y azimuth zero, +z up';
scene.metadata.array_axes = '+x horizontal, +y boresight, +z vertical';
scene.metadata.azimuth_convention = 'azimuth=atan2(x,y); positive from +y toward +x';
scene.metadata.position_unit = 'm';
scene.metadata.rotation_convention = 'R_global_from_array=Rz(-azimuth)*Rx(-downtilt)';
scene.metadata.user_mode = cfg.users.mode;
end

function users = generateSceneUsers(cfg)
scene = prepareScene(cfg);
bs = scene.bs;
R = bs.R_global_from_array;
K = cfg.K;
switch cfg.users.mode
    case 'specified'
        % 固定用户坐标每次 MC 复用；小尺度 CDL 系数仍可独立重新抽样。
        positions = cfg.users.positions_m;
    case 'random'
        range = cfg.users.radius_range_m;
        if strcmp(cfg.users.sampling,'area')
            radius = sqrt(range(1)^2+(range(2)^2-range(1)^2)*rand(K,1));
        else
            radius = range(1)+(range(2)-range(1))*rand(K,1);
        end
        azimuth = cfg.bs.azimuth_deg+cfg.users.azimuth_offsets_deg(1)+diff(cfg.users.azimuth_offsets_deg)*rand(K,1);
        x = bs.position_m(1)+radius.*sind(azimuth);
        y = bs.position_m(2)+radius.*cosd(azimuth);
        z = (cfg.bs.ground_height_m+cfg.users.height_m)*ones(K,1);
        positions = [x,y,z];
    case 'angles'
        % 旧版角度模式严格在阵列本地坐标抽样；不会再减一次基站下倾角。
        azimuthArray = cfg.user_az_range(1)+diff(cfg.user_az_range)*rand(K,1);
        elevationArray = cfg.user_el_range(1)+diff(cfg.user_el_range)*rand(K,1);
        directionArray = [cosd(elevationArray).*sind(azimuthArray),cosd(elevationArray).*cosd(azimuthArray),sind(elevationArray)].';
        directionGlobal = R*directionArray;
        positions = bsxfun(@plus,cfg.users.angular_distance_m*directionGlobal.',bs.position_m);
    otherwise
        error('Unknown users.mode: %s',cfg.users.mode);
end
delta = bsxfun(@minus,positions,bs.position_m);
distance3d = sqrt(sum(delta.^2,2));
assert(all(isfinite(distance3d) & distance3d>0),'A UE cannot coincide with the BS position.');
if ~strcmp(cfg.users.mode,'angles')
    directionGlobal = bsxfun(@rdivide,delta,distance3d).';
    directionArray = R.'*directionGlobal;
end
users.positions_m = positions;
users.direction_global = directionGlobal;
users.direction_array = directionArray;
users.azimuth_deg = atan2d(directionGlobal(1,:),directionGlobal(2,:));
users.elevation_deg = atan2d(directionGlobal(3,:),sqrt(sum(directionGlobal(1:2,:).^2,1)));
users.azimuth_array_deg = atan2d(directionArray(1,:),directionArray(2,:));
users.elevation_array_deg = atan2d(directionArray(3,:),sqrt(sum(directionArray(1:2,:).^2,1)));
users.distance_2d_m = sqrt(sum(delta(:,1:2).^2,2));
users.distance_3d_m = distance3d;
users.bs = bs;
users.mode = cfg.users.mode;
users.coordinate_note = scene.metadata;

% 对数距离路径损耗是独立可选的研究模型，不声称为完整 3GPP 场景模型。
% d0 以内不外推：同时使用参考距离和显式距离下限裁切，避免虚增近距离增益。
d0 = cfg.pathloss.reference_distance_m;
if isempty(cfg.pathloss.reference_loss_dB)
    loss0 = 20*log10(4*pi*cfg.fc_Hz*d0/cfg.c0);
else
    loss0 = cfg.pathloss.reference_loss_dB;
end
effectiveDistance = max(distance3d,max(d0,cfg.pathloss.distance_floor_m));
modelLoss = loss0+10*cfg.pathloss.exponent*log10(effectiveDistance/d0);
users.pathloss_model_dB = modelLoss;
users.pathloss_effective_distance_m = effectiveDistance;
users.pathloss_reference_loss_dB = loss0;
users.pathloss_enabled = cfg.pathloss.enabled;
users.pathloss_dB = zeros(K,1);
users.pathloss_power_gain = ones(K,1);
if cfg.pathloss.enabled
    users.pathloss_dB = modelLoss;
    users.pathloss_power_gain = 10.^(-modelLoss/10);
end
assert(all(isfinite(users.pathloss_power_gain) & users.pathloss_power_gain>0),'Pathloss produced a zero or invalid power gain.');
% 方向和路径损耗都由同一 nominal 链路生成。
end
