function cfg = configFixedArray()
%CONFIGFIXEDARRAY 固定平面阵列的上行链路配置。
% 修改本文件后运行 simulateFixedArray；无需 Communications/Statistics Toolbox。
root = fileparts(mfilename('fullpath'));

% 阵列与部署；阵面本地 x 水平、y 法向、z 竖直。
cfg.Nh = 16; cfg.Nv = 8;
cfg.dh_lambda = 0.5; cfg.dv_lambda = 3;
cfg.fc_Hz = 6.7e9; cfg.c0 = 299792458;
cfg.K = 1; % 单天线 UE 数。
cfg.bs.position_xy_m = [0 0]; cfg.bs.ground_height_m = 0;
cfg.bs.height_m = 20; cfg.bs.azimuth_deg = 0; cfg.bs.downtilt_deg = 10;
cfg.users.mode = 'random';
cfg.users.positions_m = []; % specified 模式需提供 K x 3 全局绝对坐标（m）。
cfg.users.radius_range_m = [10 150];
cfg.users.azimuth_offsets_deg = [-60 60];
cfg.users.height_m = 1.5; cfg.users.sampling = 'area';
cfg.users.angular_distance_m = 100;
cfg.user_az_range = [-60 60]; cfg.user_el_range = [-10 10];

% OFDM：24 RB、288 个实际子载波、8.64 MHz；full 为全部频点。
cfg.ofdm.enabled = true;
cfg.ofdm.resource_mode = 'rb'; cfg.ofdm.num_rb = 24;
cfg.ofdm.subcarriers_per_rb = 12; cfg.ofdm.scs_Hz = 30e3;
cfg.ofdm.nfft = 512; cfg.ofdm.rb_start_index = [];
cfg.ofdm.rb_allocation = []; % 所有 UE 共享全部 RB；指定用户分配可用 K 个 cell。
cfg.ofdm.active_indices = -144:143;
cfg.ofdm.frequency_mode = 'full'; cfg.ofdm.num_frequency_samples = 24;
cfg.ofdm.frequency_dependent_steering = true;
cfg.ofdm.cp_mode = 'ratio'; cfg.ofdm.cp_ratio = 0.12;
cfg.ofdm.cp_samples = 61; cfg.ofdm.include_cp_overhead = true;

% 每 UE 总功率 23 dBm；改变 SCS 时也应重算单子载波噪声。
Ptx_dBm = 23;
cfg.noiseVar = 1.380649e-23*290*10^(7/10)*cfg.ofdm.scs_Hz;
cfg.power.snr_reference = 'total_over_tone_noise';
cfg.snr_dB = 10*log10(10^((Ptx_dBm-30)/10)/cfg.noiseVar);
cfg.snr_dB_list = []; % 设置为非空向量可扫描发射功率参数。
cfg.pathloss.enabled = true; cfg.pathloss.exponent = 2.2;
cfg.pathloss.reference_distance_m = 1; cfg.pathloss.distance_floor_m = 1;
cfg.pathloss.reference_loss_dB = 20*log10(4*pi*cfg.fc_Hz/cfg.c0)+35;
cfg.receiver = 'mmse-irc'; % 可改为 zf。
cfg.irc.external_covariance = []; % 接收端口固定的外部干扰协方差（可选）。
cfg.throughput_mode = 'shannon';

% 正式输入接口。没有正式输入时明确报错，绝不自动换成合成数据。
cfg.pattern.file = fullfile(root,'HW_Airview_Fre6700M_E06_1to6_12.7dB.txt');
cfg.pattern.source = 'measured'; % or virtual
cfg.pattern.virtual_hpbw_deg = [120 90]; % azimuth and elevation half-power widths
cfg.pattern.native_to_body = eye(3); cfg.pattern.alignment = 'peak';
cfg.pattern.use_phase = false; cfg.front_only = true;
cfg.cdl.profile_mat = ''; % 填含 rawCDL 或 result.rawCDL 的 MAT 即可跳过 Julia。
cfg.cdl.generator = 'generateCDL'; % 或填本地同接口的真实生成函数名称。
cfg.cdl.extra_args = {'DelaySpread',100};
cfg.cdl.eoa_convention = 'elevation'; cfg.cdl.bs_angle_side = 'arrival';
cfg.cdl.center_mode = 'power_mean'; cfg.cdl.reference_az_el_deg = [];
cfg.cdl.frequency_offset_Hz = 0;

% 默认执行三次单用户 Monte Carlo；可设置 K、MC 和随机种子。
cfg.MC = 3; cfg.seed = 20260928; cfg.progress_every = 1;
cfg.output_dir = fullfile(root,'results');
cfg.save_results = true;
cfg.output.save_per_tone = false;
cfg.output.save_realizations = false;
cfg.check.norm_tol = 1e-10;
cfg.check.zf_rcond_switch = 1e-8; cfg.check.svd_rel_tol = 0;
end
