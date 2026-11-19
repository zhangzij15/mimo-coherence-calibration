# 4T256R NR 上行相干校准基础 MATLAB 版

工程默认使用独立信道、DM-RS 和传输块尺寸实现。支持 **CDL-A、CDL-B、CDL-C、CDL-D 和 CDL-E**。用 MATLAB R2024b 运行，无需安装 5G Toolbox、Communications Toolbox、Statistics and Machine Learning Toolbox 或 Parallel Computing Toolbox。仿真使用基础 MATLAB 自带的 JVM 生成 SHA-256 缓存签名。

## 本地输入

从 Git 检出后，使用现有工作目录或交接 ZIP 中的 `awgn_bler10_t1.csv`，放在仓库根目录，与 `ul4tLinkMapping.m` 同级。这是链路的 EESM 查表输入。该表的原始来源和分发权限待核实，Git 将它作为本地输入忽略。将同一份表按原字节复制到各机器，避免改变查表数据和缓存签名。

当前交接版的 SHA-256 为 `8e93265d8b192e1c3638c8d34947a0d3e6d7dcf86c98b7ec6eeaafde8fb1c834`。单独使用 `generateCdlChannel.m` 时只需复制这个函数文件。

## 运行

在 MATLAB 中切换到仓库根目录。下面的命令执行一次 50 ms 短龄期冒烟测试，使用入选 NR 场景、−20 dB 和一个随机种子。

~~~matlab
scene = {'NR_MODULE_VERTICAL_15_COMPACT_R1'};
opts = struct('stage','screen','scenarioNames',{scene},'screenMC',1, ...
    'channelProfile','CDL-C', ...
    'outputDir',fullfile(pwd,'results_base_smoke'));
report = runNrStrict256(opts);
assert(strcmp(report.status,'complete') && report.screen.completed_mc==1);
~~~

用以下命令复验入选场景。默认筛选用 8 个 MC，独立复验用 24 个 MC 和 −24 至 −12 dB 的七个 SNR 点，另对 −20 dB 运行匹配零漂移对照。

~~~matlab
scene = {'NR_MODULE_VERTICAL_15_COMPACT_R1'};
opts = struct('scenarioNames',{scene}, ...
    'channelProfile','CDL-C', ...
    'outputDir',fullfile(pwd,'results_base_multicdl_c'));
report = runNrStrict256(opts);
assert(strcmp(report.status,'complete'));
addpath(fullfile(pwd,'analysis'));
plotNrStrictThp(fullfile(opts.outputDir,'confirm_thp.csv'),opts.outputDir);
~~~

将 channelProfile 改为 'CDL-A'、'CDL-B'、'CDL-D' 或 'CDL-E'，即可把同一场景切到对应信道。省略 outputDir 时，各信道自动写入相应目录。保持输出目录和配置即可续跑。省略 scenarioNames 可重新筛选全部 15 个 NR 场景。装有 5G Toolbox 时，可设 opts.backend = '5g-toolbox' 运行参考信道及参考 DM-RS。参考结果请使用独立输出目录。

建议把 `outputDir` 设为 `results_` 开头的目录。Git 忽略仿真输出、检查点和绘图产物。移到新电脑后使用新的输出目录，原机器的检查点包含绝对路径。Agent 的接手步骤及结果文件说明见 [AGENTS.md](AGENTS.md)。

## 在其他项目中调用信道生成器

只复制 [generateCdlChannel.m](generateCdlChannel.m) 到目标代码的 MATLAB 路径。五套参数表、极化计算、射线耦合和时间演化都在这一个文件内。

~~~matlab
link = generateCdlChannel("CDL-D", ...
    CarrierFrequencyHz=6.7e9, DelaySpreadSec=100e-9, ...
    SampleTimesSec=(0:20)*0.5e-3, Seed=42, ...
    UserPositionsLambda=[0 0 0; 0 0 0.22], ...
    BaseStationPositionsLambda=[0 0 0; 0 0.5 0], ...
    SpeedMps=0.3/3.6, TravelAnglesDeg=[0 0]);
~~~

pathGains 的维度是 **时刻 × 抽头 × 发射端口 × 接收端口**，默认上行。设置 LinkDirection="downlink" 可获得同一传播实现的上下行互易方向。delay spread、抽头时延及采样时间均按秒输入或输出。两端阵元位置用载频波长作单位，列顺序为 x、y、z。Seed 控制内部随机流。D/E 采用标准表格给定的 LOS 功率，LosPhaseRad 可设置零时刻参考相位。

下例将第一个时刻的抽头系数合成为指定基带频偏下的频响，H 的维度为接收端口 × 发射端口 × 子载波。

~~~matlab
fHz = (-24:23)*30e3;
g = reshape(link.pathGains(1,:,:,:), numel(link.pathDelays_s), []);
delayPhase = exp(-1i*2*pi*link.pathDelays_s*fHz);
nTx = size(link.pathGains,3);
nRx = size(link.pathGains,4);
H = permute(reshape(double(g).'*delayPhase,nTx,nRx,[]),[2 1 3]);
~~~

## 实现范围

| 模块 | 基础 MATLAB 实现 |
|---|---|
| CDL 信道 | TR 38.901 表 7.7.1-1 至 7.7.1-5。A/B/C 为 23/23/24 个 NLOS 抽头，D/E 为 13/14 个抽头，首抽头分别使用 13.3/22 dB 的簇内 K 因子。每个簇生成 20 条散射射线，D/E 的首抽头额外包含相干 LOS 分量。采用单极化阵列和窄带空间转向，项目适配层另外应用 gNB 虚拟方向图。 |
| PUSCH DM-RS | 采用 CP-OFDM、Mapping B、Type 1 前置双符号和 Rank 1–4。按绝对 CRB 偏移生成连续的 288 RB 序列。 |
| TBS | 按 TS 38.214 计算单码字传输块尺寸，覆盖 1–4 层。现有 EESM 与 AWGN BLER 门限保持同一统计口径。 |
| SRS 与校准 | 调度器按 SRS 实际时刻获取信道快照。保留 NR 码本、DM-RS 信道估计、白噪声下的单用户 MMSE-IRC、RF 漂移、反馈补偿、调度重选及 Oracle 对照。 |
| 结果 | 生成逐 MC 数据、配对置信区间、候选状态、端口误差审计和四条短 SRS 龄期 THP–SNR 曲线。 |

当前任务的物理层口径是 **SRS 时刻信道快照与 EESM 吞吐率估算**。因此这版的验证范围是相干校准链路仿真。逐采样 SRS 波形估计、PUSCH 软解调、LDPC 译码与按 CRC 统计 BLER 属于后续波形级工作。

独立 CDL 依据公开参数重新生成衰落。相同种子保证本后端内部可重复，跨后端按信道统计和链路结论比较。THP 仍以十个短龄期时隙的有效比特除以 50 ms 计算，主 THP 未扣除校准开销。

## 验证

~~~matlab
results = runtests(pwd);
assert(~any([results.Failed]));
~~~

在没有 5G Toolbox 的环境，参考工具箱的对照用例会按假设跳过。基础链路和独立实现的测试仍会执行。

若只准备验证独立 CDL 函数，可以直接运行 `runtests('test_generate_cdl_channel.m')`。完整链路测试需要本地输入表。

## 可选的 5G Toolbox 对照

装有 5G Toolbox 时，可从仓库根目录使用以下命令在本地生成路径系数、随机信道统计量和基站空间相关收敛结果。

~~~matlab
addpath(fullfile(pwd,'analysis'));
compareCdlMatchedPhases;
compareCdlToToolbox(80,30);
checkCdlABsCovariance;
~~~

输出写入 `analysis/cdl_toolbox_comparison/`。固定射线对照显式共享耦合和初相。给两个后端传相同 `Seed` 时，各自的随机实现可重复，瞬时复系数会不同。

- [CDL 信道参数](https://www.etsi.org/deliver/etsi_tr/138900_138999/138901/18.00.00_60/tr_138901v180000p.pdf#page=73)
- [NR DM-RS 与码本](https://www.etsi.org/deliver/etsi_ts/138200_138299/138211/17.06.00_60/ts_138211v170600p.pdf)
- [TBS 计算](https://www.etsi.org/deliver/etsi_ts/138200_138299/138214/17.08.00_60/ts_138214v170800p.pdf)
