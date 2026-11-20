# 4T256R NR 上行相干校准工程交接

接手后先阅读本文件和 [README.md](README.md)。主入口是 `runNrStrict256.m`。用户如果指定了场景、信道或种子，以用户指定值为准，并在结果里记录实际配置。

## 环境与首次运行

在 MATLAB 中将当前目录切到本文件所在的仓库根目录。已验证环境为 MATLAB R2024b，启动时保留 MATLAB 自带的 JVM，用于缓存签名。默认后端为 `base-cdl`，只用基础 MATLAB。新检出需从现有工作目录或交接 ZIP 复制 `awgn_bler10_t1.csv` 到仓库根目录，来源、文件位置和校验值见 README 的本地输入一节。工程有同名的历史版本，切换时只把当前仓库和它的 `analysis` 加入 MATLAB 路径。

先运行以下单种子、单 SNR 的 50 ms 冒烟测试。命令适用于 MATLAB 命令窗口或 `matlab -batch`。

```matlab
assert(isfile('awgn_bler10_t1.csv') && isfile('runNrStrict256.m'));
disp(version('-release'));
disp(which('runNrStrict256'));
opts = struct('stage','screen','scenarioNames',{{'NR_MODULE_VERTICAL_15_COMPACT_R1'}}, ...
    'backend','base-cdl','channelProfile','CDL-C','screenMC',1, ...
    'outputDir',fullfile(pwd,'results_handoff_smoke'));
report = runNrStrict256(opts);
assert(strcmp(report.status,'complete'));
assert(isscalar(report.screen) && strcmp(report.screen.status,'complete') && ...
    report.screen.completed_mc==1);
```

用 `results_handoff_smoke/run.log` 看进度，用 `screen_status.csv` 检查完成数和错误。`results = runtests(pwd)` 可运行工程自测。检查 `results.Failed`，其中依赖 `nrCDLChannel`、`nrPUSCHDMRS` 或 `nrTBS` 的工具箱对照会在缺少对应函数时按测试假设跳过。

## 复验已有场景并绘图

使用与冒烟测试不同的输出目录。以下命令筛选用 8 个 MC，使用独立的 24 个 MC 在 −24 至 −12 dB 的七档 SNR 上复验，并在 −20 dB 加做同配置零漂移对照。

```matlab
opts = struct('scenarioNames',{{'NR_MODULE_VERTICAL_15_COMPACT_R1'}}, ...
    'backend','base-cdl','channelProfile','CDL-C', ...
    'outputDir',fullfile(pwd,'results_handoff_cdl_c'));
report = runNrStrict256(opts);
assert(strcmp(report.status,'complete'));
addpath(fullfile(pwd,'analysis'));
files = plotNrStrictThp(fullfile(opts.outputDir,'confirm_thp.csv'),opts.outputDir);
```

- `confirm_status.csv` 给出完成情况和错误。`confirm_summary.csv` 汇总损失、收益及配对置信区间，`confirm_mc.csv` 保留逐种子数据。`confirm_thp.csv` 提供四条曲线的数据，绘图输出同名的 `nr_strict_thp_vs_snr.png` 和 `.pdf`。
- 迁移到新电脑时从新输出目录启动，旧检查点记录了原机器上的绝对路径。本机续跑时使用相同的输出目录和配置，已有筛选可用相同参数和 `stage='confirm'` 继续复验。换源代码、信道、后端或种子时用新的输出目录，便于追溯。
- `channelProfile` 可以切换为 `CDL-A`、`CDL-B`、`CDL-D` 或 `CDL-E`。省略 `scenarioNames` 可重新筛选所有 15 个 NR 场景。做正式结论时，保留筛选和复验的独立种子，并在同一 MC 内对四组方案做配对。

## 复用与核对

- `generateCdlChannel.m` 是可单独复制的五模型信道生成器，包含最新的天顶角折回修复。`pathGains` 的维度为时刻 × 抽头 × 发射端口 × 接收端口。时间用秒，阵元位置用载频波长。示例见 README。
- `simulateUL4T.m` 包含链路、NR 码本、接收机和校准逻辑。`nrStrictScenarios.m` 定义候选场景。`ul4tBaseDMRS.m`、`ul4tBaseTBS.m` 和 `ul4tNativeCDL.m` 是基础后端主要模块。
- 复现实验的口径是单用户、4T256R、6.7 GHz、288 RB、30 kHz SCS、MMSE-IRC 和 NR 码本。主结果累计实际 SRS age ≤20 ms 的有效比特，并除以整个 50 ms 仿真时长。RF 误差按端口对检查，界限为 40° 和 4 dB。四条主曲线依次是理想相干、非理想相干、反馈校准和调度校准，调度校准在代码中称 `Reselect`。链路使用 SRS 时刻的信道快照和 EESM 吞吐率估算。
- 安装 5G Toolbox 时，可用 `opts.backend='5g-toolbox'` 跑可选参考后端，或按 README 的命令运行 `analysis` 下的三个信道对照脚本。输出默认写入本地 `analysis/cdl_toolbox_comparison/`。相同 `Seed` 可重现各自后端的信道。逐复系数对照时需要对齐射线耦合和初相。

## Git 提交与同步

- 同步远端为 `origin`，对应 GitHub 私有仓库 [zhangzij15/mimo-coherence-calibration](https://github.com/zhangzij15/mimo-coherence-calibration)。推送前用 `git remote get-url --push origin` 核对目标。
- 完成修改后运行相关测试。测试通过且用户明确批准本次提交时，先提交本次改动，再推送当前分支到上述私有仓库，并核对远端提交号。用户对本次提交的批准同时授权这次推送。用户当次明确要求仅在本地提交时，按当次要求执行。其余修改保留在本地工作区，等待用户批准。
- 用 `git diff --cached --name-only` 核对暂存文件，并在推送前检查将同步的提交。提交范围为本次改动的源码、测试、绘图和对照脚本及使用说明。`awgn_bler10_t1.csv` 和仿真结果始终保存在本机，包括结果 CSV、图、日志和 `.mat` 检查点。遵守 `.gitignore` 的排除规则。
- 正式仿真将 `git rev-parse --short HEAD` 返回的版本号与结果一起保存在本地，方便另一台机器复现。
