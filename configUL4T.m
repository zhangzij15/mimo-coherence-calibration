function cfg = configUL4T(mode)
%CONFIGUL4T Single-user NR CP-OFDM, two or four UE ports and a fixed BS array.
if nargin<1, mode = 'base'; end
mode = validatestring(mode,{'base','local','remote','toolbox'});
if strcmp(mode,'local'), mode = 'base'; end
cfg = configFixedArray();
cfg.Nh = 8;
cfg.Nv = 8;
cfg.dh_lambda = .5;
cfg.dv_lambda = 1.5;
cfg.pattern.file = '';
cfg.pattern.source = 'virtual';
cfg.cdl.backend = 'base-cdl';
cfg.cdl.toolbox_profile = 'CDL-C';
cfg.cdl.toolbox_delay_spread_s = 100e-9;
cfg.pathloss.enabled = false;
cfg.ofdm.scs_Hz = 30e3;
cfg.ofdm.num_rb = 240; % 2880 active subcarriers, 86.4 MHz at 30 kHz SCS
cfg.ofdm.nfft = 5120;
cfg.ofdm.frequency_dependent_steering = false;
cfg.ofdm.cp_ratio = 144/2048; % shortest NR normal CP / useful symbol
cfg.MC = 20;
cfg.seed = 20261113;
cfg.ul.rx_noise_power = 10.^(-(-20:5:-10)/10); % per Rx antenna and complex RE
cfg.ul.num_tx = 2;
cfg.ul.rank = 1;
cfg.ul.ue_layout = 'handheld'; % TR 38.901 Release 19, or 'upa'
cfg.ul.ue_antenna_indices = []; % default [1 7] for 2T, [1 7 3 5] for 4T
cfg.ul.tx_positions_lambda = [0 0 0; 0 0 .5; .5 0 0; .5 0 .5]; % 'upa' only
cfg.ul.ue_speed_kmh = 3;
cfg.ul.channel_aging_enabled = true;
cfg.ul.dmrs_filter_cache_enabled = true; % MC-local reuse for equal layer-power priors
cfg.ul.tdd_pattern = 'DSUDD';
cfg.ul.block_ms = 5;
cfg.ul.duration_ms = 200;
cfg.ul.srs_period_ms = 40; % minimum time between SRS updates at S blocks
cfg.ul.nid = 1;
cfg.ul.nscid = 0;
cfg.ul.prb_start_crb = 0;
cfg.ul.codebook_unit_power = true;
cfg.ul.selection_noise_power = 1; % fixed reference for rank > 1 TPMI selection
cfg.ul.coherence_enabled = true;
cfg.ul.coherence_model = 'uniform'; % 'uniform' or 'drift'
cfg.ul.coherence_dB = 4; % per-port uniform half-width in dB
cfg.ul.coherence_deg = 40; % per-port uniform half-width in degrees
cfg.ul.coherence_drift = struct( ...
    'sigmaPhi0Deg',20,'sigmaA0dB',0.8, ...
    'sigmaDeltaFHz',1.2,'rhoFreqCommon',0.4, ...
    'qPhiDeg2PerSec',18000,'rhoPhiCommon',0.4, ...
    'sigmaADynDb',1.0,'tauA',20e-3,'rhoACommon',0.2);
cfg.ul.normalize_pusch_energy = false; % legacy receive-energy normalization
cfg.ul.hold_tx_power = true; % match the nominal transmit precoder's total power
cfg.ul.mapping = 'eesm'; % 5G-LENA MCS Table 1 calibration
cfg.ul.sem_beta = 1; % scalar or MCS 0:28 vector, linear scaling inside log(1+beta*SINR)
cfg.ul.show_figure = true;
cfg.ul.save_figure = true;
cfg.ul.thp_scope = 'frame'; % 'frame' includes TDD duty, 'uplink' reports U-block rate
cfg.ul.cases = struct('precoder',{'svd','svd','nr','nr'}, ...
    'coherence_error',{false,true,false,true}, ...
    'sets',{{'full'},{'full'},{'full'},{'full'}});
cfg.ul.calibration = struct('num_re',16,'layout','frequency', ...
    'prefix_symbols',[],'pilot_power',1,'innovation_variance',0.25);
root = fileparts(mfilename('fullpath'));
cfg.output_dir = fullfile(root,['results_' char(mode)]);
if any(strcmp(mode,{'remote','toolbox'}))
    cfg.cdl.backend = '5g-toolbox';
    cfg.ul.show_figure = false;
end
if any(strcmp(mode,{'base','toolbox'}))
    cfg.pattern.virtual_hpbw_deg = [120 90];
    cfg.ul.num_tx = 4;
    cfg.ul.coherence_model = 'drift';
    cfg.ul.calibration.num_re = 128;
    cfg.ul.calibration.prefix_symbols = 1;
    cfg.ul.save_figure = false;
    cfg.output_dir = fullfile(root,['results_matlab_' mode]);
end
if strcmp(mode,'base')
    cfg.Nh = 16;
    cfg.Nv = 16;
    cfg.ofdm.num_rb = 288;
end
end
