function scenarios = nrStrictScenarios()
%NRSTRICTSCENARIOS NR-codebook loads with a globally bounded RF response.
% The two PA module groups are fixed before channel and seed draws.
moduleByColumn = [-.5;-.5;.5;.5];
moduleByRow = [-.5;.5;-.5;.5];
base = struct('name','NR_MODULE_20_R1','spacing_lambda',.35,'speed_kmh',.3, ...
    'rank',1,'tau_ms',.5,'phase20_deg',39.8,'amp20_dB',3.95, ...
    'direction',moduleByColumn,'initial_offset_deg',20,'initial_jitter_deg',5, ...
    'sets',{{'full'}},'selectable',true);
scenarios = repmat(base,1,15);

scenarios(2).name = 'NR_MODULE_30_R1';
scenarios(2).initial_offset_deg = 30;
scenarios(3).name = 'NR_MODULE_25_COMPACT_R1';
scenarios(3).spacing_lambda = .25;
scenarios(3).initial_offset_deg = 25;
scenarios(4).name = 'NR_MODULE_25_WIDE_R1';
scenarios(4).spacing_lambda = .5;
scenarios(4).initial_offset_deg = 25;
scenarios(5).name = 'NR_MODULE_25_REVERSED_R1';
scenarios(5).direction = -moduleByColumn;
scenarios(5).initial_offset_deg = 25;
scenarios(6).name = 'NR_MODULE_25_ROW_R1';
scenarios(6).direction = moduleByRow;
scenarios(6).initial_offset_deg = 25;
scenarios(7).name = 'NR_MODULE_25_ALL_R1';
scenarios(7).initial_offset_deg = 25;
scenarios(7).sets = {'non','partial','full'};
scenarios(8).name = 'NR_MODULE_25_ALL_R2';
scenarios(8).initial_offset_deg = 25;
scenarios(8).sets = {'non','partial','full'};
scenarios(8).rank = 2;
scenarios(9).name = 'NR_COMPACT_FAST_IID_R1';
scenarios(9).direction = [];
scenarios(9).initial_offset_deg = NaN;
scenarios(9).initial_jitter_deg = 0;

scenarios(10).name = 'NR_MODULE_0_CONTROL';
scenarios(10).initial_offset_deg = 0;
scenarios(10).selectable = false;
scenarios(11).name = 'NR_MODULE_OPPOSING_CONTROL';
scenarios(11).initial_offset_deg = -25;
scenarios(11).selectable = false;
scenarios(12).name = 'NR_HANDHELD_REFERENCE';
scenarios(12).spacing_lambda = NaN;
scenarios(12).tau_ms = 1.5;
scenarios(12).phase20_deg = 38.2;
scenarios(12).amp20_dB = 3.8;
scenarios(12).direction = [];
scenarios(12).initial_offset_deg = NaN;
scenarios(12).initial_jitter_deg = 0;
scenarios(12).selectable = false;

% The vertical module axis was identified using a separate, noiseless
% wideband channel diagnostic. The link screens still use fresh channel seeds.
scenarios(13).name = 'NR_MODULE_VERTICAL_0_R1';
scenarios(13).direction = moduleByRow;
scenarios(13).initial_offset_deg = 0;
scenarios(14).name = 'NR_MODULE_VERTICAL_10_COMPACT_R1';
scenarios(14).direction = moduleByRow;
scenarios(14).spacing_lambda = .25;
scenarios(14).initial_offset_deg = 10;
scenarios(15).name = 'NR_MODULE_VERTICAL_15_COMPACT_R1';
scenarios(15).direction = moduleByRow;
scenarios(15).spacing_lambda = .22;
scenarios(15).initial_offset_deg = 15;
end
