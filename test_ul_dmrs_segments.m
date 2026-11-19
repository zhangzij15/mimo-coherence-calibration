function tests = test_ul_dmrs_segments
%TEST_UL_DMRS_SEGMENTS Check Toolbox DM-RS parity and 4T256R integration.
%   results = runtests('test_ul_dmrs_segments.m');
tests = functiontests(localfunctions);
end

function setup(testCase)
assumeNotEmpty(testCase,which('nrPUSCHDMRS'),'Optional 5G Toolbox reference.');
end

function testSegmentedGridMatchesSingleToolboxGrid(testCase)
% A zero-based sequence restart in the second window changes its symbols.
cfg = configUL4T('toolbox');
cfg.ofdm.num_rb = 240;
cfg.ul.prb_start_crb = 7;
cfg.ul.dmrs_segments_rb = [120 120];
cfg = prepareFixedArrayConfig(cfg);
for rank = [1 2 4]
    cfg.ul.rank = rank;
    for startSymbol = [0 1]
        [segmented,groups] = ul4tToolboxDMRS(cfg,11,startSymbol);
        expected = fullGridReference(cfg,11,startSymbol);
        verifyEqual(testCase,groups,2);
        verifyEqual(testCase,segmented,expected,'AbsTol',1e-12);
    end
end
end

function test288Rb256RxUsesEveryTone(testCase)
cfg = configUL4T('toolbox');
cfg.Nh = 16;
cfg.Nv = 16;
cfg.ofdm.num_rb = 288;
cfg.ul.dmrs_segments_rb = [144 144];
cfg.ul.duration_ms = 10.5; % One PUSCH keeps the full-dimension test short.
cfg.ul.ue_speed_kmh = 0;
cfg.ul.rx_noise_power = 10^(18/10);
cfg.ul.cases = cfg.ul.cases(3:4);
cfg.MC = 1;
cfg.save_results = false;
result = simulateUL4T(cfg);
verifyEqual(testCase,result.cfg.N,256);
verifyEqual(testCase,result.cfg.ofdm.num_tones,3456);
verifyEqual(testCase,result.cfg.ofdm.active_bandwidth_Hz,103680000);
verifyEqual(testCase,result.channel_diagnostics(1).selected_path_gains_size(3:4),[4 256]);
verifyEqual(testCase,result.resource_grid_definition.dmrs_segments_rb,[144 144]);
verifyTrue(testCase,result.resource_grid_definition.is_extended_research_grid);
verifyTrue(testCase,all(isfinite(result.fixed_beta_eff_sinr_dB(:))));
verifyTrue(testCase,all(isfinite(result.thp_frame_Mbps(:))));
end

function test288RbRequiresExplicitSegments(testCase)
cfg = configUL4T('toolbox');
cfg.ofdm.num_rb = 288;
verifyError(testCase,@() simulateUL4T(cfg),'UL4T:DMRSGrid');
cfg.ul.dmrs_segments_rb = [144 143];
verifyError(testCase,@() simulateUL4T(cfg),'UL4T:DMRSGrid');
cfg.ul.dmrs_segments_rb = [276 12];
verifyError(testCase,@() simulateUL4T(cfg),'UL4T:DMRSGrid');
end

function pilot = fullGridReference(cfg,slot,startSymbol)
carrier = nrCarrierConfig('SubcarrierSpacing',30,'NSizeGrid',240, ...
    'NStartGrid',7,'NSlot',slot);
pusch = nrPUSCHConfig('NSizeBWP',240,'NStartBWP',7,...
    'PRBSet',0:239,'NumLayers',cfg.ul.rank,'TransmissionScheme','nonCodebook',...
    'MappingType','B','SymbolAllocation',[startSymbol 14-startSymbol]);
ports = {0,[0 1],[0 1 4],[0 2 4 6]};
pusch.DMRS = nrPUSCHDMRSConfig('DMRSConfigurationType',1,'DMRSLength',2, ...
    'DMRSAdditionalPosition',0,'NumCDMGroupsWithoutData',2, ...
    'DMRSPortSet',ports{cfg.ul.rank},'NIDNSCID',1,'NSCID',0);
symbols = nrPUSCHDMRS(carrier,pusch);
indices = nrPUSCHDMRSIndices(carrier,pusch);
grid = complex(zeros(2880,14,cfg.ul.rank));
grid(indices(:)) = 10^(3/20)*symbols(:)./abs(symbols(:));
pilot = grid(:,startSymbol+(1:2),:);
end
