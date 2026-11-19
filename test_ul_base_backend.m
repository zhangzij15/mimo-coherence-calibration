function tests = test_ul_base_backend
tests = functiontests(localfunctions);
end

function testStandaloneDefaultUsesConfirmedLink(testCase)
plan = runNrStrict256(struct('dryRun',true, ...
    'scenarioNames',{{'NR_MODULE_VERTICAL_15_COMPACT_R1'}}));
cfg = plan.plan.screen_configs{1};
testCase.verifyEqual(cfg.cdl.backend,'base-cdl');
testCase.verifyEqual(cfg.Nh*cfg.Nv,256);
testCase.verifyEqual(cfg.ofdm.num_rb,288);
testCase.verifyEmpty(cfg.ul.dmrs_segments_rb);
testCase.verifyEqual(cfg.fc_Hz,6.7e9);
testCase.verifyEqual({cfg.ul.cases.calibration}, ...
    {'none','none','feedback','reselect','oracle'});
end

function testStandaloneLinkPowerAndZeroDrift(testCase)
plan = runNrStrict256(struct('dryRun',true,'screenMC',1, ...
    'scenarioNames',{{'NR_MODULE_VERTICAL_15_COMPACT_R1'}}));
cfg = plan.plan.screen_configs{1};
testCase.assertEqual(cfg.cdl.backend,'base-cdl');
cfg.ul.coherence_drift.bounded.phase20_deg = 0;
cfg.ul.coherence_drift.bounded.amp20_dB = 0;
r = simulateUL4T(cfg);
testCase.verifyEqual(r.channel_diagnostics(1).backend,'base-cdl');
testCase.verifyEqual(r.resource_grid_definition.num_tones,3456);
testCase.verifyEqual(r.schedule.srs_ms,5);
testCase.verifyEqual(nnz(r.schedule.srs_age_ms<=20),10);
testCase.verifySize(r.good_bits,[1 20 1 5]);
testCase.verifyEqual(r.good_bits(:,:,:,1),r.good_bits(:,:,:,2));
testCase.verifyEqual(r.fixed_beta_eff_sinr_dB(:,:,:,1), ...
    r.fixed_beta_eff_sinr_dB(:,:,:,2),'AbsTol',1e-10);
testCase.verifyLessThan(max(abs(r.pusch_tx_power(:)-1)),1e-10);
end

function testScenarioSearchAcceptsAllCdlProfiles(testCase)
profiles = {'CDL-A','CDL-B','CDL-C','CDL-D','CDL-E'};
for k = 1:numel(profiles)
    report = runNrStrict256(struct('dryRun',true, ...
        'channelProfile',profiles{k}, ...
        'scenarioNames',{{'NR_MODULE_VERTICAL_15_COMPACT_R1'}}));
    cfg = report.plan.screen_configs{1};
    testCase.verifyEqual(cfg.cdl.toolbox_profile,profiles{k});
    testCase.verifyEqual(report.plan.channel_profile,profiles{k});
end
end

function testLosPathUsesTheExistingPuschReceiver(testCase)
report = runNrStrict256(struct('dryRun',true,'screenMC',1, ...
    'channelProfile','CDL-E', ...
    'scenarioNames',{{'NR_MODULE_VERTICAL_15_COMPACT_R1'}}));
cfg = report.plan.screen_configs{1};
cfg.ul.duration_ms = 10.5;
cfg.Nh = 4; cfg.Nv = 4;
r = simulateUL4T(cfg);
testCase.verifyEqual(r.channel_diagnostics(1).profile,'CDL-E');
testCase.verifyEqual(r.channel_diagnostics(1).first_cluster_k_factor_dB,22,'AbsTol',1e-12);
testCase.verifyEqual(r.resource_grid_definition.num_tones,3456);
testCase.verifySize(r.good_bits,[1 1 1 5]);
testCase.verifyTrue(all(isfinite(r.fixed_beta_eff_sinr_dB),'all'));
end
