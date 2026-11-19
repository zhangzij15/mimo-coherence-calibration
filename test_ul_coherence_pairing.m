function tests = test_ul_coherence_pairing
%TEST_UL_COHERENCE_PAIRING Base MATLAB regression with one PUSCH per MC.
% Run: results = runtests('test_ul_coherence_pairing.m');
%      assert(all([results.Passed]),'UL coherence pairing regression failed.');
% Retain the full 240-RB grid; shorten duration rather than bandwidth.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
cfg = configUL4T('base');
cfg.Nh = 8; cfg.Nv = 8; cfg.ofdm.num_rb = 240;
cfg.MC = 2;
cfg.seed = 20261117;
cfg.progress_every = 2;
cfg.ul.duration_ms = 10.5; % SRS at 5 ms and one PUSCH at 10 ms.
cfg.ul.rx_noise_power = 10^(14/10);
cfg.ul.show_figure = false;
cfg.ul.save_figure = false;
cfg.save_results = false;
cfg.users.mode = 'random';
zeroFields = {'sigmaPhi0Deg','sigmaA0dB','sigmaDeltaFHz', ...
    'qPhiDeg2PerSec','sigmaADynDb'};
for k = 1:numel(zeroFields)
    cfg.ul.coherence_drift.(zeroFields{k}) = 0;
end
testCase.TestData.rank1 = simulateUL4T(cfg);

rank4Cfg = cfg;
rank4Cfg.ul.rank = 4; % Changes DMRS/noise draw count in the first MC.
testCase.TestData.rank4 = simulateUL4T(rank4Cfg);

upaCfg = cfg;
upaCfg.ul.ue_layout = 'upa'; % Changes the number of candidate RF locations.
upaCfg.ul.ue_antenna_indices = [];
testCase.TestData.upa = simulateUL4T(upaCfg);

secondCfg = cfg;
secondCfg.MC = 1;
secondCfg.seed = cfg.seed+1;
testCase.TestData.secondMCAlone = simulateUL4T(secondCfg);

calCfg = cfg;
calCfg.MC = 1;
calCfg.ul.cases = calibrationCases();
testCase.TestData.calibration = simulateUL4T(calCfg);
end

function testDefaultZeroRFHardwareAndLinkArePaired(testCase)
% Catches independent per-case initial RF phases without calibration cases.
result = testCase.TestData.rank1;
testCase.verifyEqual(numel(result.schedule.ul_ms),1);
testCase.verifyEqual(result.cfg.ofdm.num_rb,240);
testCase.verifyEqual(result.coherence_gain,ones(size(result.coherence_gain)), ...
    'AbsTol',1e-12);
for pair = [1 3]
    testCase.verifyEqual(result.srs_initial_phase_rad(:,:,:,pair), ...
        result.srs_initial_phase_rad(:,:,:,pair+1),'AbsTol',1e-12);
    testCase.verifyEqual(result.pusch_effective_energy(:,:,pair), ...
        result.pusch_effective_energy(:,:,pair+1),'RelTol',1e-12);
    testCase.verifyEqual(result.eff_sinr_dB(:,:,:,pair), ...
        result.eff_sinr_dB(:,:,:,pair+1),'AbsTol',1e-10);
    testCase.verifyEqual(result.good_bits(:,:,:,pair), ...
        result.good_bits(:,:,:,pair+1));
end
end

function testLaterMCGeometryIsIndependentOfRankAndLayout(testCase)
% Catches previous-MC rank/layout-dependent RNG consumption changing geometry.
reference = testCase.TestData.rank1;
testCase.verifyEqual(testCase.TestData.rank4.ue_position_m, ...
    reference.ue_position_m,'AbsTol',1e-12);
testCase.verifyEqual(testCase.TestData.upa.ue_position_m, ...
    reference.ue_position_m,'AbsTol',1e-12);
% Checkpointing one MC at seed+index-1 must reproduce the corresponding MC.
second = testCase.TestData.secondMCAlone;
testCase.verifyEqual(second.ue_position_m,reference.ue_position_m(2,:), ...
    'AbsTol',1e-12);
testCase.verifyEqual(second.eff_sinr_dB,reference.eff_sinr_dB(2,:,:,:), ...
    'AbsTol',1e-10);
end

function testFixedBetaDiagnosticIsFiniteAndPaired(testCase)
% Catches an omitted/nonfinite diagnostic or using a different RF/noise case.
result = testCase.TestData.rank1;
testCase.assertTrue(isfield(result,'fixed_beta_eff_sinr_dB'), ...
    'A separate fixed-beta effective SINR diagnostic is required.');
testCase.assertTrue(isfield(result,'fixed_beta_linear'));
testCase.verifyEqual(result.fixed_beta_linear,1);
fixed = result.fixed_beta_eff_sinr_dB;
testCase.verifySize(fixed,size(result.eff_sinr_dB));
testCase.verifyTrue(all(isfinite(fixed(:))));
for pair = [1 3]
    testCase.verifyEqual(fixed(:,:,:,pair),fixed(:,:,:,pair+1), ...
        'AbsTol',1e-10);
end
end

function testCalibrationDiagnosticsUseKnownZeroRFReference(testCase)
% Catches gauge-dependent Oracle errors and reporting pilot energy as EPRE.
result = testCase.TestData.calibration;
required = {'estimation_phase_error_deg','estimation_amplitude_error_dB', ...
    'pilot_epre_min','pilot_epre_max','pilot_total_energy'};
for k = 1:numel(required)
    testCase.assertTrue(isfield(result.calibration,required{k}), ...
        ['Missing calibration diagnostic: ' required{k}]);
end
phase = result.calibration.estimation_phase_error_deg;
amplitude = result.calibration.estimation_amplitude_error_dB;
testCase.verifyTrue(all(isfinite(phase(:))) && all(isfinite(amplitude(:))));
testCase.verifyEqual(phase(:,:,:,:,3),zeros(size(phase(:,:,:,:,3))), ...
    'AbsTol',1e-10); % Oracle is calibrated case 3.
testCase.verifyEqual(amplitude(:,:,:,:,3),zeros(size(amplitude(:,:,:,:,3))), ...
    'AbsTol',1e-10);
testCase.verifyEqual(result.calibration.pilot_epre_min(:,:,:,1:2), ...
    ones(size(result.calibration.pilot_epre_min(:,:,:,1:2))),'AbsTol',1e-12);
testCase.verifyEqual(result.calibration.pilot_epre_max(:,:,:,1:2), ...
    ones(size(result.calibration.pilot_epre_max(:,:,:,1:2))),'AbsTol',1e-12);
testCase.verifyEqual(result.calibration.pilot_total_energy(:,:,:,1:2), ...
    128*ones(size(result.calibration.pilot_total_energy(:,:,:,1:2))), ...
    'AbsTol',1e-10);
testCase.verifyTrue(all(isnan(result.calibration.pilot_total_energy(:,:,:,3)), 'all'));
end

function cases = calibrationCases()
cases = repmat(struct('precoder','nr','coherence_error',true, ...
    'sets',{{'full'}},'calibration','none'),1,5);
modes = {'none','none','feedback','reselect','oracle'};
for k = 1:5
    cases(k).coherence_error = k~=1;
    cases(k).calibration = modes{k};
end
end
