function tests = test_ul_dmrs_cache
%TEST_UL_DMRS_CACHE Real cached/uncached regressions on the full 240-RB grid.
% Run: results = runtests('test_ul_dmrs_cache.m');
%      assert(all([results.Passed]),'DMRS cache regression failed.');
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
cfg = configUL4T('base');
cfg.Nh = 8; cfg.Nv = 8; cfg.ofdm.num_rb = 240;
cfg.MC = 2;
cfg.seed = 20261117;
cfg.progress_every = 2;
cfg.ul.duration_ms = 10.5; % One PUSCH per MC, without reducing bandwidth.
cfg.ul.rx_noise_power = 10^(14/10);
cfg.ul.show_figure = false;
cfg.ul.save_figure = false;
cfg.save_results = false;
cfg.users.mode = 'random';
cfg.ul.cases = fullCases();
for rankValue = [1 4]
    cfg.ul.rank = rankValue;
    cfg.ul.dmrs_filter_cache_enabled = true;
    testCase.TestData.(sprintf('cachedRank%d',rankValue)) = simulateUL4T(cfg);
    cfg.ul.dmrs_filter_cache_enabled = false;
    testCase.TestData.(sprintf('uncachedRank%d',rankValue)) = simulateUL4T(cfg);
end

% Non/partial rank-4 codebooks have unequal row support per layer, so a
% reselected precoder under unequal estimated RF amplitudes needs another K.
cfg.MC = 1;
cfg.ul.rank = 4;
cfg.ul.dmrs_filter_cache_enabled = true;
cfg.ul.cases = repmat(struct('precoder','nr','coherence_error',true, ...
    'sets',{{'non','partial'}},'calibration','none'),1,2);
cfg.ul.cases(2).calibration = 'reselect';
testCase.TestData.unequalPower = simulateUL4T(cfg);
end

function testRankOneCacheMatchesUncachedReceiver(testCase)
% Catches reuse of the wrong Wiener filter or an omitted/unused cache.
verifyEquivalent(testCase,testCase.TestData.cachedRank1, ...
    testCase.TestData.uncachedRank1);
end

function testRankFourCacheMatchesUncachedReceiver(testCase)
% Covers the multi-layer filter and its per-column power key.
verifyEquivalent(testCase,testCase.TestData.cachedRank4, ...
    testCase.TestData.uncachedRank4);
end

function testSecondMCUsesItsOwnPropagationPrior(testCase)
% A cache surviving the MC boundary would reuse MC1's different pRx.
cached = testCase.TestData.cachedRank4;
uncached = testCase.TestData.uncachedRank4;
p1 = cached.channel_diagnostics(1).cluster_pattern_power;
p2 = cached.channel_diagnostics(2).cluster_pattern_power;
testCase.verifyGreaterThan(max(abs(p1-p2)),1e-8, ...
    'The fixture must exercise distinct MC propagation priors.');
testCase.verifyEqual(cached.fixed_beta_eff_sinr_dB(2,:,:,:), ...
    uncached.fixed_beta_eff_sinr_dB(2,:,:,:),'AbsTol',1e-10);
testCase.assertTrue(isfield(cached,'dmrs_filter_cache'), ...
    'Per-MC cache statistics are required to verify invalidation.');
testCase.verifyEqual(cached.dmrs_filter_cache.misses,ones(2,1));
end

function testUnequalLayerPowerInvalidatesTheCache(testCase)
% Catches replacing the actual per-layer power key with a hardcoded 1/r.
result = testCase.TestData.unequalPower;
w = result.calibration.precoder(:,:,1,1,1,1);
gain = result.calibration.estimated_gain(:,1,1,1,1);
initial = exp(1i*result.srs_initial_phase_rad(:,1,1,2));
predicted = bsxfun(@times,gain.*initial,w);
predicted = predicted*(norm(w,'fro')/norm(predicted,'fro'));
power = sum(abs(predicted).^2,1);
testCase.verifyGreaterThan(max(power)-min(power),1e-8, ...
    'The fixture must produce genuinely unequal predicted layer powers.');
testCase.assertTrue(isfield(result,'dmrs_filter_cache'));
testCase.verifyEqual(result.dmrs_filter_cache.misses,2);
testCase.verifyEqual(result.dmrs_filter_cache.hits,0);
end

function verifyEquivalent(testCase,cached,uncached)
testCase.verifyEqual(cached.cfg.ofdm.num_rb,240);
testCase.verifyEqual(numel(cached.schedule.ul_ms),1);
testCase.verifyEqual(cached.ue_position_m,uncached.ue_position_m,'AbsTol',1e-12);
testCase.verifyEqual(cached.fixed_beta_eff_sinr_dB, ...
    uncached.fixed_beta_eff_sinr_dB,'AbsTol',1e-10);
testCase.verifyEqual(cached.dmrs_nmse,uncached.dmrs_nmse,'AbsTol',1e-11);
testCase.verifyEqual(cached.good_bits,uncached.good_bits);
testCase.assertTrue(isfield(cached,'dmrs_filter_cache') && ...
    isfield(uncached,'dmrs_filter_cache'), ...
    'The simulator must report actual cache hits and filter solves per MC.');
testCase.verifyEqual(cached.dmrs_filter_cache.misses,ones(2,1));
testCase.verifyEqual(cached.dmrs_filter_cache.hits,9*ones(2,1));
testCase.verifyEqual(uncached.dmrs_filter_cache.hits,zeros(2,1));
testCase.verifyEqual(uncached.dmrs_filter_cache.misses,10*ones(2,1));
end

function cases = fullCases()
cases = repmat(struct('precoder','svd','coherence_error',true, ...
    'sets',{{'full'}},'calibration','none'),1,10);
families = {'svd','nr'};
modes = {'none','none','feedback','reselect','oracle'};
for f = 1:2
    for k = 1:5
        c = (f-1)*5+k;
        cases(c).precoder = families{f};
        cases(c).coherence_error = k~=1;
        cases(c).calibration = modes{k};
    end
end
end
