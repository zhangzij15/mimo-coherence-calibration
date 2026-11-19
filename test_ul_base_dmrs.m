function tests = test_ul_base_dmrs
%TEST_UL_BASE_DMRS Cover the native pilot generator and optional Toolbox parity.
tests = functiontests(localfunctions);
end

function testScramblingAndAbsoluteCrbReference(testCase)
cfg = pilotConfig(1,4,7,65535,1,[]);
[pilot,groups] = ul4tBaseDMRS(cfg,19,1);

% TS 38.211 sequence independently evaluated with 31-bit shift registers.
symbols = ([1 -1 1 1 -1 -1; -1 1 -1 -1 -1 -1].' + ...
    1i*[1 1 1 1 -1 -1; 1 1 -1 1 1 -1].')*10^(3/20)/sqrt(2);
expected = complex(zeros(12,2,4));
expected(1:2:end,:,1) = symbols;
expected(2:2:end,:,2) = symbols;
expected(1:2:end,:,3) = symbols.*[1 -1];
expected(2:2:end,:,4) = symbols.*[1 -1];
verifyEqual(testCase,groups,2);
verifyEqual(testCase,pilot,expected,'AbsTol',1e-12);
end

function testNative288RbIsContinuousAndPortOrthogonal(testCase)
cfg = pilotConfig(288,4,9,42,0,[144 144]);
pilot = ul4tBaseDMRS(cfg,7,0);
verifySize(testCase,pilot,[3456 2 4]);
for group = 1:2
    window = cfg;
    window.ofdm.num_rb = 144;
    window.ofdm.num_tones = 1728;
    window.ul.prb_start_crb = cfg.ul.prb_start_crb+144*(group-1);
    window.ul.dmrs_segments_rb = [];
    expected = ul4tBaseDMRS(window,7,0);
    verifyEqual(testCase,pilot((1:1728)+1728*(group-1),:,:),expected);
end
flattened = reshape(pilot,[],4);
verifyEqual(testCase,flattened'*flattened,3456*10^(3/10)*eye(4),'AbsTol',1e-8);

cfg.ul.rank = 2;
pilot = ul4tBaseDMRS(cfg,7,1);
ratio = pilot(1:2:end,:,2)./pilot(1:2:end,:,1);
verifyEqual(testCase,ratio,repmat(repmat([1;-1],864,1),1,2),'AbsTol',1e-12);
end

function testReferenceParityIfToolboxAvailable(testCase)
assumeTrue(testCase,~isempty(which('nrPUSCHDMRS')), ...
    'Toolbox reference is optional for this test.');
configs = {pilotConfig(5,1,0,1,0,[]), ...
    pilotConfig(240,1,13,417,1,[]), ...
    pilotConfig(288,1,9,65535,1,[100 88 100])};
slots = [0 8 19];
for k = 1:numel(configs)
    cfg = configs{k};
    for rank = 1:4
        cfg.ul.rank = rank;
        for startSymbol = [0 1]
            [actual,groups] = ul4tBaseDMRS(cfg,slots(k),startSymbol);
            [expected,expectedGroups] = ul4tToolboxDMRS(cfg,slots(k),startSymbol);
            verifyEqual(testCase,groups,expectedGroups);
            verifyEqual(testCase,actual,expected,'AbsTol',1e-12);
        end
    end
end
end

function cfg = pilotConfig(rb,rank,crb,nid,nscid,segments)
cfg.ofdm = struct('num_rb',rb,'num_tones',12*rb,'scs_Hz',30e3);
cfg.ul = struct('rank',rank,'prb_start_crb',crb,'nid',nid,'nscid',nscid, ...
    'dmrs_segments_rb',segments);
end
