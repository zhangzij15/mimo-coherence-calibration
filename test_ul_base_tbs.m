function tests = test_ul_base_tbs
tests = functiontests(localfunctions);
end

function testQuantizationAndEmptyAllocation(t)
% The five active cases cover both quantizers and the segmentation branches.
verifyEqual(t,ul4tBaseTBS(2,1,1,69,.5),64);
verifyEqual(t,ul4tBaseTBS(2,1,16,120,3824/3840),3824);
verifyEqual(t,ul4tBaseTBS(2,1,16,120,1),3840);
verifyEqual(t,ul4tBaseTBS(2,1,64,120,.25),3848);
verifyEqual(t,ul4tBaseTBS(4,4,52,120,.48),48168);
verifyEqual(t,ul4tBaseTBS([2 4],1,0,144,[.12 .32]),[0 0]);
verifyEqual(t,ul4tBaseTBS([2 4],1,288,0,[.12 .32]),[0 0]);
end

function testProjectLinkMappingWithoutToolbox(t)
cfg = struct('cdl',struct('backend','base-matlab'), ...
    'ofdm',struct('num_rb',288),'ul',struct('rank',1,'mapping','eesm'));
main = ul4tLinkMapping('prepare',cfg,2,0);
prefix = ul4tLinkMapping('prepare',cfg,2,1);
verifyEqual(t,main.data_re_per_rb,144);
verifyEqual(t,prefix.data_re_per_rb,132);
verifyNumElements(t,main.bits,29);
verifyEqual(t,main.bits(1),9744);
verifyEqual(t,prefix.bits(1),8976);
verifyTrue(t,all(main.bits>=prefix.bits));
verifyEqual(t,main.bg(1),2);
verifyEqual(t,main.num_cb(1),3);
end

function testAgainstToolboxWhenInstalled(t)
assumeNotEmpty(t,which('nrTBS'));
[modulation,qm,rates] = mcsTable1;
for layers = 1:4
    for rb = [1 12 52 240 275 288]
        for re = [12 132 144 164]
            expected = arrayfun(@(i)nrTBS(modulation{qm(i)/2},layers,rb,re,rates(i)),1:numel(rates));
            actual = ul4tBaseTBS(qm,layers,rb,re,rates);
            verifyEqual(t,actual,expected,sprintf('layers=%d RB=%d RE=%d',layers,rb,re));
        end
    end
end
cfg = struct('cdl',struct('backend','5g-toolbox'), ...
    'ofdm',struct('num_rb',288),'ul',struct('rank',1,'mapping','eesm'));
for layers = 1:4
    cfg.ul.rank = layers;
    for prefix = 0:1
        actual = ul4tLinkMapping('prepare',cfg,2,prefix);
        expected = arrayfun(@(i)nrTBS(modulation{qm(i)/2},layers,288, ...
            actual.data_re_per_rb,rates(i)),1:numel(rates));
        verifyEqual(t,actual.bits,expected);
    end
end
end

function [modulation,qm,rates] = mcsTable1
modulation = {'QPSK','16QAM','64QAM'};
qm = [2*ones(1,10),4*ones(1,7),6*ones(1,12)];
rates = [120 157 193 251 308 379 449 526 602 679 340 378 434 490 553 ...
    616 658 438 466 517 567 616 666 719 772 822 873 910 948]/1024;
end
