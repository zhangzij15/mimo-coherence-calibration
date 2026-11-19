function tests = test_ul_bounded_rf_drift
%TEST_UL_BOUNDED_RF_DRIFT Verify per-pair drift envelopes and trajectory semantics.
tests = functiontests(localfunctions);
end

function testSlewRespectsEveryTwentyMillisecondWindow(testCase)
time = [0:0.5:135,135.1,135.17];
for group = ["balanced","spread"]
    p = settings('slew',group);
    [phase,amp] = draw(5,time,1:4,p,28849);
    verifyLessThanOrEqual(testCase,maxPairChange(phase(1:4,:),time,20),30+1e-10);
    verifyLessThanOrEqual(testCase,maxPairChange(amp(1:4,:),time,Inf),3+1e-10);
    % The model retains slope after 20 ms rather than clamping or resetting.
    verifyGreaterThan(testCase,maxPairChange(phase(1:4,:),time,Inf),60);
end
end

function testBoundedModeRespectsEveryLag(testCase)
time = 0:0.5:235;
for group = ["balanced","spread"]
    p = settings('bounded',group);
    [phase,amp] = draw(4,time,1:4,p,717);
    phaseChange = maxPairChange(phase,time,Inf);
    ampChange = maxPairChange(amp,time,Inf);
    verifyGreaterThan(testCase,phaseChange,0.5);
    verifyGreaterThan(testCase,ampChange,0.05);
    verifyLessThanOrEqual(testCase,phaseChange,30+1e-10);
    verifyLessThanOrEqual(testCase,ampChange,3+1e-10);
end
end

function testSeedAndPortSelection(testCase)
ports = [7 2 6 4];
time = 0:0.5:94;
p = settings('slew','balanced');
[phase,amp] = draw(8,time,ports,p,2701);
[repeatedPhase,repeatedAmp] = draw(8,time,ports,p,2701);
[differentPhase,differentAmp] = draw(8,time,ports,p,2702);
verifySize(testCase,phase,[8 numel(time)]);
verifySize(testCase,amp,[8 numel(time)]);
verifyEqual(testCase,phase,repeatedPhase);
verifyEqual(testCase,amp,repeatedAmp);
verifyNotEqual(testCase,phase(ports,:),differentPhase(ports,:));
verifyNotEqual(testCase,amp(ports,:),differentAmp(ports,:));
verifyEqual(testCase,phase([1 3 5 8],:),zeros(4,numel(time)));
verifyEqual(testCase,amp([1 3 5 8],:),zeros(4,numel(time)));
end

function testContinuousTimeDoesNotDependOnSamplingDates(testCase)
time = [0 4.999 5 5.001 20 39.999 40 40.001 54.999 55 55.001 85];
coarse = [1 3 6 9 12];
closePairs = [2 3;6 7;9 10];
for mode = ["slew","bounded"]
    p = settings(mode,'balanced');
    [phase,amp] = draw(4,time,1:4,p,18887);
    [coarsePhase,coarseAmp] = draw(4,time(coarse),1:4,p,18887);
    verifyEqual(testCase,phase(:,coarse),coarsePhase,'AbsTol',1e-12);
    verifyEqual(testCase,amp(:,coarse),coarseAmp,'AbsTol',1e-12);
    verifyLessThan(testCase,max(abs(phase(:,closePairs(:,1))-phase(:,closePairs(:,2))),[],'all'),0.02);
    verifyLessThan(testCase,max(abs(amp(:,closePairs(:,1))-amp(:,closePairs(:,2))),[],'all'),0.005);
end
end

function testLoadModeIsCausalAndRetainsThermalState(testCase)
time = [0 1 2 3 4 5 6 100];
p = settings('load','balanced');
p.tau_ms = 1.5;
p.txActive = logical([0 1 1 0 0 1 0 0]);
[phase,amp] = draw(4,time,1:4,p,2026);
verifyEqual(testCase,phase(:,1:2),zeros(4,2));
verifyEqual(testCase,amp(:,1:2),zeros(4,2));
verifyGreaterThan(testCase,norm(phase(:,3)),0);
verifyGreaterThan(testCase,norm(amp(:,3)),0);
% The analytic two-slot heating state is (1 + exp(-1/tau)) times
% the one-slot state. The following idle slot decays with exp(-1/tau).
decay = exp(-1/1.5);
verifyEqual(testCase,phase(:,4),(1+decay)*phase(:,3),'AbsTol',1e-12);
verifyEqual(testCase,amp(:,5),decay*amp(:,4),'AbsTol',1e-12);
verifyEqual(testCase,phase(:,5),decay*phase(:,4),'AbsTol',1e-12);
verifyGreaterThan(testCase,norm(phase(:,7)),norm(phase(:,6)));
verifyLessThan(testCase,norm(phase(:,8)),1e-12);
verifyLessThanOrEqual(testCase,maxPairChange(phase,time,Inf),30+1e-10);
verifyLessThanOrEqual(testCase,maxPairChange(amp,time,Inf),3+1e-10);
end

function testNearZeroAndZeroLimitsStaySmall(testCase)
time = 0:0.5:60;
for mode = ["slew","bounded","load"]
    p = settings(mode,'spread');
    if mode=="load", p.txActive = time>=6 & time<38; end
    p.phase20_deg = 0.5;
    p.amp20_dB = 0.1;
    [phase,amp] = draw(4,time,1:4,p,29);
    verifyGreaterThan(testCase,maxPairChange(phase,time,20),0);
    verifyGreaterThan(testCase,maxPairChange(amp,time,20),0);
    verifyLessThanOrEqual(testCase,maxPairChange(phase,time,20),0.5+1e-10);
    verifyLessThanOrEqual(testCase,maxPairChange(amp,time,20),0.1+1e-10);

    p.phase20_deg = 0;
    p.amp20_dB = 0;
    [phase,amp] = draw(4,time,1:4,p,29);
    verifyEqual(testCase,phase,zeros(4,numel(time)));
    verifyEqual(testCase,amp,zeros(4,numel(time)));
end
end

function p = settings(mode,group)
p = struct('mode',mode,'phase20_deg',30,'amp20_dB',3,'grouping',group);
end

function [phase,amp] = draw(nPorts,time,selectedPorts,p,seed)
oldState = rng;
restoreRng = onCleanup(@() rng(oldState));
rng(seed,'twister');
[phase,amp] = ul4tBoundedRfDrift(nPorts,time,selectedPorts,p);
end

function maxChange = maxPairChange(data,time,maxLag)
withinLag = abs(time(:)-time(:).')<=maxLag+1e-10;
maxChange = 0;
pairs = nchoosek(1:size(data,1),2);
for k = 1:size(pairs,1)
    delta = data(pairs(k,1),:) - data(pairs(k,2),:);
    difference = abs(delta(:)-delta(:).');
    maxChange = max(maxChange,max(difference(withinLag)));
end
end
