function tests = test_ul_nr_rf_initial_phase
%TEST_UL_NR_RF_INITIAL_PHASE Check fixed RF module directions and cold phases.
tests = functiontests(localfunctions);
end

function testFixedAxisMapsToSelectedPortsAndPreservesGlobalEnvelope(testCase)
time = 0:0.5:45;
ports = [7 2 6 4];
p = loadSettings(time);
p.initial_offset_deg = 24;
[phase,amp,cold] = sampleRf(8,time,ports,p);
testCase.verifyEqual(cold(ports),[-12;-12;12;12]);
testCase.verifyEqual(cold([1 3 5 8]),zeros(4,1));
testCase.verifyEqual(phase(ports,1),zeros(4,1));
testCase.verifyEqual(phase(ports(1),:),phase(ports(2),:));
testCase.verifyEqual(amp(ports(3),:),amp(ports(4),:));
testCase.verifyGreaterThan(phase(ports(3),40)-phase(ports(1),40),30);
testCase.verifyLessThan(phase(ports(1),40),0);
testCase.verifyEqual(phase([1 3 5 8],:),zeros(4,numel(time)));
testCase.verifyEqual(amp([1 3 5 8],:),zeros(4,numel(time)));
for pair = nchoosek(ports,2).'
    phasePair = phase(pair(1),:)-phase(pair(2),:);
    ampPair = amp(pair(1),:)-amp(pair(2),:);
    testCase.verifyLessThanOrEqual(max(phasePair)-min(phasePair),p.phase20_deg);
    testCase.verifyLessThanOrEqual(max(ampPair)-min(ampPair),p.amp20_dB);
end
end

function testColdOffsetAndJitterLeaveRfTrajectoryPaired(testCase)
time = 0:0.5:10;
p = loadSettings(time);
[referencePhase,referenceAmp,defaultCold] = sampleRf(4,time,1:4,p);
p.initial_offset_deg = -24;
[phase,amp,cold] = sampleRf(4,time,1:4,p);
testCase.verifyEqual(defaultCold,zeros(4,1));
testCase.verifyEqual(phase,referencePhase);
testCase.verifyEqual(amp,referenceAmp);
testCase.verifyEqual(cold,[12;12;-12;-12]);
p.initial_jitter_deg = 2;
[phase,amp,jittered] = sampleRf(4,time,1:4,p);
[~,~,repeated] = sampleRf(4,time,1:4,p);
testCase.verifyEqual(jittered,repeated);
testCase.verifyGreaterThan(max(abs(jittered-cold)),0);
testCase.verifyLessThanOrEqual(abs(jittered-cold),2*ones(4,1));
testCase.verifyEqual(phase,referencePhase);
testCase.verifyEqual(amp,referenceAmp);
end

function testColdPhaseReachesEveryNrCaseAndKeepsDefaultSeed(testCase)
previous = rng;
cleanup = onCleanup(@() rng(previous));
cfg = configUL4T('base');
cfg.Nh = 8; cfg.Nv = 8; cfg.ofdm.num_rb = 240;
cfg.MC = 1;
cfg.seed = 20261117;
cfg.ul.duration_ms = 10.5;
cfg.ul.rx_noise_power = 100;
cfg.ul.cases = cfg.ul.cases(3:4);
cfg.save_results = false;
cfg.ul.coherence_drift.bounded = rmfield(loadSettings(0:0.5:10),'txActive');
reference = simulateUL4T(cfg);

rng(mod(double(cfg.seed)+2^28,2^32),'twister');
legacyCold = 2*pi*rand(8,1);
legacyCold = legacyCold([1 7 3 5]);
testCase.verifyEqual(reference.srs_initial_phase_rad,repmat(legacyCold,1,1,1,2));

cfg.ul.coherence_drift.bounded.initial_offset_deg = -24;
crossed = simulateUL4T(cfg);
coldPhase = (pi/180)*[12;12;-12;-12];
testCase.verifyEqual(crossed.srs_initial_phase_rad,repmat(coldPhase,1,1,1,2), ...
    'AbsTol',1e-12);
testCase.verifyEqual(crossed.coherence_gain,reference.coherence_gain);
testCase.verifyEqual(crossed.coherence_drift,reference.coherence_drift);
end

function p = loadSettings(time)
p = struct('mode','load','grouping','balanced','tau_ms',0.5, ...
    'phase20_deg',38.2,'amp20_dB',3.8,'direction',[-0.5;-0.5;0.5;0.5], ...
    'txActive',time>=1 & time<20);
end

function [phase,amp,cold] = sampleRf(nPorts,time,ports,p)
previous = rng;
cleanup = onCleanup(@() rng(previous));
rng(1244,'twister');
[phase,amp,cold] = ul4tBoundedRfDrift(nPorts,time,ports,p);
end
