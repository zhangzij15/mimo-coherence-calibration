function [phaseDeg,ampDb,initialPhaseDeg] = ul4tBoundedRfDrift(nPorts,time_ms,selectedPorts,p)
%UL4TBOUNDEDRFDRIFT Generate continuous RF trajectories for four selected ports.
% Return absolute phase (degrees) and power gain (dB) on the supplied time grid.
% The caller controls the random stream and subtracts its actual SRS references.
% phase20_deg and amp20_dB cap the change in relative errors of any selected
% port pair over 20 ms. Bounded and load modes maintain those caps at every lag.
% The slew-mode phase can exceed its 20 ms cap at longer lags.
% Load mode requires a logical p.txActive value at each time-grid point.
% It uses piecewise-constant excitation and returns states at slot starts.
% An optional p.direction specifies port weights in selectedPorts order.
% The optional third output uses p.initial_offset_deg times those weights
% and independent uniform +/- p.initial_jitter_deg per selected port.
% Both initial phase options default to zero and leave drift unchanged.
% These models probe sensitivity under different drift mechanisms.

validateattributes(nPorts,{'numeric'},{'real','finite','scalar','integer','>=',4}, ...
    mfilename,'nPorts');
validateattributes(time_ms,{'numeric'},{'real','finite','vector','nonempty'}, ...
    mfilename,'time_ms');
time_ms = double(time_ms(:).');
assert(time_ms(1)==0 && all(diff(time_ms)>0),'ULBoundedRf:TimeGrid', ...
    'RF times must start at zero and increase strictly.');
validateattributes(selectedPorts,{'numeric'}, ...
    {'real','finite','vector','integer','numel',4,'>=',1,'<=',nPorts}, ...
    mfilename,'selectedPorts');
assert(numel(unique(selectedPorts))==4,'ULBoundedRf:Ports', ...
    'Select four distinct RF ports.');

fields = {'mode','grouping','phase20_deg','amp20_dB'};
assert(isstruct(p) && isscalar(p) && all(isfield(p,fields)), ...
    'ULBoundedRf:Parameters','Set mode, grouping, phase20_deg and amp20_dB.');
mode = string(p.mode);
grouping = string(p.grouping);
assert(isscalar(mode) && any(mode == ["slew","bounded","load"]), ...
    'ULBoundedRf:Mode','Set mode to slew, bounded or load.');
assert(isscalar(grouping) && any(grouping == ["balanced","spread"]), ...
    'ULBoundedRf:Grouping','Set grouping to balanced or spread.');
validateattributes(p.phase20_deg,{'numeric'}, ...
    {'real','finite','scalar','nonnegative','<=',40},mfilename,'p.phase20_deg');
validateattributes(p.amp20_dB,{'numeric'}, ...
    {'real','finite','scalar','nonnegative','<=',4},mfilename,'p.amp20_dB');

if isfield(p,'direction')
    validateattributes(p.direction,{'numeric'}, ...
        {'real','finite','vector','numel',4,'>=',-0.5,'<=',0.5}, ...
        mfilename,'p.direction');
    direction = double(p.direction(:));
else
    switch grouping
        case "balanced"
            direction = [-0.5;-0.5;0.5;0.5];
            direction = direction(randperm(4));
        case "spread"
            raw = randn(4,1);
            limits = [min(raw),max(raw)];
            if limits(1)==limits(2)
                direction = zeros(4,1);
            else
                direction = (raw-mean(limits))/diff(limits);
            end
    end
end
% Vary stress level across realizations independently of the channel.
phaseScale = double(p.phase20_deg)*(0.83+0.17*rand)*direction;
ampScale = double(p.amp20_dB)*(0.83+0.17*rand)*direction;
if mode=="load"
    state = loadState(time_ms,p);
    phaseShape = state;
    ampShape = state;
else
    % The sinusoidal half-range makes amplitude bounded at every lag.
    ampPeriod_ms = 60+40*rand;
    ampShape = 0.5*sin(2*pi*time_ms/ampPeriod_ms+2*pi*rand);
    if mode=="slew"
        phaseShape = time_ms/20;
    else
        phasePeriod_ms = 22+8*rand;
        phaseShape = 0.5*sin(2*pi*time_ms/phasePeriod_ms+2*pi*rand);
    end
end

phaseDeg = zeros(nPorts,numel(time_ms));
ampDb = zeros(nPorts,numel(time_ms));
phaseDeg(selectedPorts,:) = phaseScale*phaseShape;
ampDb(selectedPorts,:) = ampScale*ampShape;
if nargout>=3
    initialPhaseDeg = zeros(nPorts,1);
    if isfield(p,'initial_offset_deg')
        validateattributes(p.initial_offset_deg,{'numeric'}, ...
            {'real','finite','scalar'},mfilename,'p.initial_offset_deg');
        initialPhaseDeg(selectedPorts) = double(p.initial_offset_deg)*direction;
    end
    if isfield(p,'initial_jitter_deg')
        validateattributes(p.initial_jitter_deg,{'numeric'}, ...
            {'real','finite','scalar','nonnegative'},mfilename,'p.initial_jitter_deg');
        if p.initial_jitter_deg>0
            initialPhaseDeg(selectedPorts) = initialPhaseDeg(selectedPorts)+ ...
                double(p.initial_jitter_deg)*(2*rand(4,1)-1);
        end
    end
end
end

function state = loadState(time_ms,p)
assert(isfield(p,'txActive'),'ULBoundedRf:LoadSchedule', ...
    'Load mode requires a txActive vector on the RF time grid.');
validateattributes(p.txActive,{'logical'}, ...
    {'vector','numel',numel(time_ms)},mfilename,'p.txActive');
tau_ms = 1.5;
if isfield(p,'tau_ms'), tau_ms = p.tau_ms; end
validateattributes(tau_ms,{'numeric'},{'real','finite','scalar','positive'}, ...
    mfilename,'p.tau_ms');
relax = -expm1(-diff(time_ms)/double(tau_ms));
state = zeros(size(time_ms));
for k = 2:numel(time_ms)
    % The past slot drives the state. Events such as SRS do not reset it.
    state(k) = state(k-1)+(double(p.txActive(k-1))-state(k-1))*relax(k-1);
end
end
