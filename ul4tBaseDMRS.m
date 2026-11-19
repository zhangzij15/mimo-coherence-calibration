function [pilot,groupCount] = ul4tBaseDMRS(cfg,slotInFrame,startSymbol)
%UL4TBASEDMRS CP-OFDM PUSCH Type 1 DM-RS on a contiguous allocation.
%   Uses normal CP, mapping B, a front-loaded DM-RS pair, and two reserved
%   CDM groups. PILOT is indexed by local tone, DM-RS symbol, and layer.
%   TS 38.211 §§5.2.1, 6.4.1.1.1.1 and 6.4.1.1.3 define the sequence
%   and resource mapping. Baseline scrambling has dmrs-Uplink disabled.
if nargin < 3, startSymbol = 0; end
validateattributes(slotInFrame,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(startSymbol,{'numeric'},{'scalar','integer','>=',0,'<=',12});

numRB = double(cfg.ofdm.num_rb);
rankValue = double(cfg.ul.rank);
nid = double(cfg.ul.nid);
nscid = double(cfg.ul.nscid);
firstCRB = double(cfg.ul.prb_start_crb);
validateattributes(numRB,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(rankValue,{'numeric'},{'scalar','integer','>=',1,'<=',4});
validateattributes(nid,{'numeric'},{'scalar','integer','>=',0,'<=',65535});
validateattributes(nscid,{'numeric'},{'scalar','integer','>=',0,'<=',1});
validateattributes(firstCRB,{'numeric'},{'scalar','integer','nonnegative','finite'});

groupCount = 2;
pilotAmplitude = 10^(3*(groupCount-1)/20);
numTones = 12*numRB;
pilot = complex(zeros(numTones,2,rankValue));
portSets = {0,[0 1],[0 1 4],[0 2 4 6]};
ports = portSets{rankValue};

% Without transform precoding, Point A is the frequency reference.
% Each RB contributes six sequence values per CDM group.
firstSequenceSymbol = 6*firstCRB;
lastSequenceSymbol = firstSequenceSymbol+6*numRB;
frequencyOCC = 1-2*mod((0:6*numRB-1).',2)*mod(ports,2);
timeOCC = [ones(1,rankValue); 1-2*(ports>=4)];
toneOffset = mod(floor(ports/2),2);
for symbol = 1:2
    slotSymbol = double(startSymbol)+symbol-1;
    cInit = mod(2^17*(14*double(slotInFrame)+slotSymbol+1)*(2*nid+1)+2*nid+nscid,2^31);
    bits = goldSequence(cInit,2*lastSequenceSymbol);
    bits = bits(2*firstSequenceSymbol+1:end);
    sequence = complex(1-2*bits(1:2:end),1-2*bits(2:2:end)).' / sqrt(2);
    for layer = 1:rankValue
        tones = 1+toneOffset(layer):2:numTones;
        pilot(tones,symbol,layer) = pilotAmplitude*timeOCC(symbol,layer) ...
            *sequence.*frequencyOCC(:,layer);
    end
end
end

function bits = goldSequence(cInit,numBits)
offset = 1600;
stateLength = offset+numBits;
x1 = zeros(1,stateLength);
x2 = x1;
x1(1) = 1;
x2(1:31) = bitget(uint32(cInit),1:31);
for k = 1:stateLength-31
    x1(k+31) = mod(x1(k+3)+x1(k),2);
    x2(k+31) = mod(x2(k+3)+x2(k+2)+x2(k+1)+x2(k),2);
end
bits = mod(x1(offset+1:end)+x2(offset+1:end),2);
end
