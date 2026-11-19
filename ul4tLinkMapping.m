function varargout = ul4tLinkMapping(action,varargin)
%UL4TLINKMAPPING Effective SINR and AWGN threshold selection for TB BLER 0.1.
switch action
    case 'prepare'
        varargout{1} = prepareTable(varargin{:});
    case 'evaluate'
        [varargout{1:nargout}] = evaluate(varargin{:});
    otherwise
        error('Unknown link mapping action.');
end
end

function table = prepareTable(cfg,groupCount,prefixSymbols)
% TS 38.214 Tables 5.1.3.1-1 and 5.1.3.2-1 and PUSCH section 6.1.4.2.
if nargin<3, prefixSymbols = 0; end
R = [120 157 193 251 308 379 449 526 602 679 340 378 434 490 553 ...
    616 658 438 466 517 567 616 666 719 772 822 873 910 948]/1024;
Qm = [2*ones(1,10),4*ones(1,7),6*ones(1,12)];
beta = [1.6 1.61 1.63 1.65 1.67 1.7 1.73 1.76 1.79 1.82 3.97 4.27 4.71 ...
    5.16 5.66 6.16 6.5 9.95 10.97 12.92 14.96 17.06 19.33 21.85 24.51 ...
    27.14 29.94 32.05 34.28]; % 5G-LENA EESM calibration, linear SINR.
if strcmp(cfg.ul.mapping,'sem')
    beta = cfg.ul.sem_beta(:).';
    if isscalar(beta), beta = repmat(beta,1,numel(R)); end
    assert(numel(beta)==numel(R) && all(isfinite(beta) & beta>0),'Invalid SEM beta.');
else
    assert(strcmp(cfg.ul.mapping,'eesm'),'Mapping must be sem or eesm.');
end
root = fileparts(mfilename('fullpath'));
reference = dlmread(fullfile(root,'awgn_bler10_t1.csv'),',',1,0);
NrePerRB = min(156,12*(14-prefixSymbols)-2*6*groupCount);
table.data_re_per_rb = NrePerRB;
table.bits = ul4tBaseTBS(Qm,cfg.ul.rank,cfg.ofdm.num_rb,NrePerRB,R);
table.required_dB = inf(size(R));
table.reference_cbs = zeros(size(R));
table.bg = zeros(size(R));
table.num_cb = zeros(size(R));
table.beta = beta;
table.mapping = cfg.ul.mapping;
for m = 1:numel(R)
    tbs = table.bits(m);
    [bg,cbs,C] = codeBlockSize(tbs,R(m));
    candidate = find(reference(:,1)==bg & reference(:,2)==m-1);
    assert(~isempty(candidate) && C<=size(reference,2)-3, ...
        'AWGN reference lacks this MCS, base graph or segmentation size.');
    available = reference(candidate,3);
    [~,j] = min(abs(log(available/cbs))); % nearest measured CBS on log-bit scale
    threshold = reference(candidate(j),3+C);
    if threshold<99, table.required_dB(m) = threshold; end
    table.reference_cbs(m) = available(j);
    table.bg(m) = bg;
    table.num_cb(m) = C;
end
assert(any(isfinite(table.required_dB)),'No usable AWGN thresholds.');
end

function [goodBits,mcs,eff_dB] = evaluate(sinr,table)
% One codeword across all occupied data RE/layers. SINRs are linear.
s = max(0,sinr(:));
if strcmp(table.mapping,'sem')
    eff = expm1(mean(log1p(bsxfun(@times,s,table.beta)),1))./table.beta;
else
    floorSINR = min(s);
    eff = floorSINR-table.beta.*log(mean(exp(-bsxfun(@rdivide,s-floorSINR,table.beta)),1));
end
eff_dB_all = 10*log10(max(eff,realmin));
bits = table.bits.*(eff_dB_all>=table.required_dB);
[best,i] = max(bits);
goodBits = .9*best;
mcs = i-1;
eff_dB = eff_dB_all(i);
if best==0, mcs = -1; end
end

function [bg,cbs,C] = codeBlockSize(A,R)
% TS 38.212 sections 6.2.1--6.2.3. CBS includes lifting/filler length.
B = A+16+8*(A>3824);
if A<=292 || R<=1/4 || (A<=3824 && R<=.67)
    bg = 2;
    Kcb = 3840;
    K = 10;
    Kb = 6+2*(B>192)+(B>560)+(B>640);
else
    bg = 1;
    Kcb = 8448;
    K = 22;
    Kb = K;
end
C = 1;
if B>Kcb
    C = ceil(B/(Kcb-24));
    B = B+24*C;
end
z = sort([2*2.^(0:7),3*2.^(0:7),5*2.^(0:6),7*2.^(0:5), ...
    9*2.^(0:5),11*2.^(0:5),13*2.^(0:4),15*2.^(0:4)]);
cbs = K*z(find(Kb*z>=ceil(B/C),1));
end
