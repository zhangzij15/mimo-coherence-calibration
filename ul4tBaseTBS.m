function tbs = ul4tBaseTBS(qm,nLayers,nRB,nRePerRB,codeRate)
%UL4TBASETBS NR transport block size for a single codeword.
%   QM and CODERATE are same-size arrays of bits per symbol and target rates.
%   NREPRB is the data RE budget per RB after DM-RS and other overhead.
%   See TS 38.214 sections 5.1.3.2 and 6.1.4.2.
validateattributes(qm,{'numeric'},{'real','finite','nonempty','positive','integer'});
validateattributes(nLayers,{'numeric'},{'scalar','integer','>=',1,'<=',4});
validateattributes(nRB,{'numeric'},{'scalar','real','finite','nonnegative','integer'});
validateattributes(nRePerRB,{'numeric'},{'scalar','real','finite','nonnegative','integer'});
validateattributes(codeRate,{'numeric'},{'real','finite','>=',0,'<=',1});
assert(isequal(size(qm),size(codeRate)),'UL4T:BaseTBSShape', ...
    'Modulation orders and target code rates must have matching sizes.');

numRE = min(156,double(nRePerRB))*double(nRB);
if numRE==0
    tbs = zeros(size(qm));
    return;
end
rates = double(codeRate);
informationBits = numRE*double(nLayers)*double(qm).*rates;
tbs = arrayfun(@transportBlockSize,informationBits,rates);
end

function tbs = transportBlockSize(informationBits,codeRate)
smallTBS = [24 32 40 48 56 64 72 80 88 96 104 112 120 128 136 144 152 160 ...
    168 176 184 192 208 224 240 256 272 288 304 320 336 352 368 384 408 ...
    432 456 480 504 528 552 576 608 640 672 704 736 768 808 848 888 ...
    928 984 1032 1064 1128 1160 1192 1224 1256 1288 1320 1352 1416 ...
    1480 1544 1608 1672 1736 1800 1864 1928 2024 2088 2152 2216 2280 ...
    2408 2472 2536 2600 2664 2728 2792 2856 2976 3104 3240 3368 3496 3624 3752 3824];
if informationBits<=3824
    step = 2^max(3,floor(log2(informationBits))-6);
    quantized = max(24,step*floor(informationBits/step));
    tbs = smallTBS(find(smallTBS>=quantized,1));
    return;
end

step = 2^(floor(log2(informationBits-24))-5);
quantized = max(3840,step*round((informationBits-24)/step));
if codeRate<=1/4
    blocks = ceil((quantized+24)/3816);
elseif quantized>8424
    blocks = ceil((quantized+24)/8424);
else
    blocks = 1;
end
tbs = 8*blocks*ceil((quantized+24)/(8*blocks))-24;
end
