function sinr = ul4tPostSinr(Hhat,H,noiseVar)
%UL4TPOSTSINR Score a DMRS-based MMSE receiver over all tones in batches.
% Hhat and H are receive-port x tone x layer. The true channel is used only
% for scoring, while receiver weights are derived exclusively from Hhat.
% This is algebraically identical to solving each tone separately.
[~,F,r] = size(H);
if r==1
    pHat = sum(abs(Hhat).^2,1);
    actual = sum(conj(Hhat).*H,1);
    sinr = pHat.^2./max(realmin,abs(actual-pHat).^2+noiseVar*pHat);
    return
end
A = permute(Hhat,[1 3 2]);
B = permute(H,[1 3 2]);
J = pagemtimes(A,'ctranspose',A,'none');
C = pagemldivide(J+noiseVar*eye(r),pagectranspose(A));
expectedProduct = pagemtimes(C,A);
distortion = pagemtimes(C,B);
diagonalIndices = (1:r+1:r*r).'+r*r*(0:F-1);
expected = reshape(expectedProduct(diagonalIndices),r,F);
distortion(diagonalIndices) = distortion(diagonalIndices)-expected;
residual = reshape(sum(abs(distortion).^2,2)+ ...
    noiseVar*sum(abs(C).^2,2),r,F);
sinr = abs(expected).^2./max(realmin,residual);
end
