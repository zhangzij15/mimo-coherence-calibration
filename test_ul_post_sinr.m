function tests = test_ul_post_sinr
%TEST_UL_POST_SINR Numerical equivalence of the batched DMRS receiver score.
tests = functiontests(localfunctions);
end

function testHandCalculatedDiagonalChannels(t)
H = reshape(diag([2 1]),2,1,2);
verifyEqual(t,ul4tPostSinr(H,H,2),[2;.5],'AbsTol',1e-13);
estimate = reshape(eye(2),2,1,2);
actual = reshape(diag([1.5 .5]),2,1,2);
verifyEqual(t,ul4tPostSinr(estimate,actual,.25),[2;2],'AbsTol',1e-13);
verifyEqual(t,ul4tPostSinr(zeros(2,1,2),zeros(2,1,2),1),zeros(2,1));
end

function testComplexChannelsMatchPerToneEquations(t)
previous = rng;
cleanup = onCleanup(@()rng(previous)); %#ok<NASGU>
rng(1042,'twister');
for rankValue = [1 2 4]
    for tones = [1 7]
        actual = randn(8,tones,rankValue)+1i*randn(8,tones,rankValue);
        estimate = actual+.2*(randn(size(actual))+1i*randn(size(actual)));
        for noiseVar = [.01 3]
            expected = scalarReference(estimate,actual,noiseVar);
            observed = ul4tPostSinr(estimate,actual,noiseVar);
            verifyEqual(t,observed,expected,'RelTol',1e-11,'AbsTol',1e-12);
        end
    end
end
end

function sinr = scalarReference(estimate,actual,noiseVar)
[nReceive,nTones,nLayers] = size(actual);
sinr = zeros(nLayers,nTones);
for f = 1:nTones
    A = reshape(estimate(:,f,:),nReceive,nLayers);
    B = reshape(actual(:,f,:),nReceive,nLayers);
    C = (A'*A+noiseVar*eye(nLayers))\A';
    desired = diag(C*A);
    distortion = C*B-diag(desired);
    residual = sum(abs(distortion).^2,2)+noiseVar*sum(abs(C).^2,2);
    sinr(:,f) = abs(desired).^2./max(realmin,residual);
end
end
