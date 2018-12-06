%% Pgradient tensor 2D eigen slopes
%==========================================================================
% Returns the slopes similr Plane Wave Destructors - but here the
% eigenvalue method is used instead of PWD
% Created by Anders, Jul 6. 2016
%==========================================================================

function [A] = slopes_GT_2D(stk, gradSmooth, tensorSmooth, dt, dx)


%Possible improvement: no weighting function is applied in article
tensorSmooth = gausswinN([tensorSmooth*2+1, tensorSmooth*2+1]);
gradSmooth = gausswinN([gradSmooth*2+1,gradSmooth*2+1]);


%Gradient tensor function gives the needed sums
[T] = gradientTensor2D(stk,gradSmooth,tensorSmooth);

%Eigenvalue decomposition
[e1,e2,v1,v2] = eigenvalueDecomposition2x2(T,1);

%Slope
vt = v1(:,:,1);
vx = v1(:,:,2);
A =  dt/dx .*  vt ./ vx;

end

