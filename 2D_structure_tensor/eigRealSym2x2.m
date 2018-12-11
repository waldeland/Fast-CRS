function [e1,e2,v1,v2] = eigRealSym2x2(T,sign_dim)
% A vectorized implementation of eigenvalue decomposition of a symmetric
% 2x2 matrix
if nargin<1;
    test();
    return
end

%Eigenvalue decomposition 
%https://en.wikipedia.org/wiki/Eigenvalue_algorithm
trace = T(:,:,1,1) + T(:,:,2,2) ;
deter = T(:,:,1,1) .* T(:,:,2,2)    -       T(:,:,1,2) .* T(:,:,2,1);
eA = ( trace + sqrt(trace.^2-4.*deter) )./2;
eB = ( trace - sqrt(trace.^2-4.*deter) )./2;
e1 = max( cat(4,eA,eB), [], 4);
e2 = eA + eB - e1;

% Identity
I = repmat( reshape([1,0;0,1],1,1,2,2), size(T,1), size(T,2) );

%Find V1/V2 - camly hammiton theoremn
M2 = T-bsxfun(@times,e2,I);
v2=M2(:,:,:,1);        

M1 = T-bsxfun(@times,e1,I);
v1=M1(:,:,:,1);

%Make the sign_dim positive
v2_sign = sign(v2(:,:,sign_dim));
v2 = bsxfun(@times,v2,v2_sign);

v1_sign = sign(v1(:,:,sign_dim));
v1 = bsxfun(@times,v1,v1_sign);

%Normalize to 1
norm = sqrt( v2(:,:,1).^2 + v2(:,:,2).^2 );
v2 = bsxfun(@times,v2,1./norm);

norm = sqrt( v1(:,:,1).^2 + v1(:,:,2).^2 );
v1 = bsxfun(@times,v1,1./norm);


end

function test()

%Build tensor
N1 = 100; N2 = 100;
T = randn(N1,N2,2,2);
T(:,:,1,2) = T(:,:,2,1);

%Standard method

% Output arrays
out = T(:,:,1,1).*0;
e1 = out; e2 = out; 
out = T(:,:,1:2).*0;
v1 = out;   v2 = out;


%Eigenvalue decomposition

for i = 1:size(T,1)
    for j = 1:size(T,2)
        
        t = reshape(T(i,j,:,:),2,2);
        [v,e] = eig(t);
        [e,ind] = sort(diag(e));
        e1(i,j) = e(2);
        e2(i,j) = e(1);
        %v = v(:,ind);
        
        if v(1,1) < 0
            v(:,1) = -v(:,1);
        end
        v1(i,j,:) = v(:,1);
        
        if v(1,2) < 0
            v(:,2) = -v(:,2);
        end
        v2(i,j,:) = v(:,2);
 
    end
end
% This method
[e1_,e2_,v1_,v2_] = eigRealSym2x2(T,1);

rms = @(x,y) mean(abs(x(:)-y(:)));
fprintf([ 'Unit test: eigRealSym2x2.m:' char(10)])
%fprintf([ 'verification of eig: ' num2str(rms(ver,ver.*0) ) char(10)])
fprintf([ 'E1 mean error: ' num2str(rms(e1_,e1) ) char(10)])
fprintf([ 'E2 mean error: ' num2str(rms(e2_,e2) ) char(10)])
fprintf([ 'V1 mean error: ' num2str(rms((v1_),(v1)) ) char(10)])
fprintf([ 'V2 mean error: ' num2str(rms((v2_),(v2)) ) char(10)])

end