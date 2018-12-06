function [e1,e2,e3,v1,v2,v3] = fastEig3D(T,sign_ind) 
% Given a real symmetric 3x3 matrix A, compute the eigenvalues and
% eigenvectors
% Algorithm from
% https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
% (11.11.15)
% vectorized implementation by Anders Ueland 2015
if nargin <1
   test()
   return
end

%Index of vector component that should be positive
if nargin < 2
    sign_ind = 3;
end
   
% T is given with size 3 x 3 X N
[~,~,N] = size(T);


p1 = T(1,2,:).^2 + T(1,3,:).^2 + T(2,3,:).^2;

T_is_diagonal = p1==0;
tnd = ~T_is_diagonal; %T Not Diagonal

%if (p1 == 0)
% T is diagonal.
e1(T_is_diagonal) = T(1 ,1, T_is_diagonal);
e2(T_is_diagonal) = T(2, 2, T_is_diagonal);
e3(T_is_diagonal) = T(3, 3, T_is_diagonal);

clear T_is_diagonal

%else
%q = trace(T )/3;
trace3 = ( T(1,1,tnd) + T(2,2,tnd) + T(3,3,tnd) )./3;

p2 = (T(1,1,tnd) - trace3).^2 + (T(2,2,tnd) - trace3).^2 + (T(3,3,tnd) - trace3).^2 + 2 * p1(tnd);
p = sqrt(p2 ./ 6);
%B = (1 / p) * (T - q * I)       % I is the identity matrix
clear p2 p1

I =  repmat( reshape(eye(3,3),3,3) ,1,1,N) ;  
qI = bsxfun(@times,I(:,:,tnd),trace3);
B = bsxfun(@times,1./p,T(:,:,tnd)-qI);
clear qI I

%r = det(B) / 2
%From http://www.dr-lex.be/random/matrix-inv.html
%det(a)  = a11       (a33       a22     -a32       a23     ) - a21       (a33       a12     -a32       a13     ) + a31       (a23       a12     -a22       a13     )
r        = B(1,1,:).*(B(3,3,:).*B(2,2,:)-B(3,2,:).*B(2,3,:)) - B(2,1,:).*(B(3,3,:).*B(1,2,:)-B(3,2,:).*B(1,3,:)) + B(3,1,:).*(B(2,3,:).*B(1,2,:)-B(2,2,:).*B(1,3,:));
r =r./2;

clear B

% In exact arithmetic for a symmetric matrix  -1 <= r <= 1
% but computation error can leave it slightly outside this range.
phi = r.*0;

% if (r <= -1)
%     phi = pi / 3
r_leq_neg1 = r <= -1;
phi(r_leq_neg1) = pi/3;

% elseif (r >= 1)
%     phi = 0
r_geq_1 = r >= 1;
%phi(r_geq_1) = 0; %redundant

% else
%     phi = acos(r) / 3
phi( ~r_leq_neg1 & ~r_geq_1) = acos(r( ~r_leq_neg1 & ~r_geq_1))./3;

clear r_leq_neg1 r_geq_1 r



% the eigenvalues satisfy eig3 <= eig2 <= eig1
e1(tnd) = trace3 + 2 .* p .* cos(phi);
e3(tnd) = trace3 + 2 .* p .* cos(phi + (2.*pi/3));
e2(tnd) = 3 .* squeeze(trace3) - e1(tnd)' - e3(tnd)' ;    % since trace(A) = eig1 + eig2 + eig3

clear q p phi tnd

%Compute eigenvalues:
% v1 = last column in (A-e2I)(A-e3I)

if nargout>3
    v1 = getVa(e2,e3,T,sign_ind, N);
end

if nargout>4
    v2 = getVa(e1,e3,T,sign_ind, N);
end

if nargout>5
    v3 = getVa(e1,e2,T,sign_ind, N);
end





end

function v = getVa(eb,ec,T,sign_ind, N)

    I =  repmat( reshape(eye(3,3),3,3) ,1,1,N) ;    
    
    T_m_ebI = T - bsxfun(@times,I,reshape(eb,1,1,N));
    
    T_m_ecI = T - bsxfun(@times,I,reshape(ec,1,1,N));
    
    %vs = mtimesx(T_m_ebI,T_m_ecI);
    vs = cross(T_m_ebI,T_m_ecI);
    v = vs(1,:,:);
    
    v = reshape(v,3,N);
    eps = 0.0000000;
    v = bsxfun(@times, v, 1./ sqrt( v(1,:).^2 + v(2,:).^2 + v(3,:).^2 +eps) );
    
    %Force z-component to be positive
    flip_sign = double(sign(v(sign_ind,:)) > 0);
    v = bsxfun(@times,v, sign(flip_sign - .5) ); 

end

%Vectorized cross (?) multiplication  (replaces mtimesx-function to avoid dependencies on compiled files)
function C = cross(A,B);

dim = 3;
C = A.*0;

%Loop through output elements
for i = 1:dim;
    for j = 1:dim;
        C(i,j,:) = 0;
        for k = 1:dim
        
                 C(i,j,:) = C(i,j,:) + A(i,k,:).*B(k,j,:);
        end
    end
end
end

function test()

%Build tensor
N1 = 10; N2 = 10; N3 = 10;
T = randn(3,3,N1*N2*N3);
T(1,2,:) = T(2,1,:); T(1,3,:) = T(3,1,:); T(2,3,:) = T(3,2,:);

%Standard method
e1 = zeros(N1*N2*N3,1); e2 = zeros(N1*N2*N3,1); e3 = zeros(N1*N2*N3,1);
v1 = zeros(3,N1*N2*N3); v2 = zeros(3,N1*N2*N3); v3 = zeros(3,N1*N2*N3);
ver = zeros(N1,N2,N3);
for i=1:N1*N2*N3
    
            size(T);
            t = reshape(T(:,:,i),3,3);
            [v,e] = eig(t);
            ver(i) = sum(sum( abs(t*v - v*e) ));
            [e,ind] = sort(diag(e));
            e1(i) = e(3);
            e2(i) = e(2);
            e3(i) = e(1);
            %v = v(:,ind);
            v1(:,i) = v(:,3);
            v2(:,i) = v(:,2);
            v3(:,i) = v(:,1);
            

end

% This method
[e1_,e2_,e3_,v1_,v2_,v3_] = fastEig3D(T);

rms = @(x,y) mean(abs(x(:)-y(:)));
fprintf([ 'Unit test: fastEig3D.m:' char(10)])
fprintf([ 'verification of eig: ' num2str(rms(ver,ver.*0) ) char(10)])
fprintf([ 'E1 mean error: ' num2str(rms(e1_,e1) ) char(10)])
fprintf([ 'E2 mean error: ' num2str(rms(e2_,e2) ) char(10)])
fprintf([ 'E3 mean error: ' num2str(rms(e3_,e3) ) char(10)])
fprintf([ 'V1 mean error (absolute values): ' num2str(rms(abs(v1_),abs(v1)) ) char(10)])
fprintf([ 'V2 mean error (absolute values): ' num2str(rms(abs(v2_),abs(v2)) ) char(10)])
fprintf([ 'V2 mean error (absolute values): ' num2str(rms(abs(v3_),abs(v3)) ) char(10)])

end