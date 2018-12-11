function [slopes, curvatures, coherency] = fastCRS(data, sigma_g, sigma_T, d)
%% Description

%Input data should always have temporal axis first. For example:
% T x X  - 2D post stack CRS
% T x X x H  - 2D pre stack CRS
% T x X x Y  - 3D post stack CRS
% T x X x Y x Hx x Hy - 3D pre stack CRS


%Output slope will be on the form Nt x Nx x ... x Ny x Ns
%Output curvatures will be on the form Nt x Nx x Ny x Ns x Ns


%%%%%% EXAMPLE USAGE %%%%%%
%dt = 1; dx = 1; dy = 1;
%grad_size = 1; % corresponds to windowsize ~ 7 pixels  (1*3*2 + 1) with gaussian smoothing function
%smooth_size=5; % corresponds to windowsize ~ 31 pixels (5*3*2 + 1) with gaussian smoothing function

%[slopes, curvatures, coherency] = fastCRS(data, grad_size, smooth_size, [dt,dx,dy]); 

%For 2D (input on the form T x X ):
%dt/dx = slopes
%d2t/dx2 = curvatures

%For 3D (input on the form T x X x Y ):
%dt/dx = slopes(:,:,:,1)
%dt/dy = slopes(:,:,:,2)
%d2t/dx2 = curvatures(:,:,:,1,1)
%d2t/dy2 = curvatures(:,:,:,2,2)
%d2t/dxdy = curvatures(:,:,:,1,2) / = curvatures(:,:,:,2,1)

                                %    1   2   3    4
%For 5D (input on the form T x X x Y x Hx x Hy ):
%dt/dx = slopes(:,:,:,1)
%dt/dy = slopes(:,:,:,2)
%dt/dhx = slopes(:,:,:,3)
%dt/dhy = slopes(:,:,:,4)
%d2t/dx2 = curvatures(:,:,:,1,1)
%d2t/dy2 = curvatures(:,:,:,2,2)
%d2t/dhxdx = curvatures(:,:,:,3,1)
%d2t/dydhy = curvatures(:,:,:,2,4)
% etc.... 


%%
%Determine dimension of data
global Ns dim N
dim = numel(size(data));
Ns = size(data); %Number of pixels pr dimension
N = prod(Ns); %Number of pixels in total

%Default arguments
if ~exist('d','var'); d=ones(dim,1); end
if ~exist('sigma_g','var'); sigma_g=1; end
if ~exist('sigma_T','var'); sigma_T=5; end

%Ensure that d is column
d = reshape(d,numel(d),1);

% Get tensor fields
[T,xT,XT]  = quadraticGradientTensor(data, sigma_g, sigma_T, nargout~=1);

% Since all the following operations are the same (and independent) for all
% spatial coordinates, we reshape the tensors such that the spatial
% dimmensions are represented in only one dimention - the last one. This
% makes it trival to handle ND images
T  = reshape(T, dim,dim,N);
if nargout~=1
    xT = reshape(xT,dim,dim,dim,N);
    XT = reshape(XT,dim,dim,dim,dim,N);
end

%% Get slopes
% Do eigenvalue decomposition of GST
[e,v] = eigenvalue_decomposition(T);

%Make sure that t-component in eigenvector 1 is positive (to avoid phase jumps)
v(:,1) = standardize_eigenvectors(v(:,1),1);
% and than the second component in the other vectors is positve
v(:,2:end) = standardize_eigenvectors(v(:,2:end),2);

%Compute GST coherency 
coherency = e(1,:)./sum(e,1); %It is possible to derive coherency from QST as well 

%Compute slopes (dx/dt, dy/dt etc.)
slopes = squeeze( bsxfun(@times, - bsxfun(@times, v(2:end,1,:), 1./v(1,1,:)) , d(1)./d(2:end) ) );

%... hack to fix Matlab wierdness ...
if dim == 2; slopes = slopes'; end

%If the user only need the slopes
if nargout==1;
    slopes = permute(slopes, [ 2, 1]);
    slopes = reshape(slopes, cat(2, Ns, dim-1) ); 
    return; 
end

%% Define surface normal vector 
%Remove physical dimentions from slopes
q = bsxfun(@times,slopes, d(2:end)./d(1));

% Matricies with 1 and zero    
one = ones(1,N);
zero = zeros(1,N);

% Surface normal vector has one in t-component and slopes in the other components
eu = cat(1, one, -q);
eu = normalize_vec(eu);


%% Extract curvatures

%Matrix for all curvatures
curvatures = zeros(dim-1,dim-1,N);

%Get all second order derivative (not mixed derivatives)
evs = {};
for i = 1:dim-1
    
    %Create second unit vector for this curvature (see Waldeland et al. 2018 - Geophysics)
    before = repmat(zero,i-1,1); %number of zeros before the one
    after = repmat(zero,dim-i-1,1); %number of zeros after the one
    ev  = [q(i,:); before; one; after;];
    ev  = normalize_vec(ev);
    
    %Save ev to construct rotated coordinate system later
    evs{i} = ev;
    
    %Get curvature with units
    curvatures(i,i,:) = extract_curvature(xT,XT,eu,ev) .* d(1)./d(i+1).^2;
end

%Mixed derivatives (rotated coordinate systems (vec(o), vec(p) etc) )
for i = 1:dim-1
    for j = i+1:dim-1
        ev = evs{i};
        ew = evs{j};
                
        eo = normalize_vec( ev+ew );
        ep = normalize_vec( ev-ew );
        
        ko =  extract_curvature(xT,XT,eu,eo);
        kp =  extract_curvature(xT,XT,eu,ep);
        
        curvatures(i,j,:) = -(kp-ko)/2 .* d(1)/( d(i)*d(j));
        
        %Symmetric
        curvatures(j,i,:) = curvatures(i,j,:);
        
    end
end

%% Post process arrays

%Make data axes the first dimensions
curvatures = permute(curvatures, [3, 1, 2]);
slopes = permute(slopes, [ 2, 1]);

%Reshape to fit image dimensions
curvatures = reshape(curvatures, cat(2, Ns, dim-1,dim-1) );
coherency = reshape(coherency, cat(2, Ns, 1) );
slopes = reshape(slopes, cat(2, Ns, dim-1) );



%l1 = fufu ;
%l2 = fvfv-vfufv.^2/vvfufu;
% coh(iy,ix) = (l1)./(l1+l2);


end

%% Extract curvature 
%this function implements the equations in Bakker 2002
function kappa = extract_curvature(xT,XT,eu,ev)
global Ns dim

    %  ev' . (XT o eu x ev') ;
    nom =  dot( ev, circ(xT, cross(eu, ev) ) );
    %  ev' . (XT o eu x eu') . ev;
    denom =   dot( dot( eu, circ(XT, cross(ev, ev))), eu );
    
    kv = (nom./denom);
    
    %Rotation
    q = get_slope(ev);
    kappa = kv.*(1+q.^2).^(3/2);
    
end

%% Ensure that eigenvectors has positive component along axis
function v = standardize_eigenvectors(v,axis)
    v = bsxfun(@times, v, sign(v(axis,:)) );
end

%% Eigenvaule decompostion
function [e,v] = eigenvalue_decomposition(T)
global N dim

%If 3D we use a vectorized faster implementation
if dim == 3;
    [e1,e2,e3,v1,v2,v3] = fastEig3D(T,1) ;
    e = cat(1,e1,e2,e3);
    v = cat(2,reshape(v1,3,1,N) ,reshape(v2,3,1,N), reshape(v3,3,1,N));
    return
end

e = zeros(dim,N);
v = zeros(dim,dim,N);

for i=1:N;
    %Eigenvalue decomp
    [vectors,values] = eig(T(:,:,i));
    values = diag(values);
    
    %Sort on eigenvalues
    [values,sorting] = sort(values,'descend');
    vectors = vectors(:,sorting);
    
    %Insert into arrays
    e(:,i) = values; 
    v(:,:,i) = vectors;
end
end

%% vectorized calculation of cross product of vectors
function A = cross(b,c)
global N dim
A = zeros(cat(2,dim,dim,N));

for k = 1:size(b,1)
    for l = 1:size(c,1)
        A(k,l,:) = b(k,:).*c(l,:);
    end
end
end

%% vectorized calculation of dot product between 1d or 2d matrix and 1d matrix
function A = dot(C_,b_)
global N dim
%Matrix first
if numel(C_) > numel(b_)
    C = C_;
    b = b_;
else
    b = C_;
    C = b_;
end

% If  C is vector (assuming b is vector)
if numel(size(C)) == 2
    A = zeros(1,N);
    
    for l = 1:dim
        A(1,:) = reshape( A(1,:) +  b(l,:) .* C(l,:) , 1, prod(N));
    end
    return
end

% If C is matrix (assuming b is vector)
if numel(size(C)) == 3
    
    A = zeros(dim,N);
    for k = 1:dim
        for l = 1:dim
            A(k,:) =  squeez(A(k,:)) +  squeez(b(l,:)) .* squeez(C(l,k,:)) ;
        end
    end
    
    return
end

error('Error in dot-function')
end

%% vectorized calculation of the circ operator
function A = circ(B,C)
global N dim


% Use this for XT circ ...
if numel(size(B)) == 4 + 1 
    
    A = zeros(cat(2,dim,dim,N));
    
    for k = 1:dim
        for l = 1:dim
            
            for i = 1:dim
                for j = 1:dim
                    A(k,l,:) = squeez(A(k,l,:)) +   squeez(B(i,j,k,l,:)).* squeez(C(i,j,:)) ;
                end
            end
            
        end
    end
    
else % And this for xT circ ...
 
    A = zeros(dim,N);
    
    for l = 1:dim
        A(l,:) = 0;
        for i = 1:dim
            for k = 1:dim
                A(l,:) = squeez(A(l,:)) + squeez(B(l,k,i,:)).*squeez(C(i,k,:));
            end
        end
    end
    
end
end

%% Normalize to make unit vector
function vec = normalize_vec(vec);
 vec = bsxfun(@times,vec, sqrt(sum(vec.^2,1)));
end

%% Get slope from unit vector
function q = get_slope(ev)
q = ev(1,:)./sqrt(sum( ev(2:end,:).^2, 1));
end

%% True squeeze function, since MATLAB consider both column and row vectors for 1D
function s = squeez(s)
s = squeeze(s);
if numel(size(s)) == 2
    s = s(:);
end
end
