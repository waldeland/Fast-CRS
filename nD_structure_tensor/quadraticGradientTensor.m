function [T,xT,XT]  = quadraticGradientTensor(data, sigma_g, sigma_T, compute_QST)
if ~exist('compute_QST','var'); compute_QST = true; end; %Switch to not compute xT and XT if only GST is needed
break_of_sigma = 3.; %number og std with of gaussian functions


%Determine dimension of data
dim = numel(size(data));
Ns = size(data); %Number of pixels/voxels/...



%Make sigma_g
if numel(sigma_g) == 1 %User has selected same size for all window sides
    sigma_g = repmat(sigma_g,1,dim);
end

%Make sigma_T
if numel(sigma_T) == 1 %User has selected same size for all window sides
    sigma_T = repmat(sigma_T,1,dim);
end



%% Make gradient filters
filtersize = ceil(break_of_sigma*sigma_g);

%Make components of filter
grad_smooth_filt = {};
for d = 1:dim
    %Make grid
    grid = -filtersize(d):filtersize(d); 
    %Reshape to correct dimension
    grid = reshape(grid, [repmat(1,1,d-1),length(grid), repmat(1,1,dim-d)]);
    %Gaussian function
    win =  1/(2 * pi * sigma_g(d).^2)  * exp(grid.^2 / (-2 * sigma_g(d).^2));
    grad_smooth_filt{d} = win;
end

%Add derivative component and filter to get gradients
gs = {}; %gradient images
for d = 1:dim
    %Make grid
    grid = -filtersize(d):filtersize(d); 
    %Reshape to correct dimension
    grid = reshape(grid, [repmat(1,1,d-1),length(grid), repmat(1,1,dim-d)]);
    
    %Make copy of smoothing window
    filter = grad_smooth_filt;
    
    %Get derivative of given dimention
    filter{d} =  1/(sigma_g(d)^2)* grid .* filter{d} ;
    
    gs{d} = fast_conv(data,filter); 
    
    ndgrid_input{d} = -filtersize(d):filtersize(d);
end


%% Grids and smoothing window for QST/GST
filtersize = break_of_sigma*sigma_T;
smooth_win = {};
grids = {};
for d = 1:dim
    %Make grid
    grid = -filtersize(d):filtersize(d);
    %Reshape to correct dimension
    grid = reshape(grid, [repmat(1,1,d-1),length(grid), repmat(1,1,dim-d)]);
    grids{d} = grid;
    %Gaussian window
    win = 1/(2 * pi * sigma_T(d).^2) * exp(grid.^2/(-2 * sigma_T(d) * sigma_T(d)));
    smooth_win{d} = win;
end


%% Output arrays (we put the i,j,k,l dimentions first, then the spatial dimentions contained in Ns-variable )
                % k,l
T  = zeros(cat(2,dim,dim,Ns),'single');

if nargout>1 && compute_QST
                        %i,j,k,l
    XT = zeros(cat(2,dim,dim,dim,dim,Ns));
                        %i,k,l
    xT = zeros(cat(2,dim,dim,dim,Ns));
else
    xT = [];
    XT = [];
end

%Compute tensor elements
ind = 1;
for k = 1:dim
    for l = k:dim
        
        %Gradient term
        FkFl = gs{k}.*gs{l};
        
        % GST
        out =  fast_conv(FkFl, smooth_win,'same');
        T(k,l,:) = reshape(out,1,1,numel(out));
        T(l,k,:) = reshape(out,1,1,numel(out));
            
        if nargout>1 && compute_QST
            for i = 1:dim
                filter = smooth_win;
                filter{i} = filter{i}.*grids{i};
                
                %Linear elements
                out =  fast_conv(FkFl, filter,'same');
                xT(i,k,l,:) = reshape(out,1,1,1,numel(out));
                xT(i,l,k,:) = reshape(out,1,1,1,numel(out));
                
                %Quadratic elements
                for j = i:dim
                    filter = smooth_win;
                    filter{i} = filter{i}.*grids{i};
                    filter{j} = filter{j}.*grids{j};
                
                    out =  fast_conv(FkFl, filter,'same');
                    XT(i,j,k,l,:) = reshape(out,1,1,1,1,numel(out));
                    XT(i,j,l,k,:) = reshape(out,1,1,1,1,numel(out));
                    XT(j,i,k,l,:) = reshape(out,1,1,1,1,numel(out));
                    XT(j,i,l,k,:) = reshape(out,1,1,1,1,numel(out));
                end
                
            end
        end
        
        fprintf(['Created tensor element ' num2str(ind) ' of ' num2str(sum(1:dim)) char(10)])
        ind = ind +1 ;
    end
end
%%


end

%Compute with separable convolutions (assuming f is struct with 1D filters)
function y = fast_conv(x,f, type)
if nargin == 2; type = 'same'; end
y = x;
for i = 1:length(f)
    y = convn(y,f{i},type);
end

end