function [data, true_a, true_b, true_c, true_d, true_e, true_t, x, h, N] = toy_data_2D(a,b,c,d,e,  x0, h0, t0, N )
%% Create synthetic data
if nargin==0
    %Form a surface with these CRS parameters
    a = 0.1;
    b = .02;
    c = -.03;
    d = -.2;
    e = .02;
    
    %Size of volume
    N = 25;
    
    %Apex
    x0 = 0;
    h0 = 0;
    t0 = 100;
    
end

%grid
[x,h] = ndgrid(-N:N, -N:N);



%deltas
dx = x-x0;
dh = h-h0;

%CRS traveltime
t_CRS = sqrt(   (t0+ a.*dx + d.*dh).^2 + t0.*b.*dx.^2 + t0.*c.*dh.^2 + 2.*t0.*e.*dx.*dh);

%Display
%clf
%surf(x,h,t_CRS)
%axis equal

%Fit polynomial to surface to obtain analytical CRS parameters
f1 = fit([x(:),h(:)],t_CRS(:),'poly22');

syms x_ h_
f = f1.p00 + f1.p10*x_ + f1.p01*h_ + f1.p20*x_^2 + f1.p11*x_*h_ + f1.p02*h_^2;
dfdx = diff(f,x_);
dfdh = diff(f,h_);
dfdx2 = diff(diff(f,x_),x_);
dfdh2 = diff(diff(f,h_),h_);
dfdhx = diff(diff(f,h_),x_);

true_a = feval(matlabFunction(dfdx),x,h);
true_d = feval(matlabFunction(dfdh),x,h);
true_c = feval(matlabFunction(dfdh2));
true_b = feval(matlabFunction(dfdx2));
true_e = feval(matlabFunction(dfdhx));


%Make synthetic data with gauss wavelet:
min_CRS_t = round(min(t_CRS(:)));
max_CRS_t = round(max(t_CRS(:)));

data = zeros(round(max_CRS_t - min_CRS_t+30), 2.*N+1, 2.*N+1);

g=gausswinN([7,7,7]);

t_CRS = t_CRS - min_CRS_t + 10;
for i = 1:size(t_CRS,1);
    for j = 1:size(t_CRS,2);
        t = t_CRS(i,j);
        
        for i_ = 1:size(g,3)
            for j_ = 1:size(g,2)
                for k_ = 1:size(g,1)
                    
                    di = i+i_-ceil(size(g,3)/2);
                    dj = j+j_-ceil(size(g,2)/2);
                    dk = t+k_-ceil(size(g,1)/2);
                    
                    if di >0 && di < size(data,2);
                        if dj >0 && dj < size(data,3);
                            if dk >1 && dk < size(data,1)-1;
                                
                                t_low = floor(dk);
                                t_high = ceil(dk);
                                t_weight=t_high-dk;
                                
                                data(t_low,di,dj) = data(t_low,di,dj) + g(k_,i_,j_).*t_weight;
                                data(t_high,di,dj) = data(t_high,di,dj) + g(k_,i_,j_).*(1-t_weight);
                            end
                        end
                    end
                end
                
            end
        end
        
    end
end

true_t = round(t_CRS);



