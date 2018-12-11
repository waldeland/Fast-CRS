function [T,xhy,xhx] = getGradientTensor2D(data, preSmoothWin, postSmoothWin)
%Make gradient filters
grad = [-1,0,1];

%Apply pre-smoothing
hy =  convn(preSmoothWin,grad(:));
hx =  convn(preSmoothWin,grad(:)');

%Normalize filters
hx = hx./sum(abs(hx(:)));
hy = hy./sum(abs(hy(:)));

%Do convolution 
xhx = convn(data,hx,'same');
xhy = convn(data,hy,'same');

%Make tensor elements
y2 = xhy.^2;
x2 = xhx.^2;
xy = xhx.*xhy;

%Smooth tensor elements
y2 = convn(y2,postSmoothWin,'same');
x2 = convn(x2,postSmoothWin,'same');
xy = convn(xy,postSmoothWin,'same');

%Return tensor
A = cat(3,y2,xy);
B = cat(3,xy,x2);
clear x2 y2 xy

T = cat(4,A,B);
clear A B 

%Remve NaNs
T(isnan(T)) = 0;

end