%% Gradient tensor 2D
%==========================================================================
% Returns the gradient tensor for each point in the data
% Created by Anders, Jul 4. 2016
%==========================================================================

function [T,xhy,xhx] = getGradientTensor2D(data, preSmooth, postSmooth)
%Make gradient filters
grad = [-1,0,1];


%Apply pre-smoothing
if numel(postSmooth) > 2
    %User has provided window function
    preSmoothWin = preSmooth;
else
    %Construct smoothing window from parameters using gaussian window
    if length(preSmooth) == 1
        preSmooth = repmat(preSmooth,2);
    end
    
    % One dimention
    preSmoothWin = gausswin( preSmooth(1)*2+1 );
    % Two dimentions
    preSmoothWin = preSmoothWin * gausswin( preSmooth(2)*2+1 )';
    % Normalize
    preSmoothWin = preSmoothWin ./ sum(abs( preSmoothWin(:) ));
end

if preSmooth == 0;
    preSmoothWin = zeros(3,3);
    preSmoothWin(2,2) = 1;
end

hy =  convn(preSmoothWin,grad(:));
hx =  convn(preSmoothWin,grad(:)');


%Normalize filters
hx = hx./sum(abs(hx(:)));
hy = hy./sum(abs(hy(:)));

%Do convolution 

%fprintf('%d/%d Computing gradient x \n',currentN,totalN); currentN = currentN+1;
xhx = convn(data,hx,'same');
%fprintf('%d/%d Computing gradient y \n',currentN,totalN); currentN = currentN+1;
xhy = convn(data,hy,'same');


%Make tensor elements
y2 = xhy.^2;
x2 = xhx.^2;
xy = xhx.*xhy;


%Smooth tensor elements
if numel(postSmooth) > 2
    %User has provided window function
    postSmoothWin = postSmooth;
else
    %Construct smoothing window from parameters using gaussian window
    if length(postSmooth) == 1
        postSmooth = repmat(postSmooth,2);
    end
    
    % One dimention
    postSmoothWin = gausswin( postSmooth(1)*2+1 );
    % Two dimentions
    postSmoothWin = postSmoothWin * gausswin( postSmooth(2)*2+1 )';
    % Normalize
    postSmoothWin = postSmoothWin ./ sum(abs( postSmoothWin(:) ));
end

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