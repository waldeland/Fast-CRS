%% imlim
%==========================================================================
% finds the limits for a imagesc plot based on percentile
% Created by Anders, Jul 5. 2016
%==========================================================================
function lim = imlim(data,percentage,sym)
if nargin <2
    percentage = 1;
end

if 100-percentage < percentage
    percentage = 100-percentage;
end
%Remove Nans
data = data(:);
data = data(~isnan(data));
data = data(~isinf(data));


lim = [prctile(data(:),percentage), prctile(data(:),100-percentage)];

for i = 1:2
    if isnan(lim(i))
        lim(i)=0;
    end
end


if lim(2) <= lim(1); 
    lim(2) = lim(2) +  lim(2) -lim(1) + .01;
end

if nargin>2 & sym
    lim = [-max(abs(lim)),max(abs(lim))];
end


