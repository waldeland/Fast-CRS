function win = gausswinN(s)
%1D window
win1 = gausswin(s(1)); 
win = win1;

%2D window
if numel(s)>1
    win2 = gausswin(s(2)); 
    win = win1*win2'; 
end

%3D window
if numel(s)>2
    win3 = gausswin(s(3));  
    win = repmat(win,1,1,s(3)); 
    win = bsxfun(@times,win,reshape(win3,1,1,s(3))); 
end

%4D window
if numel(s)>3
    win4 = gausswin(s(4));  
    win = repmat(win,1,1,1,s(4)); 
    win = bsxfun(@times,win,reshape(win4,1,1,1,s(4))); 
end

%5D window
if numel(s)>4
    win5 = gausswin(s(5));  
    win = repmat(win,1,1,1,1,s(5)); 
    win = bsxfun(@times,win,reshape(win5,1,1,1,1,s(5))); 
end

if numel(s)>5
    error('Only up to 5D windows supported. Code can easily be extended though...')
end

%Normalize
win = win./sum(win(:));
end