function xy = coor2D(obj, ind)
% returns x and y coordinates of channels in 2D plane
% FORMAT coor2D(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel

% this is a preliminary version 

if isfield(obj.channels, 'X_plot2D')
    x = cat(2, obj.channels.X_plot2D);
    y = cat(2, obj.channels.Y_plot2D);
    
    xy = [x; y];
else 
    xy = [];
end

if nargin > 1
    xy = xy(:, ind);
end