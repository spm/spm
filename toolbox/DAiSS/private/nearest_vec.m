function [ index, vec, dist ] = nearest_vec( vec_array, vec_to_find )

% locate bilateral coordinate
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging

% Mark Woolrich
% $Id: nearest_vec.m 8307 2022-08-26 11:00:54Z george $

%--------------------------------------------------------------------------

dists=(sqrt(sum((vec_array-repmat(vec_to_find,size(vec_array,1),1)).^2,2)));

[dist,index]=min(dists);

vec=vec_array(index,:);

end
