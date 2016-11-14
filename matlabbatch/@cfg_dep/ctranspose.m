function b = ctranspose(a)

% function b = ctranspose(a)
% Transpose a configuration dependency
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2016 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: ctranspose.m 6929 2016-11-14 13:07:31Z guillaume $

rev = '$Rev: 6929 $'; %#ok

s = size(a);
b = reshape(a,fliplr(s));
