function res = nconditions(obj)
% Method for getting the number of unique conditions in the file
% FORMAT res = nconditions(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id $

res = length(conditions(obj));