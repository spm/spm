function res = nconditions(obj)
% Method for getting the number of unique conditions in the file
% FORMAT res = nconditions(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: nconditions.m 1125 2008-01-30 12:12:18Z vladimir $

res = length(conditions(obj));