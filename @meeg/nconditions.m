function res = nconditions(this)
% Method for getting the number of unique conditions in the file
% FORMAT res = nconditions(obj)
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: nconditions.m 5025 2012-10-31 14:44:13Z vladimir $

res = numel(condlist(this));
