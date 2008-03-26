function res = nconditions(obj)
% Method for getting the number of unique conditions in the file
% FORMAT res = nconditions(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: nconditions.m 1236 2008-03-20 18:15:33Z stefan $

res = size(unique(conditions(obj), 'rows'),1);