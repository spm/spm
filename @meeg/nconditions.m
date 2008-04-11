function res = nconditions(obj)
% Method for getting the number of unique conditions in the file
% FORMAT res = nconditions(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: nconditions.m 1373 2008-04-11 14:24:03Z spm $

res = length(unique(conditions(obj)));
