function res = nconditions(this)
% Method for getting the number of unique conditions in the file
% FORMAT res = nconditions(obj)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


res = numel(condlist(this));
