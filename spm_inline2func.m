function [h] = spm_inline2func(f)
% Convert an inline object to a function handle
% FORMAT [h] = spm_inline2func(f)
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_inline2func.m 5664 2013-10-01 18:39:05Z spm $


% input argument list
%--------------------------------------------------------------------------
names = argnames(f);
args  = names{1};
for i = 2:numel(names)
    args = [args ',' names{i}];
end
h     = eval(['@(' args ')' formula(f)]);
