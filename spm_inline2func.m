function [h] = spm_inline2func(f)
% Convert an inline object to a function handle
% FORMAT [h] = spm_inline2func(f)
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2013-2022 Wellcome Centre for Human Neuroimaging


% input argument list
%--------------------------------------------------------------------------
names = argnames(f);
args  = names{1};
for i = 2:numel(names)
    args = [args ',' names{i}];
end
h     = eval(['@(' args ')' formula(f)]);
