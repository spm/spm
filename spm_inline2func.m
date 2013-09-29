function [h] = spm_inline2func(f)
% cconverts an inline object to a function handle
% FORMAT [h] = spm_inline2func(f)
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_speye.m 1131 2008-02-06 11:17:09Z spm $


% input argument list
%--------------------------------------------------------------------------
names = argnames(f);
args  = names{1};
for i = 2:length(names)
    arglist = [args ',' names{i}];
end
h     = eval(['@(' args ')' char(f)]);


