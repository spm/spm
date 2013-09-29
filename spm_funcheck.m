function [h] = spm_funcheck(f)
% cconverts strings and inline objects to function handles
% FORMAT [h] = spm_funcheck(f)
%
% f   - filename, character expression or inline function
% h   - corresponding function handle
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_speye.m 1131 2008-02-06 11:17:09Z spm $


% create function handle
%==========================================================================

% if f is already a function handle
%--------------------------------------------------------------------------
if isa(f,'function_handle')
    h     = f;
    
% if f is filename or expression
%--------------------------------------------------------------------------
elseif isa(f,'char')
    if exist(f,'builtin') || exist(f,'file');
        h = str2func(f);
    else
        h = spm_funcheck(inline(f));
    end
    
% if f is an inline object
%--------------------------------------------------------------------------
elseif isa(f,'inline')
    names = argnames(f);
    args  = names{1};
    for i = 2:length(names)
        args = [args ',' names{i}];
    end
    h     = eval(['@(' args ')' char(f)]);
end



