function [h] = spm_funcheck(f)
% Convert strings and inline objects to function handles
% FORMAT [h] = spm_funcheck(f)
%
% f   - filename, character expression or inline function
% h   - corresponding function handle
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2013-2022 Wellcome Centre for Human Neuroimaging


%-Create function handle
%==========================================================================

% if f is already a function handle
%--------------------------------------------------------------------------
if isa(f,'function_handle') || isempty(f)
    h     = f;
    
% if f is filename or expression
%--------------------------------------------------------------------------
elseif isa(f,'char')
    if exist(f,'builtin') || exist(f,'file')
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
    h     = eval(['@(' args ')' formula(f)]);
end
