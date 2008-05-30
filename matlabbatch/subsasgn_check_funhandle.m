function sts = subsasgn_check_funhandle(val)

% function sts = subsasgn_check_funhandle(val)
% Return true if val is either empty, or a function or function handle.
% One could also check for nargin == 1 and nargout == 1.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_check_funhandle.m 1764 2008-05-30 13:09:40Z volkmar $

rev = '$Rev: 1764 $'; %#ok

sts = isempty(val) || isa(val, 'function_handle') || ...
    (ischar(val) && (any(exist(val) == 2:6) || ~isempty(which(val))));
if sts && isa(val, 'function_handle')
    try
        f = functions(val);
    catch
        % fail silently if "functions" evaluation fails, keep sts == true
        return;
    end;
    try
        if ~strcmp(f.type,'anonymous')
            % file name should not be empty for file functions
            sts = ~isempty(f.file) || ~isempty(which(f.function));
        end;
    catch
        % fail silently if something is wrong with f, keep sts == true
        return;
    end;
end;
        
        
        