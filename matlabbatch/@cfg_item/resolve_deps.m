function [val sts] = resolve_deps(item, cj)

% function [val sts] = resolve_deps(item, cj)
% Resolve dependencies for an cfg item. This is a generic function that
% returns the contents of item.val{1} if it is an array of cfg_deps. If
% there is more than one dependency, they will be resolved in order of
% appearance. The returned val will be the concatenation of the values of
% all dependencies. A warning will be issued if this concatenation fails
% (which would happen if resolved dependencies contain incompatible
% values).
% If any of the dependencies cannot be resolved, val will be empty and sts
% false.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: resolve_deps.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

try
    val1 = cell(size(item.val{1}));
    for k = 1:numel(item.val{1})
        % Outputs are stored in .jout field of cfg_exbranch, which is
        % neither included in .src_exbranch nor in .src_output substruct
        val1{k} = subsref(cj, [item.val{1}(k).src_exbranch, substruct('.','jout'), item.val{1}(k).src_output]);
    end;
    sts = true;
catch
    % some dependencies are not yet resolved
    val = [];
    sts = false;
    return;
end;
if sts
    % All items resolved, try concatenation
    try
        % try concatenation along 1st dim
        val = cat(1, val1{:});
    catch
        % try concatenation along 2nd dim
        try
            val = cat(2, val1{:});
        catch
            % all concatenations failed, display warning
            warning('matlabbatch:resolve_deps:concat','Dependencies resolved, but incompatible values.');
            l = lasterror;
            fprintf('%s\n',l.message);
            fprintf('In item %s:\n', subsref(item, substruct('.','name')));
            for k = 1:numel(item.val{1})
                fprintf('Dependency %d: %s\n', k, item.val{1}(k).sname);
                disp(val1{k});
            end;
            % reset val and sts
            val = [];
            sts = false;
            return;
        end;
    end;
end;
% all collected, check subsasgn validity
if sts
    [sts val] = subsasgn_check(item, substruct('.','val','{}',{1}), val);
end;
if ~sts
    warning('matlabbatch:resolve_deps:subsasgn','Item ''%s'': Dependencies resolved, but not suitable for this item.', subsref(item, substruct('.','name')));
    return;
end;
