function [sts, val] = subsasgn_check(item,subs,val)

% function [sts, val] = subsasgn_check(item,subs,val)
% Perform assignment checks for .num field. Checks for .val field could
% include filtering and numel checks, if inputs are not passed as
% reference.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_check.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

sts = true;
switch subs(1).subs
    case {'num'}
        sts = subsasgn_check_num(val);
    case {'filter'}
        try
            val = cellstr(val);
        catch
            sts = false;
        end
    case {'val'}
        % val{1} should be a cellstr or a cfg_dep
        sts = iscell(val) && (isempty(val) || isempty(val{1}) || ...
                              isa(val{1}, 'cfg_dep') || iscellstr(val{1}));
        if ~sts
            cfg_message('matlabbatch:checkval', ...
                    '%s: Value must be either empty, a cellstr or a cfg_dep object.', subsasgn_checkstr(item,subs));
        end
        if sts && ~isempty(val) 
            if iscellstr(val{1})
                % do filtering and .num checks
                % this is already done in interactive mode, but not in batch
                % mode (e.g. after resolve_deps).
                if strcmpi(item.filter, 'dir')
                    % don't filter dirs - they would be filtered by
                    % item.ufilter only, which is ignored here
                    val1 = val{1};
                else
                    % don't filter for item.ufilter - this may have been
                    % overridden by user interface
                    [val1, sts1] = cfg_getfile('filter',val{1},item.filter,'.*');
                end
                if numel(val1) < item.num(1)
                    sts = false;
                    cfg_message('matlabbatch:checkval', ...
                                ['%s: Number of matching files (%d) less than ' ...
                                 'required (%d).'], ...
                                subsasgn_checkstr(item,subs), numel(val1), item.num(1));
                elseif numel(val1) > item.num(2)
                    cfg_message('matlabbatch:checkval', ...
                                ['%s: Number of matching files larger than ' ...
                                 'max allowed, keeping %d/%d files.'], ...
                                subsasgn_checkstr(item,subs), item.num(2), numel(val{1}));
                    val1 = val1(1:item.num(2));
                end
                val = {val1};
            elseif isa(val{1}, 'cfg_dep')
                % {val{1}.tgt_spec} does not seem to work properly
                cval = val{1};
                % Check dependency match
                sts2 = cellfun(@(cspec)match(item,cspec),{cval.tgt_spec});
                if ~all(sts2)
                    cfg_message('matlabbatch:checkval', ...
                        '%s: Dependency does not match.', subsasgn_checkstr(item,subs));
                end
                val{1} = val{1}(sts2);
                sts = any(sts2);
            end
        end
end
