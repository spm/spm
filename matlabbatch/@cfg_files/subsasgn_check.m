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
% $Id: subsasgn_check.m 1862 2008-06-30 14:12:49Z volkmar $

rev = '$Rev: 1862 $'; %#ok

sts = true;
switch subs(1).subs
    case {'num'}
        sts = subsasgn_check_num(val);
    case {'val'}
        % val{1} should be a cellstr or a cfg_dep
        sts = iscell(val) && (isempty(val) || isempty(val{1}) || ...
                              isa(val{1}, 'cfg_dep') || iscellstr(val{1}));
        if ~sts
            cfg_message('matlabbatch:checkval', ...
                    '%s: Value must be either empty, a cellstr or a cfg_dep object.', subsasgn_checkstr(item,subs));
        end
        if ~isempty(val) && iscellstr(val{1})
            % do filtering and .num checks
            % this is already done in interactive mode, but not in batch
            % mode (e.g. after resolve_deps).
            if strcmpi(item.filter,'image') || strcmpi(item.filter,'nifti')
                typ = 'extimage';
            else
                typ = item.filter;
            end;
            % don't filter for item.ufilter - this may have been
            % overridden by user interface
            [val{1} sts1] = cfg_getfile('filter',val{1},typ,'.*',Inf);
            if numel(val{1}) < item.num(1)
                sts = false;
                cfg_message('matlabbatch:checkval', ...
                        ['%s: Number of matching files (%d) less than ' ...
                         'required (%d).'], ...
                        subsasgn_checkstr(item,subs), numel(val{1}), item.num(1));
            elseif numel(val{1}) > item.num(2)
                cfg_message('matlabbatch:checkval', ...
                        ['%s: Number of matching files larger than ' ...
                         'max allowed, keeping %d/%d files.'], ...
                        subsasgn_checkstr(item,subs), item.num(2), numel(val{1}));
                val{1} = val{1}(1:item.num(2));
            end
        end
end
