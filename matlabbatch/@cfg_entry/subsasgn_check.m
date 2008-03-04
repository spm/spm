function [sts val] = subsasgn_check(item,subs,val)

% function [sts val] = subsasgn_check(item,subs,val)
% Perform validity checks for cfg_entry inputs. Does not yet support
% evaluation of inputs.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_check.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

sts = true;
checkstr = sprintf('Item ''%s'', field ''%s''', subsref(item,substruct('.','name')), subs(1).subs);
switch subs(1).subs
    case {'num'}
        sts = subsasgn_check_num(val);
    case {'val'}
        % perform validity checks
        if numel(subs) == 1
            % pass a cell of items
            if ~iscell(val)
                warning('matlabbatch:cfg_entry:subsasgn_check:iscell', '%s: Value must be a cell.', checkstr);
                sts = false;
                return;
            end;
            if numel(val) > 1
                warning('matlabbatch:cfg_entry:subsasgn_check:iscell', '%s: Value must be a single cell.', checkstr);
                sts = false;
                return;
            end;
            if isempty(val)
                val = {};
            else
                % check whether val{1} is a valid element
                sts = valcheck(item,val{1});
            end;
        elseif numel(subs) == 2
            % pass either a cell of items or a single item
            if numel(subs(2).subs) > 1 || subs(2).subs{1} > 1
                warning('matlabbatch:cfg_entry:subsasgn_check:iscell', '%s: Value must be a single item assigned to val{1}.', checkstr);
                sts = false;
                return;
            end;
            switch(subs(2).type)
                case {'()'}
                    sts = valcheck(item,val{1});
                case {'{}'}
                    sts = valcheck(item,val);
            end;
        end;
end;

function sts = valcheck(item,val)
% taken from spm_jobman/stringval
% spm_eeval goes into GUI
sts = true;
checkstr = sprintf('Item ''%s'', field ''val''', subsref(item,substruct('.','name')));
if ~isa(val,'cfg_dep')
    switch item.strtype
        case {'s'}
            if ~ischar(val)
                warning('matlabbatch:cfg_entry:subsasgn_check:strtype', '%s: Item must be a string.', checkstr);
                sts = false;
            end;
        case {'s+'}
            warning('matlabbatch:cfg_entry:subsasgn_check:strtype', '%s: FAILURE: Cant do s+ yet', checkstr);

        otherwise
            ind = item.num>0 & isfinite(item.num);
            csz = size(val);
            if numel(csz) == 2
                % also try transpose for 2D arrays
                cszt = size(val');
            else
                cszt = csz;
            end;
            if numel(item.num) ~= numel(csz)
                warning('matlabbatch:cfg_entry:subsasgn_check:notdim', '%s: Dimension mismatch (required %d, present %d).', checkstr, numel(item.num), numel(csz));
                sts = false;
                return;
            end;
            if any(item.num(ind)-csz(ind)) && any(item.num(ind)-cszt(ind))
                warning('matlabbatch:cfg_entry:subsasgn_check:notsize', '%s: Size mismatch (required [%s], present [%s]).', checkstr, num2str(item.num), num2str(csz))
                sts = false;
                return
            end;
    end;
end;