function [sts, val] = subsasgn_check(item,subs,val)

% function [sts, val] = subsasgn_check(item,subs,val)
% Perform validity checks for cfg_entry inputs. Does not yet support
% evaluation of inputs.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_check.m 1409 2008-04-15 10:22:46Z volkmar $

rev = '$Rev: 1409 $';

sts = true;
checkstr = sprintf('Item ''%s'', field ''%s''', subsref(item,substruct('.','name')), subs(1).subs);
switch subs(1).subs
    case {'num'}
        % special num treatment - num does describe the dimensions of
        % input in cfg_entry items, not a min/max number
        sts = isnumeric(val) && (isempty(val) || numel(val)>=2 && all(val(:) >= 0));
        if ~sts
             warning('matlabbatch:cfg_entry:subsasgn_check:num', ...
                     '%s: Value must be empty or a vector of non-negative numbers with at least 2 elements', ...
                     checkstr);
        end;
    case {'val'}
        % perform validity checks - subsasgn_check should be called with
        % a cell containing one item
        if ~iscell(val)
            warning('matlabbatch:cfg_entry:subsasgn_check:iscell', '%s: Value must be a cell.', checkstr);
            sts = false;
            return;
        end;
        if isempty(val)
            val = {};
        else
            % check whether val{1} is a valid element
            [sts vtmp] = valcheck(item,val{1});
            val{1} = vtmp;
        end;
    case {'strtype'}
        strtypes = {'s','e','f','n','w','i','r','c','x','p'};
        sts = isempty(val) || (ischar(val) && ...
                               any(strcmp(val, strtypes)));
        if ~sts
            warning('matlabbatch:cfg_entry:subsasgn_check:strtype', ...
                    '%s: Value must be a valid strtype.', checkstr);
        end;
end;

function [sts, val] = valcheck(item,val)
% taken from spm_jobman/stringval
% spm_eeval goes into GUI
sts = true;
checkstr = sprintf('Item ''%s'', field ''val''', subsref(item,substruct('.','name')));
if ~isa(val,'cfg_dep')
    switch item.strtype
        case {'s'}
            if ~ischar(val)
                warning('matlabbatch:cfg_entry:subsasgn_check:strtype', ...
                        '%s: Item must be a string.', checkstr);
                sts = false;
            else
                [sts val] = numcheck(item,val);
            end;
        case {'s+'}
            warning('matlabbatch:cfg_entry:subsasgn_check:strtype', ...
                    '%s: FAILURE: Cant do s+ yet', checkstr);
        case {'f'}
            % test whether val is a function handle or a name of an
            % existing function
            sts = isempty(val) || isa(val, 'function_handle') || ...
                  (ischar(val) && any(exist(val) == 2:6));
            if ~sts
                warning('matlabbatch:cfg_entry:subsasgn_check:strtype', ...
                        '%s: Item must be a function handle or function name.', ...
                        checkstr);
            end;
        case {'n'}
            tol = 4*eps;
            sts = isempty(val) || (isnumeric(val) && all(val(:) >= 1) && ...
                                   all(abs(round(val(:))-val(:)) <= tol));
            if ~sts
                warning('matlabbatch:cfg_entry:subsasgn_check:strtype', ...
                        '%s: Item must be an array of natural numbers.', checkstr);
                return;
            end;
            [sts val] = numcheck(item,val);
        case {'i'}
            tol = 4*eps;
            sts = isempty(val) || (isnumeric(val) && ...
                                   all(abs(round(val(:))-val(:)) <= tol));
            if ~sts
                warning('matlabbatch:cfg_entry:subsasgn_check:strtype', ...
                        '%s: Item must be an array of integers.', checkstr);
                return;
            end;
            [sts val] = numcheck(item,val);
        case {'r'}
            sts = isempty(val) || (isnumeric(val) && all(isreal(val(:))));
            if ~sts
                warning('matlabbatch:cfg_entry:subsasgn_check:strtype', ...
                        '%s: Item must be an array of real numbers.', checkstr);
                return;
            end;
            [sts val] = numcheck(item,val);
        case {'w'}
            tol = 4*eps;
            sts = isempty(val) || (isnumeric(val) && all(val(:) >= 0) && ...
                                   all(abs(round(val(:))-val(:)) <= tol));
            if ~sts
                warning('matlabbatch:cfg_entry:subsasgn_check:strtype', ...
                        '%s: Item must be an array of whole numbers.', checkstr);
                return;
            end;
            [sts val] = numcheck(item,val);
        otherwise
            % only do size check for other strtypes
            [sts val] = numcheck(item,val);
    end;
end;

function [sts, val] = numcheck(item,val)
checkstr = sprintf('Item ''%s'', field ''val''', subsref(item,substruct('.','name')));
% allow arbitrary size, if num field is empty
sts = true;
csz = size(val);
if ~isempty(item.num)
    if item.strtype == 's' && numel(item.num) == 2
        % interpret num field as [min max] # elements
        sts = item.num(1) <= numel(val) && numel(val) <= item.num(2);
        if ~sts
            warning('matlabbatch:cfg_entry:subsasgn_check:notsize', ...
                    '%s: Size mismatch (required [%s], present [%s]).', ...
                    checkstr, num2str(item.num), num2str(csz));
        end;
    else
        ind = item.num>0 & isfinite(item.num);
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
        if any(item.num(ind)-csz(ind))
            if any(item.num(ind)-cszt(ind))
                warning('matlabbatch:cfg_entry:subsasgn_check:notsize', ...
                        '%s: Size mismatch (required [%s], present [%s]).', ...
                        checkstr, num2str(item.num), num2str(csz));
                sts = false;
                return
            else
                val = val';
                warning('matlabbatch:cfg_entry:subsasgn_check:transp', ...
                        '%s: Value transposed to match required size [%s].', ...
                        checkstr, num2str(item.num));
            end;
        end;
    end;
end;