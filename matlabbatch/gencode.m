function [str, tag, cind] = gencode(item, tag, tagctx)

% GENCODE  Generate code to recreate any MATLAB struct/cell variable.
% For any MATLAB variable, this function generates a .m file that
% can be run to recreate it. Classes can implement their class specific
% equivalent of gencode with the same calling syntax. By default, classes
% are treated similar to struct variables.
%
% [str, tag, cind] = gencode(item, tag, tagctx)
% Input arguments:
% item - MATLAB variable to generate code for (the variable itself, not its
%        name)
% tag     - optional: name of the variable, i.e. what will be displayed left
%           of the '=' sign. This can also be a valid struct/cell array
%           reference, like 'x(2).y'. If not provided, inputname(1) will be
%           used.
% tagctx  - optional: variable names not to be used (e.g. keywords,
%           reserved variables). A cell array of strings.
% Output arguments:
% str  - cellstr containing code lines to reproduce the input variable
% tag  - name of the generated variable (equal to input tag)
% cind - index into str to the line where the variable assignment is coded
%        (usually 1st line for non-object variables)
%
% See also GENCODE_RVALUE, GENCODE_SUBSTRUCT, GENCODE_SUBSTRUCTCODE.
%
% This code has been developed as part of a batch job configuration
% system for MATLAB. See  
%      http://sourceforge.net/projects/matlabbatch
% for details about the original project.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: gencode.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

if nargin < 2
    tag = inputname(1);
end;
if nargin < 3
    tagctx = {};
end
if isempty(tag)
    tag = genvarname('val', tagctx);
end;
% Item count
cind = 1;

% try to generate rvalue code
[rstr, sts] = gencode_rvalue(item);
if sts
    lvaleq = sprintf('%s = ', tag);
    if numel(rstr) == 1
        str{1} = sprintf('%s%s;', lvaleq, rstr{1});
    else
        str = cell(size(rstr));
        indent = {repmat(' ', 1, numel(lvaleq)+1)};
        str{1} = sprintf('%s%s', lvaleq, rstr{1});
        str(2:end-1) = strcat(indent, rstr(2:end-1));
        str{end} = sprintf('%s%s;', indent{1}, rstr{end});
        if numel(str) > 10
            % add cell mode comment to structure longer output
            str = [{'%%'} str(:)' {'%%'}];
        end
    end
else   
    switch class(item)
        case 'char'
            str = {};
            szitem = size(item);
            subs = gensubs('()', {':',':'}, szitem(3:end));
            for k = 1:numel(subs)
                substag = gencode_substruct(subs{k}, tag);
                str1 = gencode(subsref(item, subs{k}), substag{1}, tagctx);
                str  = [str(:)' str1(:)'];
            end
        case 'cell'
            str = {};
            szitem = size(item);
            subs = gensubs('{}', {}, szitem);
            for k = 1:numel(subs)
                substag = gencode_substruct(subs{k}, tag);
                str1 = gencode(subsref(item, subs{k}), substag{1}, tagctx);
                str  = [str(:)' str1(:)'];
            end
        case 'struct'
            str = gencode_structobj(item, tag, tagctx);
        otherwise
            if isobject(item) || ~(isnumeric(item) || islogical(item))
                % This branch is hit for objects without a gencode method
                try
                    % try to generate code in a struct-like fashion
                    str = gencode_structobj(item, tag, tagctx);
                catch
                    % failed - generate a warning in generated code and
                    % warn directly
                    str = {sprintf('warning(''%s: No code generated for object of class %s.'')', tag, class(item))};
                    if any(exist('cfg_message') == 2:6)
                        cfg_message('matlabbatch:gencode:unknown', ...
                                    '%s: Code generation for objects of class ''%s'' must be implemented as object method.', tag, class(item));
                    else
                        warning('gencode:unknown', ...
                                '%s: Code generation for objects of class ''%s'' must be implemented as object method.', tag, class(item));
                    end
                end
            elseif issparse(item)
                % recreate sparse matrix from indices
                [tmpi, tmpj, tmps] = find(item);
                [stri, tagi, cindi] = gencode(tmpi);
                [strj, tagj, cindj] = gencode(tmpj);
                [strs, tags, cinds] = gencode(tmps);
                str = [stri(:)' strj(:)' strs(:)'];
                cind = cind + cindi + cindj + cinds;
                str{end+1} = sprintf('%s = sparse(tmpi, tmpj, tmps);', tag);
            else
                str = {};
                szitem = size(item);
                subs = gensubs('()', {':',':'}, szitem(3:end));
                for k = 1:numel(subs)
                    substag = gencode_substruct(subs{k}, tag);
                    str1 = gencode(subsref(item, subs{k}), substag{1}, tagctx);
                    str  = [str(:)' str1(:)'];
                end
            end
    end
end

function subs = gensubs(type, initdims, sz)
% generate a cell array of subscripts into trailing dimensions of
% n-dimensional arrays. Type is the subscript type (either '()' or '{}'),
% initdims is a cell array of leading subscripts that will be prepended to
% the generated subscripts and sz contains the size of the remaining
% dimensions.

% deal with special case of row vectors - only add one subscript in this
% case
if numel(sz) == 2 && sz(1) == 1 && isempty(initdims)
    ind = 1:sz(2);
else
    % generate index array, rightmost index varying fastest
    ind = 1:sz(1);
    for k = 2:numel(sz)
        ind = [kron(ind, ones(1,sz(k))); kron(ones(1,size(ind,2)), 1:sz(k))];
    end;
end;

subs = cell(1,size(ind,2));
% for each column of ind, generate a separate subscript structure
for k = 1:size(ind,2)
    cellind = num2cell(ind(:,k));
    subs{k} = substruct(type, [initdims(:)' cellind(:)']);
end;

function str = gencode_structobj(item, tag, tagctx)

% Create code for a struct array. Also used as fallback for object
% arrays, if the object does not provide its own gencode implementation.

citem = class(item);
% try to figure out fields/properties that can be set
if isobject(item) && exist('metaclass','builtin')
    mobj = metaclass(item);
    % Only create code for properties which are
    % * not dependent or dependent and have a SetMethod
    % * not constant
    % * not abstract
    % * have public SetAccess
    sel = cellfun(@(cProp)(~cProp.Constant && ...
        ~cProp.Abstract && ...
        (~cProp.Dependent || ...
        (cProp.Dependent && ...
        ~isempty(cProp.SetMethod))) && ...
        strcmp(cProp.SetAccess,'public')),mobj.Properties);
    fn = cellfun(@(cProp)subsref(cProp,substruct('.','Name')),mobj.Properties(sel),'uniformoutput',false);
else
    % best guess
    fn = fieldnames(item);
end
if isempty(fn)
    if isstruct(item)
        str{1} = sprintf('%s = struct([]);', tag);
    else
        str{1} = sprintf('%s = %s;', tag, citem);
    end
elseif isempty(item)
    if isstruct(item)
        fn = strcat('''', fn, '''', ', {}');
        str{1} = sprintf('%s = struct(', tag);
        for k = 1:numel(fn)-1
            str{1} = sprintf('%s%s, ', str{1}, fn{k});
        end
        str{1} = sprintf('%s%s);', str{1}, fn{end});
    else
        str{1} = sprintf('%s = %s.empty;', tag, citem);
    end        
elseif numel(item) == 1
    if isstruct(item)
        str = {};
    else
        str{1} = sprintf('%s = %s;', tag, citem);
    end
    for l = 1:numel(fn)
        str1 = gencode(item.(fn{l}), sprintf('%s.%s', tag, fn{l}), tagctx);
        str  = [str(:)' str1(:)'];
    end
else
    str = {};
    szitem = size(item);
    subs = gensubs('()', {}, szitem);
    for k = 1:numel(subs)
        if ~isstruct(item)
            str{end+1} = sprintf('%s = %s;', gencode_substruct(subs{k}, tag), citem);
        end
        for l = 1:numel(fn)
            csubs = [subs{k} substruct('.', fn{l})];
            substag = gencode_substruct(csubs, tag);
            str1 = gencode(subsref(item, csubs), substag{1}, tagctx);
            str  = [str(:)' str1(:)'];
        end
    end
end
