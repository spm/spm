function [str, tag, cind, ccnt] = gencode(item, tag, stoptag, tropts)

% function [str, tag, cind, ccnt] = gencode(item, tag, stoptag, tropts)
% Generate code to recreate any MATLAB struct/cell variable. Classes need
% to implement their class specific equivalent of gencode with the same
% calling syntax.
%
% Input arguments:
% item - MATLAB variable to generate code for (the variable itself, not its
%        name)
% tag  - optional: name of the variable, i.e. what will be displayed left
%        of the '=' sign. This can also be a valid struct/cell array
%        reference, like 'x(2).y'. If not provided, inputname(1) will be
%        used.
% stoptag - optional: tag which is used if struct/cell traversal stopped
%        due to tropts limitations. If not provided, tag will be used.
% tropts  - optional: traversal options - struct with fields
% .stopspec - (only used for matlabbatch objects)
% .dflag    - (only used for matlabbatch objects)
% .clvl     - current level in variable - level is increased if fields of
%             structures or cell items are traversed
% .mlvl     - maximum level to generate code for - range 1 (top level only)
%             to Inf (all levels)
% .cnt      - (only used for matlabbatch objects)
% .mcnt     - (only used for matlabbatch objects)
%
% Output arguments:
% str  - cellstr containing code lines to reproduce 
% tag  - name of the generated variable (equal to input tag)
% cind - index into str to the line where the variable assignment is coded
%        (usually 1st line for non-object variables)
% ccnt - item count (outside matlabbatch objects always 1)
%
% For scalar, 1D or 2D char, numeric or cell arrays whose contents can be
% written as a MATLAB array, the helper function GENCODE_RVALUE will be
% called. This function can also be used on its own. It will produce code
% to generate the array, but without a left hand side assignment.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: gencode.m 2657 2009-01-27 16:24:01Z volkmar $

rev = '$Rev: 2657 $'; %#ok

if nargin < 2
    tag = inputname(1);
end;
if nargin < 3
    stoptag = tag;
end;
if nargin < 4
    tropts = cfg_tropts({{}},1,inf,1,inf,true);
end;

%% Get variable name
% Check whether to generate code
if tropts.clvl > tropts.mlvl
    % Stopping - tag based on stoptag, tag of item and expected new item count
    tag = genvarname(sprintf('%s_%s_0001', stoptag, tag));
    str = {};
    cind = [];
    ccnt = 0;
    return;
else
    % Tag based on item count
    if isempty(tag)
        tag = genvarname(sprintf('val_%04d', tropts.cnt));
    end;
end;
% Item count
ccnt = 1;
cind = 1;

% try to generate rvalue code
[rstr sts] = gencode_rvalue(item);
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
                [str1 tag1 cind1 ccnt1] = gencode(subsref(item, subs{k}), substag{1}, stoptag, tropts);
                str = {str{:} str1{:}};
            end
        case 'cell'
            str = {};
            szitem = size(item);
            subs = gensubs('{}', {}, szitem);
            tropts.clvl = tropts.clvl + 1;
            for k = 1:numel(subs)
                substag = gencode_substruct(subs{k}, tag);
                [str1 tag1 cind1 ccnt1] = gencode(subsref(item, subs{k}), substag{1}, stoptag, tropts);
                str = {str{:} str1{:}};
            end
        case 'struct'
            fn = fieldnames(item);
            if isempty(fn)
                str{1} = sprintf('%s = struct([]);', tag);
            elseif isempty(item)
                fn = strcat('''', fn, '''', ', {}');
                str{1} = sprintf('%s = struct(', tag);
                for k = 1:numel(fn)-1
                    str{1} = sprintf('%s%s, ', str{1}, fn{k});
                end
                str{1} = sprintf('%s%s);', str{1}, fn{end});
            elseif numel(item) == 1
                str = {};
                tropts.clvl = tropts.clvl + 1;
                for l = 1:numel(fn)
                    [str1 tag1 cind1 ccnt1] = gencode(item.(fn{l}), sprintf('%s.%s', tag, fn{l}), stoptag, tropts);
                    str = {str{:} str1{:}};
                end
            else
                str = {};
                szitem = size(item);
                subs = gensubs('()', {}, szitem);
                tropts.clvl = tropts.clvl + 1;
                for k = 1:numel(subs)
                    for l = 1:numel(fn)
                        csubs = [subs{k} substruct('.', fn{l})];
                        substag = gencode_substruct(csubs, tag);
                        [str1 tag1 cind1 ccnt1] = gencode(subsref(item, csubs), substag{1}, stoptag, tropts);
                        str = {str{:} str1{:}};
                    end
                end
            end
        otherwise
            if isobject(item) || ~(isnumeric(item) || islogical(item))
                str = {sprintf('warning(''%s: No code generated for object of class %s.'')', tag, class(item))};
                % Objects need to have their own gencode method implemented
                cfg_message('matlabbatch:gencode:unknown', ...
                            '%s: Code generation for objects of class ''%s'' must be implemented as object method.', tag, class(item));
            elseif issparse(item)
                % recreate sparse matrix from indices
                [tmpi tmpj tmps] = find(item);
                [stri tagi cindi ccnti] = gencode(tmpi);
                [strj tagj cindj ccntj] = gencode(tmpj);
                [strs tags cinds ccnts] = gencode(tmps);
                str = {stri{:} strj{:} strs{:}};
                cind = cind + cindi + cindj + cinds;
                str{end+1} = sprintf('%s = sparse(tmpi, tmpj, tmps);', tag);
            else
                str = {};
                szitem = size(item);
                subs = gensubs('()', {':',':'}, szitem(3:end));
                for k = 1:numel(subs)
                    substag = gencode_substruct(subs{k}, tag);
                    [str1 tag1 cind1 ccnt1] = gencode(subsref(item, subs{k}), substag{1}, stoptag, tropts);
                    str = {str{:} str1{:}};
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

% for each column of ind, generate a separate subscript structure
for k = 1:size(ind,2)
    cellind = num2cell(ind(:,k));
    subs{k} = substruct(type, {initdims{:} cellind{:}});
end;

