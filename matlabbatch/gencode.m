function [str, tag, cind, ccnt] = gencode(item, tag, stoptag, tropts)

% function [str, tag, cind, ccnt] = gencode(item, tag, stoptag, tropts)
% Generate code to recreate any MATLAB struct/cell variable. Classes need
% to implement their class specific equivalent of gencode.
%
% Traversal options
% struct with fields
% stopspec - match spec to stop code generation (not used here)
% dflag    - (not used here)
% clvl     - current level in tree - level is increased if fields of
%            structures or cell items are traversed
% mlvl     - maximum level to generate - range 1 (top level only) to
%            Inf (all levels)
% cnt      - item count - not used outside cfg_item objects
% mcnt     - (not evaluated here)
% 
% To generate code for numerical arrays, MATLABs mat2str function is used
% with precision 18. This should be sufficient to recreate most numerical
% arrays without precision loss, however there may be rounding errors.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: gencode.m 1783 2008-06-03 10:48:24Z volkmar $

rev = '$Rev: 1783 $'; %#ok

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
            if isobject(item)
                % Objects need to have their own gencode method implemented
                error('matlabbatch:gencode:object', 'Code generation for objects of class ''%s'' must be implemented as object method.', class(item));
            elseif ~(isnumeric(item) || islogical(item))
                error('matlabbatch:gencode:unknown', 'Code generation for objects of class ''%s'' not implemented.', class(item));
            end
            if issparse(item)
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

