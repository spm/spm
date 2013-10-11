function [str, tag, cind, ccnt] = gencode_item(item, tag, tagctx, stoptag, tropts)

% function [str, tag, cind, ccnt] = gencode_item(item, tag, tagctx, stoptag, tropts)
% Generate code to recreate a cfg_(m)choice item. This code does not deal with
% arrays of cfg_items, such a configuration should not exist with the
% current definition of a configuration tree.
%
% Traversal options
% struct with fields
% stopspec - match spec to stop forced setting of eflag
% dflag    - if set to true, don't create code for .val children (code
%            for .val field is created)
% clvl     - current level in tree
% mlvl     - maximum level to force settings - range 1 (top level only) to
%            Inf (all levels)
% cnt      - item count - used for unique tags
% mcnt     - (not evaluated here)
%
% This function is identical for cfg_choice and cfg_mchoice classes.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: gencode_item.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

%% Parent object
% Generate generic object
itropts = tropts;
istoptag = stoptag;
if tropts.dflag
    % don't descend into .val tree of cfg_item
    itropts.clvl = 1;
    itropts.mlvl = 1;
    istoptag     = '';
end;
[str, tag, cind, ccnt] = gencode_item(item.cfg_item, tag, tagctx, istoptag, itropts);
tagctx = [tagctx {tag}];
% Check whether to generate code - ccnt == 0 means that generic object did
% not return code
if (tropts.clvl > tropts.mlvl || (~isempty(tropts.stopspec) && match(item, tropts.stopspec))) || ccnt == 0
    str = {};
    cind = [];
    ccnt = 0;
    return;
end;
% Reclassify generic object
str{cind} = sprintf('%s         = %s;', tag, class(item));
%% Values
% Generate values field
if numel(item.values) > 0
    % Traverse values{:} tree
    cstr = {};
    % Update clvl
    ctropts = tropts;
    ctropts.clvl = ctropts.clvl + 1;
    ctropts.cnt  = ctropts.cnt + ccnt;
    ctag = cell(size(item.values));
    for k = 1:numel(item.values)
        % tags are used as variable names and need to be unique in the
        % context of this .values tag. This includes the item's tag itself
        % and the tags of its immediate children.
        ctag{k} = genvarname(subsref(item.values{k}, substruct('.','tag')), ...
                             tagctx);
        [ccstr, ctag{k}, ccind, cccnt] = gencode_item(item.values{k}, ctag{k}, ...
                                              tagctx, stoptag, ctropts);
        tagctx = [tagctx ctag(k)];
        if ~isempty(ccstr)
            % Child has returned code
            cstr = [cstr(:)' ccstr(:)'];
            ctropts.cnt = ctropts.cnt + cccnt;
            ccnt = ccnt + cccnt;
        end;
    end;
    % Update position of class definition
    cind = cind+numel(cstr);
    % Prepend code of children
    str = [cstr(:)' str(:)'];
    str{end+1} = sprintf('%s.values  = {%s};', tag, sprintf('%s ', ctag{:}));
end;
