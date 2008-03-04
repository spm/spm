function [str tag cind ccnt] = gencode(item, tag, stoptag, tropts)

% Generate code to recreate a cfg_choice item. This code does not deal with
% arrays of cfg_items, such a configuration should not exist with the
% current definition of a configuration tree.
%
% Traversal options
% struct with fields
% stopspec - match spec to stop forced setting of eflag
% dflag    - (not used here)
% clvl     - current level in tree
% mlvl     - maximum level to force settings - range 1 (top level only) to
%            Inf (all levels)
% cnt      - item count - used for unique tags
% mcnt     - (not evaluated here)
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: gencode.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

%% Parent object
% Generate generic object
[str tag cind ccnt] = gencode(item.cfg_item, tag, stoptag, tropts);
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
    % Traverse values{:} tree, if items are cfg_items
    cstr = {};
    % Update clvl
    tropts.clvl = tropts.clvl + 1;
    tropts.cnt  = tropts.cnt + ccnt;
    ctag = {};
    for k = 1:numel(item.values)
        % tags are used as variable names and need to be unique in the
        % context of this .values tag. This includes the item's tag itself
        % and the tags of its immediate children.
        ctag{k} = genvarname(subsref(item.values{k}, substruct('.','tag')), ...
                             {ctag{:} tag});
        [ccstr ctag{k} ccind cccnt] = gencode(item.values{k}, ctag{k}, stoptag, tropts);
        if ~isempty(ccstr)
            % Child has returned code
            cstr = {cstr{:} ccstr{:}};
            tropts.cnt = tropts.cnt + cccnt;
            ccnt = ccnt + cccnt;
        end;
    end;
    % Update position of class definition
    cind = cind+numel(cstr);
    % Prepend code of children
    str = {cstr{:} str{:}};
    str{end+1} = sprintf('%s.values  = {%s};', tag, sprintf('%s ', ctag{:}));
end;
