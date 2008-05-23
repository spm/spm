function [str, tag, cind, ccnt] = gencode(item, tag, tagctx, stoptag, tropts)

% function [str, tag, cind, ccnt] = gencode(item, tag, tagctx, stoptag, tropts)
% Generate code to recreate a generic item. This code does not deal with
% arrays of cfg_items, such a configuration should not exist with the
% current definition of a configuration tree.
%
% Traversal options
% struct with fields
% stopspec - match spec to stop code generation
% dflag    - (not used here)
% clvl     - current level in tree
% mlvl     - maximum level to generate - range 1 (top level only) to
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
% $Id: gencode.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

%% Parent object
% Generate generic object
[str tag cind ccnt] = gencode(item.cfg_item, tag, tagctx, stoptag, tropts);
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
%% Labels
% Generate labels field
str1 = gencode(item.labels, sprintf('%s.labels', tag), stoptag, tropts);
str = {str{:} str1{:}};
%% Values
% Generate values field
str1 = gencode(item.values, sprintf('%s.values', tag), stoptag, tropts);
str = {str{:} str1{:}};
%% Def
% Do not create deprecated def field