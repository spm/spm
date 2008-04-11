function [str, tag, cind, ccnt] = gencode(item, tag, tagctx, stoptag, tropts)

% function [str, tag, cind, ccnt] = gencode(item, tag, tagctx, stoptag, tropts)
% Generate code to recreate a cfg_exbranch item. This code first generates
% code for the parent cfg_branch item and adds code for its own fields.
% Note that function references will be broken if they refer to a local
% function in the original config file. This code does not deal with arrays
% of cfg_items, such a configuration should not exist with the current
% definition of a configuration tree.
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
% $Id: gencode.m 1366 2008-04-11 10:24:17Z volkmar $

rev = '$Rev: 1366 $';

%% Parent object
% Generate branch object
[str tag cind ccnt] = gencode(item.cfg_branch, tag, tagctx, stoptag, tropts);
% Check whether to generate code - ccnt == 0 means that generic object did
% not return code
if (tropts.clvl > tropts.mlvl || (~isempty(tropts.stopspec) && match(item, tropts.stopspec))) || ccnt == 0
    str = {};
    cind = [];
    ccnt = 0;
    return;
end;
% Reclassify branch object
str{cind} = sprintf('%s         = %s;', tag, class(item));
%% Generate code for other fields
funs = {'prog', 'vfiles', 'vout'};
for k = 1:numel(funs)
    if ~isempty(item.(funs{k}))
        str1 = gencode(item.(funs{k}), sprintf('%s.%s', tag, funs{k}), stoptag, tropts);
        str = {str{:} str1{:}};
    end;
end;
%% Modality
% Generate modality field
if numel(item.modality) > 0
    str1 = gencode(item.modality, sprintf('%s.modality', tag), stoptag, tropts);
    str = {str{:} str1{:}};
end;
