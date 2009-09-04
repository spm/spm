function [str, tag, cind, ccnt] = gencode_item(item, tag, tagctx, stoptag, tropts)

% function [str, tag, cind, ccnt] = gencode_item(item, tag, tagctx, stoptag, tropts)
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
% $Id: gencode_item.m 3355 2009-09-04 09:37:35Z volkmar $

rev = '$Rev: 3355 $'; %#ok

%% Parent object
% if there are function handles in .prog, .vout, .vfiles add their names to
% tagctx
if ~isempty(item.prog) && isa(item.prog, 'function_handle')
    functx = {func2str(item.prog)};
else
    functx = {};
end
if ~isempty(item.vout) && isa(item.vout, 'function_handle')
    functx{end+1} = func2str(item.vout);
end
if ~isempty(item.vfiles) && isa(item.vfiles, 'function_handle')
    functx{end+1} = func2str(item.vfiles);
end
tagctx = [tagctx functx];
% Generate branch object
[str tag cind ccnt] = gencode_item(item.cfg_branch, tag, tagctx, stoptag, tropts);
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
        str1 = gencode(item.(funs{k}), sprintf('%s.%s', tag, funs{k}), ...
                       tagctx);
        str = [str(:)' str1(:)'];
    end;
end;
%% Modality
% Generate modality field
if numel(item.modality) > 0
    str1 = gencode(item.modality, sprintf('%s.modality', tag), tagctx);
    str = [str(:)' str1(:)'];
end;
