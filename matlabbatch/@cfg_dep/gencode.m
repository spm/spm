function [str, tag, cind, ccnt] = gencode(item, tag, stoptag, tropts)

% function [str, tag, cind, ccnt] = gencode(item, tag, stoptag, tropts)
% Generate code to recreate an cfg_dep object.
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
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: gencode.m 1472 2008-04-23 12:02:44Z volkmar $

rev = '$Rev: 1472 $';

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

str = {};
fn = fieldnames(item);
for k = 1:numel(item)
    % explicitly create each item in array
    str{end+1} = sprintf('%s(%d) = %s;', tag, k, class(item));
    for l = 1:numel(fn)
        switch class(item(k).(fn{l}))
            case 'struct'
                % substruct fields
                % only create if not empty, use gencode_substructcode to
                % produce a one-line code.
                if numel(item(k).(fn{l})) >= 1
                    % force index (1) if there is exactly one entry
                    tag1 = sprintf('%s(%d).%s', tag, k, fn{l});
                    [str1 tag1] = gencode_substructcode(item(k).(fn{l}), tag1);
                    str = {str{:} str1{:}};
                end;
            otherwise
                % other field should not be indexed
                tag1 = sprintf('%s(%d).%s', tag, k, fn{l});
                [str1 tag1 ccnt1 cind1] = gencode(item(k).(fn{l}), tag1, tag1, tropts);
                str = {str{:} str1{:}};
        end;
    end;
end;