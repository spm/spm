function cfg = cfg_txtdesc2cfg(fname)

% Create a cfg_item tree from a short-hand text description
% cfg = cfg_txtdesc2cfg(fname)
% This utility reads a text file from fname and creates a configuration
% object tree according to the following grammar rules.
% 
% Each line in the file has the form
%
%  TAGNAME = RIGHTHANDSIDE
%
% where TAGNAME is a valid tag name for a cfg_item object. For each line, a
% cfg_item object with tag TAGNAME will be created. Its class is determined
% by the format of RIGHTHANDSIDE. RIGHTHANDSIDE can be one of
% 
%  (TAGNAME_1 TAGNAME_2 ... TAGNAME_N) - cfg_branch
%
%  {TAGNAME_1 TAGNAME_2 ... TAGNAME_N} - cfg_repeat
%
%  [TAGNAME_1 TAGNAME_2 ... TAGNAME_N] - cfg_mchoice
%
%  |TAGNAME_1 TAGNAME_2 ... TAGNAME_N| - cfg_choice
%
% with .val (for cfg_branch) or .values (all other cases) fields set to the
% cfg_item objects referenced by TAGNAME_1 ... TAGNAME_N.
%
%  TAGNAME_1 - same type as TAGNAME_1, but with tag TAGNAME instead of
%              TAGNAME_1
%
%  'some_matlab_code' - this MATLAB code will be evaluated. Its return
%                       value should be a cfg_* object. The tag of this
%                       object will be set to 'TAGNAME'
%
% The cfg_item object returned will be the one defined in the first line of
% the file. The depth of the substitutions is not limited, but all
% substitutions must finally be resolvable to 'some_matlab_code'.
%
% A valid description would be
%
%  toplevel = {mychoice mybranch}
%  mychoice = |myrepeat myconst|
%  myrepeat = {mymenu}
%  mybranch = (mymchoice myentry)
%  mymchoice = [myfiles myfiles1]
%  myfiles1 = myfiles
%  myfiles  = 'cfg_files'
%  myconst  = 'cfg_const'
%  myentry  = 'cfg_entry'
%  mymenu   = 'cfg_menu'
%
% The resulting object tree will need further adjustment, but it can serve
% as a good starting point for modifying program code. The sequence
%
%  cfg = cfg_txtdesc2cfg('mygrammar.txt');
%  cfgstr = gencode(cfg);
%  cfgchr = sprintf('%s\n',cfgstr{:});
%  clipboard('copy', cfgchr)
%
% will create a cfg_item tree as defined in 'mygrammar.txt', convert it
% into MATLAB code lines, print them into a newline-separated string and
% copy this string into the clipboard. From there, it can be pasted into
% any application (MATLAB editor, external editor ...) where it can be
% processed further.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_txtdesc2cfg.m 4867 2012-08-30 13:04:51Z volkmar $

rev = '$Rev: 4867 $'; 


% Read file
fid = fopen(fname, 'r');
s = textscan(fid,'%s','Delimiter','\n');

% Parse lines into left/right values, separated by equal signs and
% surronding space characters
pat = '^\s*(?<left>\w+)\s*=\s*(?<right>.+)$';
lr  = regexp(s{1}, pat, 'names');
fail = cellfun(@isempty, lr);
if any(fail)
    fprintf('The following lines could not be parsed into left/right pairs:\n')
    ifail = find(fail);
    for k = 1:numel(ifail)
        fprintf('%0*d: %s\n', floor(log10(max(ifail)))+1, ifail(k), s{1}{k});
    end
end
lr = [lr{:}]';

% Collect left values and their cfg_* expressions
lvals = [{lr.left}' cell(numel(lr),1)];

% First pass - resolve all terminal items, enclosed in single quotes and
% evaluate terminal code 
% NB: this can be the result of any MATLAB function call. This can be used
% to import entire configuration trees into a generated configuration.
rvals  = regexp({lr.right}','''([A-Za-z]\w*)''','tokens');
isleaf = ~cellfun(@isempty, rvals);
cfg_leafs = cellfun(@(rv)eval(rv{1}{1}), rvals(isleaf), 'UniformOutput',false);
lvals(isleaf, 2) = cellfun(@(cl,tag)subsasgn(cl, substruct('.','tag'),tag), cfg_leafs, lvals(isleaf,1), 'UniformOutput', false); 

type = {'cfg_branch', 'cfg_repeat', 'cfg_mchoice', 'cfg_choice','none'};

% Keep track on unresolved items
toresolve = cellfun(@isempty,lvals(:,2));
while any(toresolve)
    trind = find(toresolve);
    for k = 1:numel(trind)
        % get child tags
        ctags = regexp(lr(trind(k)).right,'([A-Za-z]\w+)','tokens');
        ctags = [ctags{:}]';
        resolved = true;
        vals = cell(size(ctags));
        for l = 1:numel(ctags)
            clval = strcmp(lvals(:,1), ctags{l});
            resolved = resolved && any(clval) && ~isempty(lvals{clval,2});
            if ~resolved
                break
            end
            vals{l} = lvals{clval,2};
        end
        if resolved
            % determine type of cfg_item
            tsel = ~cellfun(@isempty,regexp(lr(trind(k)).right, {'\(.*\)','{.*}','\[.*\]','\|.*\|'}));
            tsel(end+1)   = ~any(tsel);
            switch type{tsel}
                case 'none'
                    % variable substitution - no enclosure, single item
                    if numel(vals) == 1
                        lvals{trind(k), 2}     = vals{1};
                        lvals{trind(k), 2}.tag = lvals{trind(k), 1};
                    else
                    end
                case 'cfg_branch'
                    lvals{trind(k), 2}     = cfg_branch;
                    lvals{trind(k), 2}.tag = lvals{trind(k), 1};
                    lvals{trind(k), 2}.val = vals;
                case 'cfg_choice'
                    lvals{trind(k), 2}        = cfg_choice;
                    lvals{trind(k), 2}.tag    = lvals{trind(k), 1};
                    lvals{trind(k), 2}.values = vals;
                case 'cfg_mchoice'
                    lvals{trind(k), 2}        = cfg_mchoice;
                    lvals{trind(k), 2}.tag    = lvals{trind(k), 1};
                    lvals{trind(k), 2}.values = vals;
                case 'cfg_repeat'
                    lvals{trind(k), 2}        = cfg_repeat;
                    lvals{trind(k), 2}.tag    = lvals{trind(k), 1};
                    lvals{trind(k), 2}.values = vals;
                    if numel(vals) == 1
                        % Override child tag for single-item repeats
                        lvals{trind(k), 2}.values{1}.tag = lvals{trind(k), 1};
                    end
            end
        end
    end
    toresolvenew = cellfun(@isempty,lvals(:,2));
    if isequal(toresolvenew , toresolve)
        break;
    else
        toresolve = cellfun(@isempty,lvals(:,2));
    end
end
if any(toresolve)
    fprintf('The following lines could not be parsed:\n')
    ifail = find(toresolve);
    for k = 1:numel(ifail)
        fprintf('%0*d: %s\n', floor(log10(max(ifail)))+1, ifail(k), s{1}{k});
    end
    cfg = '';
else
    cfg = lvals{1,2};
end