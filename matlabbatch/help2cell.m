function cellhelp = help2cell(topic)
% HELP2CELL - translate help texts into cell arrays
% cellhelp = help2cell(topic)
% Create a cell array of help strings from the MATLAB help on 'topic'.
% If a line ends with a ' ', it is assumed to be continued and the next 
% line will be appended, thus creating one cell per paragraph.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: help2cell.m 3489 2009-10-20 09:49:40Z volkmar $

rev = '$Rev: 3489 $'; %#ok

try
    tmphelp = textscan(help(topic),'%s', 'delimiter',char(10));
    cellhelp = tmphelp{1}(1);
    ch = 1;
    for k = 2:numel(tmphelp{1})
        if ~isempty(cellhelp{ch}) && (cellhelp{ch}(end) == char(32))
            cellhelp{ch} = [cellhelp{ch} tmphelp{1}{k}];
        else
            ch = ch + 1;
            cellhelp{ch} = tmphelp{1}{k};
        end;
    end;
catch
    cellhelp = {sprintf('No help available for topic ''%s''.', ...
                        topic)};
end;
