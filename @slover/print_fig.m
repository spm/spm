function print_fig(obj, filename, printstr)
% Print slice overlay figure
% FORMAT print_fig(obj, filename, printstr)
%
% Input
% obj       - object
% filename  - optional filename to print to (obj.filename)
% printstr  - optional string giving print command (obj.printstr)
%
% Based on spm_figure print, and including fix from thence for ps
% printing
%__________________________________________________________________________

% Matthew Brett
% $Id: print_fig.m 6623 2015-12-03 18:38:08Z guillaume $

if nargin < 2
    filename = [];
end
if isempty(filename)
    filename = obj.printfile;
end
if nargin < 3
    printstr = '';
end
if isempty(printstr)
    printstr = obj.printstr;
end

%-Note current figure, & switch to figure to print
cF = get(0,'CurrentFigure');
set(0,'CurrentFigure',obj.figure)

%-Temporarily change all units to normalized prior to printing
% (Fixes bizzarre problem with stuff jumping around!)
%-----------------------------------------------------------------------
H  = findobj(get(obj.figure,'Children'),'flat','Type','axes');
un = cellstr(get(H,'Units'));
set(H,'Units','normalized')

%-Print
%-----------------------------------------------------------------------
try
    eval([printstr ' ' filename])
catch
    errstr = lasterr;
    tmp = [find(abs(errstr)==10),length(errstr)+1];
    str = {errstr(1:tmp(1)-1)};
    for i = 1:length(tmp)-1
        if tmp(i)+1 < tmp(i+1)
            str = [str, {errstr(tmp(i)+1:tmp(i+1)-1)}];
        end
    end
    str = [str,  '','- print command is:',['    ',printstr ' ' filename],...
        '','- current directory is:',['    ',pwd],...
        '','            * nothing has been printed *'];
    for i=1:length(str)
        disp(str{i});end
end

set(H,{'Units'},un)
set(0,'CurrentFigure',cF)
