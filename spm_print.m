function spm_print(F)
% SPM print function that executes PRINTSTR
% FORMAT spm_print(F)
% F	- [Optional] Figure to print
%	- defaults to the figure 'Tag'ged as 'Graphics'
%	- If no such figure is found, then CurrentFigure is used, if avaliable
%___________________________________________________________________________
%
% spm_print creates a footnote with details of the current
% session and evaluates the global string variable PRINTSTR
%
%__________________________________________________________________________
% %W% %E%


%-Retrieve print command
%---------------------------------------------------------------------------
global PRINTSTR
if isempty(PRINTSTR)
	PrintCmd = 'print -dps2 fig.ps';
else
	PrintCmd = PRINTSTR;
end

%-Find window to print
%---------------------------------------------------------------------------
if (nargin<1), F = findobj(get(0,'Children'),'Tag','Graphics'); end
if (isempty(F) & length(get(0,'Children'))), F = gcf; end
figure(F)

%-Create footnote
%---------------------------------------------------------------------------
NAME	        = getenv('USER');			% username
GRAPHSTR  	= ['SPM analysis - date: ' date '  user: ' NAME];
GRAPHSTR  	= GRAPHSTR(GRAPHSTR ~= 10);		% footnote

%-Create axis and text objects
%---------------------------------------------------------------------------
axes('position',[0.02 0.02 1 1],'Visible','off')
text(0,0,GRAPHSTR,'FontSize',10);

%-Print command
%---------------------------------------------------------------------------
eval(PrintCmd)
