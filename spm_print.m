function spm_print
% SPM print function that executes PRINTSTR
% FORMAT spm_print
%___________________________________________________________________________
%
% spm_print creates a footnote with details of the current
% session and evaluates the global string variable PRINTSTR
%
%__________________________________________________________________________
% %W% %E%

% create footnote
%---------------------------------------------------------------------------
global PRINTSTR

NAME	        = getenv('USER');			% username
GRAPHSTR  	= ['SPM analysis - date: ' date '  user: ' NAME];
GRAPHSTR  	= GRAPHSTR(GRAPHSTR ~= 10);		% footnote

figure(3)

% create axis and text objects
%---------------------------------------------------------------------------
axes('position',[0.02 0.02 1 1],'Visible','off')
text(0,0,GRAPHSTR,'FontSize',10);

% print command
%---------------------------------------------------------------------------
eval(PRINTSTR)
