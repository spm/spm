function spm_clf(F)
% Selective deletion of objects in the current figure window
% FORMAT spm_clf(F)
% F	- Figure number, or 'Tag' string of figure(s) to clear
%___________________________________________________________________________
%
% Deletes the children of the specified figure, other than objects with
% 'Tag' attribute = 'NoDelete'. This serves as a partial clf,
% preventing some objects (e.g. menu bar objects) from being deleted.
%
% If the current window is 'Tag'ged interactive, then the figures name
% is cleared and the pointer reset.
%
% F Defaults to the current figure, if there is one.
%
% This is just a gateway to spm_figure('Clear',F)
%__________________________________________________________________________
% %W% %E%

%-Call spm_figure
%---------------------------------------------------------------------------
if nargin==0
	spm_figure('Clear')
else
	spm_figure('Clear',F)
end
