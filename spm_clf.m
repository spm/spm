function spm_clf
% Selective deletion of objects in the current figure window
% FORMAT spm_clf
%___________________________________________________________________________
%
% deletes the children of the current figure other that objects
% with 'Tag' attribute = 'NoDelete'.   This serves as a partial
% clf, preventing some objects (e.g. menu bar objects) from being
% deleted
%
%__________________________________________________________________________
% %W% %E%


%---------------------------------------------------------------------------
h = get(gcf,'Children');
for i = 1:length(h)
	if ~strcmp(get(h(i),'Tag'),'NoDelete'); delete(h(i)); end
end
