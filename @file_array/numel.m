function t = numel(obj)
% Number of simple file arrays involved.
% _______________________________________________________________________
% %W% John Ashburner %E%

% Should be this, but it causes problems when accessing
% obj as a structure.
%t = prod(size(obj));

t  = numel(struct(obj));
