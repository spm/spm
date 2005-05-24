function t = numel(obj)
% Number of simple file arrays involved.
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: numel.m 174 2005-05-24 11:03:32Z john $


% Should be this, but it causes problems when accessing
% obj as a structure.
%t = prod(size(obj));

t  = numel(struct(obj));
