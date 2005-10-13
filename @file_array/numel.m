function t = numel(obj)
% Number of simple file arrays involved.
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: numel.m 253 2005-10-13 15:31:34Z guillaume $


% Should be this, but it causes problems when accessing
% obj as a structure.
%t = prod(size(obj));

t  = numel(struct(obj));
