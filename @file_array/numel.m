function t = numel(obj)
% Number of simple file arrays involved.
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

%
% $Id: numel.m 7147 2017-08-03 14:07:01Z spm $


% Should be this, but it causes problems when accessing
% obj as a structure.
%t = prod(size(obj));

t  = numel(struct(obj));
