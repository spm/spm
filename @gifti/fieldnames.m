function names = fieldnames(this)
% Fieldnames method for GIfTI objects
% FORMAT names = fieldnames(this)
% this   -  GIfTI object
% names  -  field names
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: fieldnames.m 2076 2008-09-10 12:34:08Z guillaume $

if numel(this) > 1, warning('Only handle scalar objects yet.'); end

pfn = {'vertices' 'faces' 'normals' 'cdata' 'mat'};

names = pfn(isintent(this,pfn));