function names = fieldnames(this)
% Fieldnames method for GIfTI objects
% FORMAT names = fieldnames(this)
% this   -  GIfTI object
% names  -  field names
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2008-2023 Wellcome Centre for Human Neuroimaging


if numel(this) > 1, warning('Only handle scalar objects yet.'); end

pfn = {'vertices','faces','normals','cdata','mat','labels','indices'};

names = unique(pfn(isintent(this,pfn)));
