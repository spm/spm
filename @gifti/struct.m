function s = struct(this)
% Struct method for GIfTI objects
% FORMAT s = struct(this)
% this   -  GIfTI object
% s      -  a structure containing public fields of the object
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: struct.m 2076 2008-09-10 12:34:08Z guillaume $

names = fieldnames(this);
values = cell(length(names), length(this(:)));

for i=1:length(names)
    [values{i,:}] = subsref(this(:), substruct('.',names{i}));
end
s  = reshape(cell2struct(values,names,1),size(this));