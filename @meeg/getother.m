function res = getother(obj, name)
% Method for retrieving stuff from the field 'other'
% FORMAT res = getother(obj, name)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel

if isfield(obj.other, name)
    eval(['res = obj.other.' name ';']);
else
    res = [];
end
