function res = pickconditions(obj, label)
% Method for returning indices of trials of a certain trial type
% FORMAT res = pickconditions(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
%$Id: conditions.m 1125 2008-01-30 12:12:18Z vladimir $

c = conditions(obj);

res = strmatch(label, c);

if isempty(res)
    res = [];
end
    