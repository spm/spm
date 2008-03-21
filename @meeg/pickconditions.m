function res = pickconditions(obj, label)
% Method for returning indices of trials of a certain trial type
% FORMAT res = pickconditions(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel

c = conditions(obj);

res = strmatch(label, c);
    