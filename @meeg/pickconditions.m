function res = pickconditions(obj, label)
% Method for returning indices of trials of a certain trial type
% FORMAT res = pickconditions(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: pickconditions.m 1373 2008-04-11 14:24:03Z spm $

c = conditions(obj);

res = strmatch(label, c);
    
