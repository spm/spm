function res = pickconditions(this, label)
% Method for returning indices of trials of a certain trial type
% FORMAT res = pickconditions(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: pickconditions.m 1490 2008-04-28 11:16:29Z vladimir $

c = conditions(this);

res = strmatch(label, c);

res = res(~reject(this, res));
    
