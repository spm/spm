function res = pickconditions(this, label)
% Method for returning indices of trials of a certain trial type.
% note that this function will also return the 'bad' trials.
% FORMAT res = pickconditions(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: pickconditions.m 2433 2008-11-03 16:26:13Z stefan $

c = conditions(this);

res = strmatch(label, c);
    
