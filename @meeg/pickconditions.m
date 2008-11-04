function res = pickconditions(this, label, bad)
% Method for returning indices of trials of a certain trial type.
% If input argument bad == 1, the function will also return trial 
% indices which are set to bad (i.e. rejected). If bad is omitted, 
% the default is to not include rejected trials.
% FORMAT res = pickconditions(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: pickconditions.m 2438 2008-11-04 11:21:19Z stefan $

if ~exist('bad', 'var')
    bad = 1;
end

c = conditions(this);

res = strmatch(label, c, 'exact');
    
if bad & ~isempty(res)
    res = res(~reject(this, res));
end