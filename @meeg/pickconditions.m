function res = pickconditions(this, label, rejectbad)
% Method for returning indices of trials of a certain trial type.
% If input argument rejectbad == 0, the function will also return trial 
% indices which are set to bad (i.e. rejected). If bad is omitted, 
% the default is not to include rejected trials.
% FORMAT res = pickconditions(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: pickconditions.m 2448 2008-11-07 16:56:07Z vladimir $

if nargin<3
    rejectbad = 1;
end

c = conditions(this);

res = strmatch(deblank(label), deblank(c), 'exact');
    
if rejectbad && ~isempty(res)
    res = res(~reject(this, res));
end