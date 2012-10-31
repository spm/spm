function res = pickconditions(this, label, rejectbad)
% Method for returning indices of trials of a certain trial type.
% If input argument rejectbad == 0, the function will also return trial 
% indices which are set to bad (i.e. rejected). If bad is omitted, 
% the default is not to include rejected trials.
% FORMAT res = pickconditions(this)
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: pickconditions.m 5025 2012-10-31 14:44:13Z vladimir $

warning_flexible('pickconditions method is deprecated. Use ''indrial'' instead');

if nargin<3
    rejectbad = 1;
end

if isa(label, 'char')
    label = {label};
end

if rejectbad
    label{end+1} = 'GOOD';
end

res = indtrial(this, label);