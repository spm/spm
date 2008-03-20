function res = conditionlabels(obj, ind)
% Method for getting the unique labels of conditions in the file
% FORMAT res = conditionlabels(obj, ind)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: nconditions.m 1125 2008-01-30 12:12:18Z vladimir $

res = strvcat(unique({obj.trials.label}));

if exist('ind') == 1
    res = deblank(res(ind,:));
end