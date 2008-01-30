function res = conditions(obj)
% Method for getting the condition labels
% FORMAT res = conditions(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
%$Id: conditions.m 1125 2008-01-30 12:12:18Z vladimir $

if obj.Nsamples>0
    res = unique({obj.trials.label});
else
    res = [];
end