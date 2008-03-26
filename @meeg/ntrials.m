function res = ntrials(obj)
% Method for getting the number of trials in the file
% FORMAT res = ntrials(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: ntrials.m 1125 2008-01-30 12:12:18Z vladimir $

res = length(obj.trials);