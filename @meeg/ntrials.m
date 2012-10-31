function res = ntrials(this)
% Method for getting the number of trials in the file
% FORMAT res = ntrials(this)
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: ntrials.m 5025 2012-10-31 14:44:13Z vladimir $

res = length(this.trials);