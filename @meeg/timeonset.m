function res = timeonset(this, newonset)
% Method for reading and setting the time onset
% FORMAT res = timeonset(this)
%        res = timeonset(this, newonset)
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: timeonset.m 5025 2012-10-31 14:44:13Z vladimir $

if nargin == 1
    res = this.timeOnset;
else
    this.timeOnset = newonset;
    res = this;
end