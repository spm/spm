function res = timeonset(this, value)
% Method for and setting the time onset
% FORMAT res = timeonset(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: timeonset.m 1254 2008-03-27 18:41:42Z vladimir $

if nargin == 1
    res = this.timeOnset;
else
    this.timeOnset = value;
    res = this;
end