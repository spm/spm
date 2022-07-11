function res = timeonset(this, newonset)
% Method for reading and setting the time onset
% FORMAT res = timeonset(this)
%        res = timeonset(this, newonset)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if nargin == 1
    res = this.timeOnset;
else
    this.timeOnset = newonset;
    res = this;
end
