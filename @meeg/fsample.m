function res = fsample(this, value)
% Method for getting and setting the sampling rate
% FORMAT res = fsample(this)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if nargin == 1
    res = this.Fsample;
else
    this.Fsample = value;
    res = this;
end