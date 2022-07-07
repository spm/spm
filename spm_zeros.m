function [X] = spm_zeros(X)
% Fill a cell or structure array with zeros
% FORMAT [X] = spm_zeros(X)
% X  - numeric, cell or structure array[s]
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


% create zeros structure
%--------------------------------------------------------------------------
X = spm_unvec(zeros(spm_length(X),1),X);
