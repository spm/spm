function [X] = spm_zeros(X)
% fills a cell or structure array with zeros
% FORMAT [X] = spm_zeros(X)
% X  - numeric, cell or structure array[s]
%__________________________________________________________________________
% Copyright (C) 2005-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_zeros.m 8183 2021-11-04 15:25:19Z guillaume $


% create zeros structure
%--------------------------------------------------------------------------
X = spm_unvec(zeros(spm_length(X),1),X);
