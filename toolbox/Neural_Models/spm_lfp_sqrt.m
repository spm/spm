function [y] = spm_lfp_sqrt(y,M)
% Feature selection for lfp and mtf (spectral) neural mass models
% FORMAT [y] = spm_lfp_sqrt(y,M)
% 
% Y -> log(y) (including cells)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_lfp_sqrt.m 1132 2008-02-06 14:12:17Z karl $

% log transform
%--------------------------------------------------------------------------
y = spm_unvec(sqrt(spm_vec(y)),y);