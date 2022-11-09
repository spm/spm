function [y] = spm_lfp_sqrt(y,M)
% Feature selection for lfp and mtf (spectral) neural mass models
% FORMAT [y] = spm_lfp_sqrt(y,M)
% 
% Y -> log(y) (including cells)
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% log transform
%--------------------------------------------------------------------------
y = spm_unvec(sqrt(spm_vec(y)),y);
