function [coh,fsd] = spm_ker2coh(ker,pst)
% computes coherence from kernels
% FORMAT [coh,fsd] = spm_ker2coh(ker,pst))
%
% ker  - first-order (Volterra) kernels
% pst  - time samples
%
% coh  - coherence
% fsd  - frequency specific delay (seconds) 
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ker2coh.m 7774 2020-01-25 18:07:03Z karl $


% coherence
%==========================================================================

% via modulation transfer function
%--------------------------------------------------------------------------
[mtf,Hz]  = spm_ker2mtf(ker,pst);
[coh,fsd] = spm_mtf2coh(mtf,Hz);