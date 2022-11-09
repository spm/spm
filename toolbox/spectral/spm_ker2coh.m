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

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% coherence
%==========================================================================

% via modulation transfer function
%--------------------------------------------------------------------------
[mtf,Hz]  = spm_ker2mtf(ker,pst);
[coh,fsd] = spm_mtf2coh(mtf,Hz);