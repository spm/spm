function [y] = spm_fs_csd(y,M)
% observer for a CSD DCM
% FORMAT [y] = spm_fs_csd(y,M)
% y      - CSD
% M      - model structure
%__________________________________________________________________________
%
% This simple removes the (redundant) lower triangular part of the cross-
% spectra
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fs_csd.m 4261 2011-03-24 16:39:42Z karl $
 
% parameterised lead field times [perturbations] of states
%--------------------------------------------------------------------------
y = y;
 
