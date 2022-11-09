function [x,M] = spm_x_nmm(P)
% Initialise a state structure for a mean field model
% FORMAT [x,M] = spm_x_nmm(P)
%
% P - parameter structure
% M - model structure
%
% x        - array of states
% x(i,j,k) - k-th state of j-th population on i-th source
%
%   population: 1 - excitatory spiny stellate cells (input cells)
%               2 - inhibitory interneurons
%               3 - excitatory pyramidal cells      (output cells)
%
%        state: 1 V  - voltage
%               2 gE - conductance (excitatory)
%               3 gI - conductance (inhibitory)
%
% M - model structure
%
% see also: spm_x_mfm
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


% get initialisation from full mean-field model
%==========================================================================
[x,M] = spm_x_mfm(P);
 
% remove dispersion and fix the covariance of the states (Cx)
%--------------------------------------------------------------------------
M.x   = x{1};
M.Cx  = x{2}(:,:,1,1);
x     = x{1};
