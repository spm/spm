function [x,M] = spm_x_nmm(P,M)
% initialises a state structure for a mean field model
% FORMAT [x,M] = spm_x_nmm(P,M)
%
% P - parameter structure
% M - model     structure
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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_x_nmm.m 2330 2008-10-10 18:23:42Z karl $
 

% get initialisation from full mean-field model
%==========================================================================
try
    [x M] = spm_x_mfm(P,M);
catch
    [x M] = spm_x_mfm(P);
end
 
% remove dispersion and fix the covariance of the states (Cx)
%--------------------------------------------------------------------------
M.x   = x{1};
M.Cx  = x{2}(:,:,1,1);
x     = x{1};
