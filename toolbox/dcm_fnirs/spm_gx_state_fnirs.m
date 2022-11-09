function [y] = spm_gx_state_fnirs(x,u,P,M)
% Neurodynamics and Hemodynamics underling DCM for fNIRS
% FORMAT [y] = spm_gx_state_fnirs(x,u,P,M)
%
% x   - state vector     (see spm_fx_fnirs)
% u   - experimental inputs 
% P   - prior of latent variables 
% M   - model structure
%
% y   - fNIRS response and copied state vector
%
% The `copied state vector' passes the first hidden variable in each region
% to the output variable y, so that neuronal activity and state variables
% can be plotted. 
%__________________________________________________________________________

% Will Penny & Sungho Tak
% Copyright (C) 2015-2022 Wellcome Centre for Human Neuroimaging
 
y = spm_gx_fnirs(x,u,P,M);

% y: optical density changes 
% x(:,1): neuronal activity 
% x(:,5): deoxyHb
% x(:,6): totalHb 
y = [y; x(:,1); x(:,5); x(:,6)]; 
y = full(y);
