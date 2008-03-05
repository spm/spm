function [x,M] = spm_x_mfm(P)
% intialises a state structure for a mean field model
% FORMAT [x,M] = spm_x_mfm(P)
%
% P - parameter structure
% x - states and covariances
% M - model structure
%
% x{1}(i,j,k)   - k-th state of i-th source in j-th population
% x{2}(i,j,k,l) - covariance of i-th and j-th state (k-th source in l-th
%                 population
%
%   population: 1 - excitatory spiny stellate cells (input cells)
%               2 - inhibitory interneurons
%               3 - excitatory pyramidal cells      (output cells)
%
%        state: 1 V  - voltage
%               2 gE - conductance (excitatory)
%               3 gI - conductance (inhibitory)
%__________________________________________________________________________

% dimensions
%--------------------------------------------------------------------------
ns   = size(P.A{1},1);                           % number of sources
np   = 3;                                        % number of populations

% create (initlaise voltage at -8mV)
%--------------------------------------------------------------------------
x{1}        = zeros(ns,np,3);
x{1}(:,:,1) = -32;
x{2}        = zeros(3,3,ns,np);
for i = 1:ns
    for j = np
        x{2}(:,:,i,j) = eye(3,3)/128;
    end
end

% steady-state solution 
%==========================================================================

% create LFP model
%--------------------------------------------------------------------------
M.f   = 'spm_fx_mfm';
M.x   = x;
M.pE  = P;
M.n   = length(spm_vec(x));
M.m   = size(P.C,2);
M.l   = size(P.C,1);

% solve for fixed point 
%--------------------------------------------------------------------------
U.u   = sparse(64,1);
U.dt  = 8/1000;
x     = spm_int_ode(P,M,U);
x     = spm_unvec(x(end,:),M.x);
M.x   = x;

