function [x] = spm_x_nmm(P)
% intialises a state structure for a mean field model
% FORMAT [x] = spm_x_nmm(P)
%
% P - parameter structure
% x - state array

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
%__________________________________________________________________________

% dimensions
%--------------------------------------------------------------------------
ns   = size(P.A{1},1);                           % number of sources
np   = 3;                                        % number of populations

% create (initlaise voltage at -8mV)
%--------------------------------------------------------------------------
x        = zeros(ns,np,3);
x(:,:,1) = -8;


% steady-state solution 
%==========================================================================

% create LFP model
%--------------------------------------------------------------------------
M.f   = 'spm_fx_nmm';
M.x   = x;
M.pE  = P;
M.n   = length(spm_vec(x));
M.m   = size(P.C,2);
M.l   = size(P.C,1);

% solve for fixed point 
%--------------------------------------------------------------------------
U.u   = sparse(64,1);
U.dt  = 8/1000;
x     = spm_int_L(P,M,U);
x     = spm_unvec(x(end,:),M.x);
