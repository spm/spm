function [x,M] = spm_x_mfm(P,GE,GI)
% initialises a state structure for a mean field model
% FORMAT [x,M] = spm_x_mfm(P,GE,GI)
%
% P  - parameter structure (encoding extrinsic connections)
% GE - extrinsic connections (excitatory)
% GI - extrinsic connections (inhibitory)
%
%
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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_x_mfm.m 2310 2008-10-06 19:20:45Z karl $
 
 
% intrinsic connections (specifying the number of populations per source)
%==========================================================================
try, GE; catch
 
    % intrinsic connections (np x np) - excitatory
    %----------------------------------------------------------------------
    GE   = [0   0   1/2;
            0   0   1;
            1   0   0  ];
end
try, GI; catch
 
    % intrinsic connections (np x np) - inhibitory
    %----------------------------------------------------------------------
    GI   = [0   1/2 0;
            0   0   0;
            0   2   0];
end
 
 
% dimensions
%--------------------------------------------------------------------------
ns   = size(P.A{1},1);                           % number of sources
np   = size(GE,1);                               % number of populations
 
% create (initialise voltage at -70mV)
%--------------------------------------------------------------------------
x{1}        = zeros(ns,np,3);
x{1}(:,:,1) = -70;
x{2}        = zeros(3,3,ns,np);
for i = 1:ns
    for j = 1:np
        x{2}(:,:,i,j) = eye(3,3)/128;
    end
end
 
% steady-state solution 
%==========================================================================
 
% create MFM model
%--------------------------------------------------------------------------
M.f   = 'spm_fx_mfm';
M.x   = x;
M.pE  = P;
M.n   = length(spm_vec(x));
M.m   = size(P.C,2);
M.l   = size(P.C,1);
M.GE  = GE;
M.GI  = GI;
 
% solve for fixed point 
%--------------------------------------------------------------------------
U.u   = sparse(16,1);
U.dt  = 16/1000;
x     = spm_int_ode(P,M,U);
x     = spm_unvec(x(end,:),M.x);
M.x   = x;
