function [f] = spm_fx_mfm(x,u,P,M)
% state equations for a mean-field model
% FORMAT [f] = spm_fx_mfm(x,u,P,M)
% x      - state vector
% x{i,j} - states of i-th source and j-th population
% 
% x{i,1} - excitatory spiny stellate cells (input cells)
% x{i,2} - inhibitory interneurons
% x{i,3} - excitatory pyramidal cells      (output cells)
%
% x{i,j}.V  - voltage
% x{i,j}.gE - conductance (excitatory)
% x{i,j}.gI - conductance (inhibitory)
%
% f    = dx(t)/dt  = f(x(t))
%
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_mfm.m 1174 2008-02-27 20:22:30Z karl $
 
% get dimensions and configure state variables
%--------------------------------------------------------------------------
x     = spm_unvec(x,M.x);
ns    = size(x,1);                            % number of sources
np    = size(x,2);                            % number of populations
nc    = length(P.A);                          % number of connection types

% exponential transform to ensure positivity constraints
%--------------------------------------------------------------------------
A{1} = exp(P.A{1})*32;
A{2} = exp(P.A{2})*16;
A{3} = exp(P.A{3})*4;
C    = exp(P.C);

% time constants (ns x np) (excitatory, inhibitory)
%--------------------------------------------------------------------------
TE   =  4/1000*exp(P.T(:,1));                 % excitatory time constants
TI   = 16/1000*exp(P.T(:,2));                 % inhibitory time constants
TE   = kron([1 1 1],TE);
TI   = kron([1 1 1],TI);

% internal input (from last population x{:,end})
%--------------------------------------------------------------------------
VL    = -32;
VE    =  16;
VI    = -128;
VR    = exp(P.R(1));
S2    = exp(P.R(2))*256;
for i = 1:ns
    for j = 1:np
        m(i,j) = spm_Ncdf(x{i,j}.V,VR,128);
    end
end

% external input (to first population x{:,1})
%--------------------------------------------------------------------------
U     = kron([1 0 0],C*u);

% G - switches on extrinsic connections (np x nc)
%--------------------------------------------------------------------------
SE    = exp(P.G(:,1));
SI    = exp(P.G(:,2));
G     = [1 0 1;
         0 1 1;
         0 1 1];

% G - switches on intrinsic connections (np x np) - excitatory
%--------------------------------------------------------------------------
GE    = [0 0 1;
         0 0 1;
         1 0 0];

     
% G - switches on intrinsic connections (np x np) - inhibitory
%--------------------------------------------------------------------------
GI    = [0 0 0;
         0 0 0;
         0 1 0];


% flow
%==========================================================================
c     = 2/1000;
gL    = 1/32;
f     = x;
for i = 1:ns
    for j = 1:np

        % intrinsic coupling
        %------------------------------------------------------------------
        E     = 0;
        I     = 0;
        for k = 1:np
            E = E + SE(i)*GE(j,k)*m(i,k);
            I = I + SI(i)*GI(j,k)*m(i,k);
        end

        % extrinsic coupling
        %------------------------------------------------------------------
        for k = 1:ns
            for l = 1:nc
                E = E + A{l}(i,k)*G(j,l)*m(k,end);
            end
        end

        % Voltage
        %------------------------------------------------------------------
        f{i,j}.V =        gL*(VL - x{i,j}.V) + ...
                   x{i,j}.gE*(VE - x{i,j}.V) + ...
                   x{i,j}.gI*(VI - x{i,j}.V) + U(i,j);

        f{i,j}.V = f{i,j}.V/c;

        % Conductances
        %------------------------------------------------------------------
        f{i,j}.gE = (E - x{i,j}.gE)/TE(i,j);
        f{i,j}.gI = (I - x{i,j}.gI)/TI(i,j);
        
    end
end

% vectorise
%--------------------------------------------------------------------------
f = spm_vec(f);


