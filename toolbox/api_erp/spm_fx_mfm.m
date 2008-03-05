function [f] = spm_fx_mfm(x,u,P,M)
% state equations for a mean-field model
% FORMAT [f] = spm_fx_mfm(x,u,P,M)
%
% x - states and covariances
%
% x{1}(i,j,k)   - k-th state of j-th population on i-th source
% x{2}(i,j,k,l) - covariance of l-th and k-th state
%
%   population: 1 - excitatory spiny stellate cells (input cells)
%               2 - inhibitory interneurons
%               3 - excitatory pyramidal cells      (output cells)
%
%        state: 1 V  - voltage
%               2 gE - conductance (excitatory)
%               3 gI - conductance (inhibitory)
%
%--------------------------------------------------------------------------
% refs:
%
% Marreiros et al (2008) Population dynamics under the Laplac assumption 
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_mfm.m 1190 2008-03-05 17:20:23Z karl $
 
% get dimensions and configure state variables
%--------------------------------------------------------------------------
ns   = size(x{1},1);                            % number of sources
np   = size(x{1},2);                            % number of populations
nc   = length(P.A);                             % number of connection types


% extrinsic connection strengths
%==========================================================================

% exponential transform to ensure positivity constraints
%--------------------------------------------------------------------------
A{1} = exp(P.A{1})/8;                         % forward
A{2} = exp(P.A{2})/16;                        % backward
A{3} = exp(P.A{3})/64;                        % lateral
C    = exp(P.C)/16;                           % subcortical

% switches on extrinsic afferent connections (np x nc)
%--------------------------------------------------------------------------
SA   = sparse([1 0 1;
               0 1 1;
               0 1 1]);
            
% intrinsic connection strengths
%==========================================================================
G    = exp(P.G);

% switches on intrinsic connections (np x np) - excitatory
%--------------------------------------------------------------------------
SE   = sparse([0 0 1;
               0 0 1;
               1 0 0]);
     
% switches on intrinsic connections (np x np) - inhibitory
%--------------------------------------------------------------------------
SI   = sparse([0 0 0;
               0 0 0;
               0 1 0]);
                

% rate constants (ns x np) (excitatory 8ms, inhibitory 16ms)
%--------------------------------------------------------------------------
KE   = exp(-P.T)*1000/8;                     % excitatory time constants
KI   = 1000/16;                              % inhibitory time constants

% internal input
%--------------------------------------------------------------------------
VL   = -32;                                 % reversal  potential leak
VE   =  16;                                 % reversal  potential excite
VI   = -128;                                % reversal  potential inhib
VR   = 0;                                   % threshold potential
VP   = exp(P.S)*128;                        % population variance (mV^2)
m    = spm_Ncdf_jdw(x{1}(:,:,1),VR,VP);     % mean population firing rate

% afferent extrinsic input
%--------------------------------------------------------------------------               
for k = 1:nc
    a(:,k) = A{k}*m(:,end);
end

% external input (to first population x{:,1})
%--------------------------------------------------------------------------
U     = kron(sparse(1,1,1,1,np),C*u);


% flow and dispersion over every (ns x np) subpopulation
%==========================================================================
CV    = 1/1000;                             % membrane capacitance
gL    = 1/16;                               % leak condunctance                        
fxx   = sparse([2 3 1 1],[1 1 2 3],-1/CV);  % df(V)/dxx
D     = sparse(diag([1/32 1/128 1/128]));   % diffusion
f     = x;                                  % flow
for i = 1:ns
    for j = 1:np

        % 1st moment - expected states
        %==================================================================
        
        % intrinsic coupling
        %------------------------------------------------------------------
        E = G(i)*SE(j,:)*m(i,:)';
        I =      SI(j,:)*m(i,:)';
        
        % extrinsic coupling (excitatory only)
        %------------------------------------------------------------------
        E =  E + SA(j,:)*a(i,:)';

        % Voltage
        %------------------------------------------------------------------
        f{1}(i,j,1) =          gL*(VL - x{1}(i,j,1)) + ...
                      x{1}(i,j,2)*(VE - x{1}(i,j,1)) + ...
                      x{1}(i,j,3)*(VI - x{1}(i,j,1)) + U(i,j);

        f{1}(i,j,1) = f{1}(i,j,1)/CV + tr(x{2}(:,:,i,j),fxx)/2;

        % Conductances
        %------------------------------------------------------------------
        f{1}(i,j,2) = (E - x{1}(i,j,2))*KE(i);
        f{1}(i,j,3) = (I - x{1}(i,j,3))*KI;
        
        % 2nd moments - covariances
        %==================================================================
        
        % df/dx
        %------------------------------------------------------------------
        Sg  = gL + x{1}(i,j,2) + x{1}(i,j,3);
        fx  = [-Sg/CV (VE - x{1}(i,j,1))/CV (VI - x{1}(i,j,1)/CV);
               0      -KE(i)              0                ;
               0       0                 -KI              ];
        
        % dCdt
        %------------------------------------------------------------------
        St            = fx*x{2}(:,:,i,j) + D;
        f{2}(:,:,i,j) = St + St';
        
    end
end

% vectorise
%--------------------------------------------------------------------------
f = spm_vec(f);

return

% trace(a*b)
%--------------------------------------------------------------------------
function x = tr(a,b);
%__________________________________________________________________________
b   = b';
x   = a(:)'*b(:);


