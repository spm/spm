function [f] = spm_fx_mfm_ensemble(x,u,P)
% state equations for a mean-field model
% FORMAT [f] = spm_fx_mfm_ensemble(x,u,P)
%
% X{i} - states and covariances of i-th particle x
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
% $Id: spm_fx_mfm_ensemble.m 1212 2008-03-14 19:08:47Z karl $
 
% get dimensions and configure state variables
%--------------------------------------------------------------------------
nn   = size(x,1);                               % number of neurons
ns   = size(x,2);                               % number of sources
np   = size(x,3);                               % number of populations
nc   = length(P.A);                             % number of connection types


% extrinsic connection strengths
%==========================================================================

% exponential transform to ensure positivity constraints
%--------------------------------------------------------------------------
A{1} = exp(P.A{1})/8;                         % forward
A{2} = exp(P.A{2})/16;                        % backward
A{3} = exp(P.A{3})/64;                        % lateral
C    = exp(P.C);                              % subcortical
 
% switches on extrinsic afferent connections (np x nc)
%--------------------------------------------------------------------------
SA   = sparse([1 0 1;
               0 1 1;
               0 1 1]);
            
% intrinsic connection strengths
%==========================================================================
G    = exp(P.G);
 
% intrinsic connections (np x np) - excitatory
%--------------------------------------------------------------------------
SE   = sparse([0   0   1/2;
               0   0   1;
               1   0   0  ]);
     
% intrinsic connections (np x np) - inhibitory
%--------------------------------------------------------------------------
SI   = sparse([0   1/2 0;
               0   0   0;
               0   2   0]);
                
 
% rate constants (ns x np) (excitatory 8ms, inhibitory 16ms)
%--------------------------------------------------------------------------
KE   = exp(-P.T)*1000/4;                     % excitatory time constants
KI   = 1000/16;                              % inhibitory time constants
 
% Voltages
%--------------------------------------------------------------------------
VL   = -70;                                  % reversal  potential leak (K)
VE   =  60;                                  % reversal  potential excite (Na)
VI   = -90;                                  % reversal  potential inhib (Cl)
VR   = -40;                                  % threshold potential
 
CV   = 8/1000;                               % membrane capacitance
GL   = 1;                                    % leak conductance

% mean population firing rate
%--------------------------------------------------------------------------
m    = mean(spm_Ncdf_jdw(x(:,:,:,1),VR,exp(-8)),1);
m    = squeeze(m);
m    = reshape(m,ns,np);

% mean afferent extrinsic input
%--------------------------------------------------------------------------               
for k = 1:nc
    a(:,k) = A{k}*m(:,end);
end

% external input (to first population x{:,1})
%--------------------------------------------------------------------------
U     = kron(sparse(1,1,1,1,np),C*u);


% flow and dispersion over every (ns x np) subpopulation
%==========================================================================
f     = x;                                  % flow
for i = 1:ns
    for j = 1:np
        for k = 1:nn
                        
            % states
            %==============================================================

            % intrinsic coupling
            %--------------------------------------------------------------
            E = G(i)*SE(j,:)*m(i,:)';
            I =      SI(j,:)*m(i,:)';

            % extrinsic coupling (excitatory only)
            %--------------------------------------------------------------
            E =  E + SA(j,:)*a(i,:)';

            % Voltage
            %--------------------------------------------------------------
            f(k,i,j,1) =         GL*(VL - x(k,i,j,1)) + ...
                         x(k,i,j,2)*(VE - x(k,i,j,1)) + ...
                         x(k,i,j,3)*(VI - x(k,i,j,1)) + U(i,j);

            f(k,i,j,1) = f(k,i,j,1)/CV;

            % Conductances
            %--------------------------------------------------------------
            f(k,i,j,2) = (E - x(k,i,j,2))*KE(i);
            f(k,i,j,3) = (I - x(k,i,j,3))*KI;
            
        end
    end
end

% vectorise
%--------------------------------------------------------------------------
f = spm_vec(f);

return


