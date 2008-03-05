function [f] = spm_fx_nmm(x,u,P,M)
% state equations for a mean-field model
% FORMAT [f] = spm_fx_nmm(x,u,P,M)
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
% f    = dx(t)/dt  = f(x(t))
%
%--------------------------------------------------------------------------
% refs:
%
% An Approximation to the Probability Integral
% J. D. Williams 
% The Annals of Mathematical Statistics, Vol. 17, No. 3. (Sep., 1946), pp. 363-365. 
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_nmm.m 1190 2008-03-05 17:20:23Z karl $
 
% get dimensions and configure state variables
%--------------------------------------------------------------------------
ns   = size(x,1);                            % number of sources
np   = size(x,2);                            % number of populations
nc   = length(P.A);                          % number of connection types


% extrinsic connection strengths
%==========================================================================

% exponential transform to ensure positivity constraints
%--------------------------------------------------------------------------
A{1} = exp(P.A{1})/2;                         % forward
A{2} = exp(P.A{2})/8;                         % backward
A{3} = exp(P.A{3})/32;                        % lateral
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
m    = spm_Ncdf_jdw(x(:,:,1),VR,VP);        % mean population firing rate

% afferent extrinsic input
%--------------------------------------------------------------------------               
for k = 1:nc
    a(:,k) = A{k}*m(:,end);
end

% external input (to first population x{:,1})
%--------------------------------------------------------------------------
U     = kron(sparse(1,1,1,1,np),C*u);


% flow
%==========================================================================
CV    = 1/1000;                             % membrane capacitance
gL    = 1/16;                               % leak condunctance                        
f     = x;
for i = 1:ns
    for j = 1:np

        % channel opening (e - excitatory and I - inhibitory)
        %==================================================================
        
        % intrinsic coupling
        %------------------------------------------------------------------
        E = G(i)*SE(j,:)*m(i,:)';
        I =      SI(j,:)*m(i,:)';
        
        % extrinsic coupling (excitatory only)
        %------------------------------------------------------------------
        E =  E + SA(j,:)*a(i,:)';

        % Voltage
        %==================================================================
        f(i,j,1) =       gL*(VL - x(i,j,1)) + ...
                   x(i,j,2)*(VE - x(i,j,1)) + ...
                   x(i,j,3)*(VI - x(i,j,1)) + U(i,j);

        f(i,j,1) = f(i,j,1)/CV;

        % Conductances
        %==================================================================
        f(i,j,2) = (E - x(i,j,2))*KE(i);
        f(i,j,3) = (I - x(i,j,3))*KI;
        
    end
end

% vectorise
%--------------------------------------------------------------------------
f = spm_vec(f);


