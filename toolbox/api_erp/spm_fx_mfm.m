function [f] = spm_fx_mfm(x,u,P,M)
% state equations for neural-mass and mean-field models
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
% Marreiros et al (2008) Population dynamics under the Laplace assumption
%
% See also:
%
% Friston KJ.
% The labile brain. I. Neuronal transients and nonlinear coupling. Philos
% Trans R Soc Lond B Biol Sci. 2000 Feb 29;355(1394):215-36. 
% 
% McCormick DA, Connors BW, Lighthall JW, Prince DA.
% Comparative electrophysiology of pyramidal and sparsely spiny stellate
% neurons of the neocortex. J Neurophysiol. 1985 Oct;54(4):782-806.
% 
% Brunel N, Wang XJ.
% What determines the frequency of fast network oscillations with irregular
% neural discharges? I. Synaptic dynamics and excitation-inhibition
% balance. J Neurophysiol. 2003 Jul;90(1):415-30.
% 
% Brunel N, Wang XJ.
% Effects of neuromodulation in a cortical network model of object working
% memory dominated by recurrent inhibition. J Comput Neurosci. 2001
% Jul-Aug;11(1):63-85.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_mfm.m 1277 2008-03-28 18:36:49Z karl $
 
% get dimensions and configure state variables
%--------------------------------------------------------------------------
if iscell(x)
    mfm = 1;                                    % mean-field model
else
    mfm = 0;
    x = {x};                                    % neural-mass model
end
ns   = size(x{1},1);                            % number of sources
np   = size(x{1},2);                            % number of populations
nc   = length(P.A);                             % number of connections
 
 
% extrinsic connection strengths
%==========================================================================
 
% exponential transform to ensure positivity constraints
%--------------------------------------------------------------------------
A{1} = exp(P.A{1})/2;                         % forward
A{2} = exp(P.A{2})/4;                         % backward
A{3} = exp(P.A{3})/4;                         % lateral
C    = exp(P.C);                              % subcortical
 
% switches on extrinsic afferent connections (np x nc)
%--------------------------------------------------------------------------
SA   = sparse([1 0 1;
               0 1 1;
               0 1 1]);
            
% intrinsic connection strengths
%==========================================================================
G    = exp(P.G);
try
    
    % get intrinsic connections from model structure
    %----------------------------------------------------------------------
    GE = M.GE;
    GI = M.GI;
    
catch

    % intrinsic connections (np x np) - excitatory
    %----------------------------------------------------------------------
    GE   = [0   0   1/2;
            0   0   1;
            1   0   0  ];

    % intrinsic connections (np x np) - inhibitory
    %----------------------------------------------------------------------
    GI   = [0   1/2 0;
            0   0   0;
            0   2   0];
end
                
 
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
fxx  = sparse([2 3 1 1],[1 1 2 3],-1/CV);    % curvature: df(V)/dxx
 
% mean-field effects
%==========================================================================
if mfm
    
    Vx  = shiftdim(x{2}(1,1,:,:));           % population variance (mV^2)
    Vx  = reshape(Vx,ns,np);                 % of voltage
 
    D   = sparse(diag([1/32 1 1]));          % diffusion
    D   = exp(P.S)*D;
    
else
    
    % neural-mass approximation to covariance of states
    %----------------------------------------------------------------------
    try
        Cx = M.Cx;
    catch
        Cx = [75    0.2    0.8;
              0.2   0.004  0  ;
              0.8   0      0.02];
    end
    Cx = exp(P.S)*Cx;
    Vx = Cx(1,1);  
    
end
 
% mean population firing and afferent extrinsic input
%--------------------------------------------------------------------------
m     = spm_Ncdf_jdw(x{1}(:,:,1),VR,Vx);  
for k = 1:nc
    a(:,k) = A{k}*m(:,end);
end
 
% Exogenous input (to first population x{:,1})
%--------------------------------------------------------------------------
U     = C*u;
 
% flow and dispersion over every (ns x np) subpopulation
%==========================================================================
f     = x;
for i = 1:ns
    for j = 1:np
 
        % 1st moment - expected states
        %==================================================================
        
        % intrinsic coupling
        %------------------------------------------------------------------
        E = G(i)*GE(j,:)*m(i,:)';
        I =      GI(j,:)*m(i,:)';
        
        % extrinsic coupling (excitatory only)
        %------------------------------------------------------------------
        E =  E + SA(j,:)*a(i,:)';
 
        % Voltage
        %------------------------------------------------------------------
        f{1}(i,j,1) =         (GL*(VL - x{1}(i,j,1)) + ...
                      x{1}(i,j,2)*(VE - x{1}(i,j,1)) + ...
                      x{1}(i,j,3)*(VI - x{1}(i,j,1)) )/CV;
                  
        % Exogenous input (U/)
        %------------------------------------------------------------------
        if j == 1
            f{1}(i,j,1) = f{1}(i,j,1) + U(i)/CV;
        end
 
        % Conductances
        %------------------------------------------------------------------
        f{1}(i,j,2) = (E - x{1}(i,j,2))*KE(i);
        f{1}(i,j,3) = (I - x{1}(i,j,3))*KI;
        
        
        % 2nd moments - covariances
        %==================================================================
        if mfm
            
            % add curvature-dependent dispersion to flow
            %--------------------------------------------------------------
            f{1}(i,j,1) = f{1}(i,j,1) + tr(x{2}(:,:,i,j),fxx)/2;
            
            % df/dx
            %--------------------------------------------------------------
            Sg  = -GL - x{1}(i,j,2) - x{1}(i,j,3);
            fx  = [Sg/CV (VE - x{1}(i,j,1))/CV (VI - x{1}(i,j,1)/CV);
                   0     -KE(i)              0                ;
                   0      0                 -KI              ];
        
            % dCdt
            %--------------------------------------------------------------
            St            = fx*x{2}(:,:,i,j) + D;
            f{2}(:,:,i,j) = St + St';
            
        else
            % fixed covariance (Cx)
            %--------------------------------------------------------------
            f{1}(i,j,1) = f{1}(i,j,1) + tr(Cx,fxx)/2;
            
        end
        
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
