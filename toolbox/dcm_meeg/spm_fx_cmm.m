function [f,J,Q] = spm_fx_cmm(x,u,P,M)
% state equations for canonical neural-mass and mean-field models
% FORMAT [f,J,Q] = spm_fx_cmm(x,u,P,M)
%
% x - states and covariances
%
% x(i,j,k)        - k-th state of j-th population of i-th source
%                   i.e., running over sources, pop. and states
%
%   population: 1 - excitatory spiny stellate cells (input cells)
%               2 - superficial pyramidal cells     (forward output cells)
%               3 - inhibitory interneurons         (intrisic interneuons)
%               4 - deep pyramidal cells            (backward output cells)
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
% $Id: spm_fx_cmm.m 5019 2012-10-26 19:32:57Z karl $
 
% get dimensions and configure state variables
%--------------------------------------------------------------------------
ns   = size(M.x,1);                      % number of sources
np   = size(M.x,2);                      % number of populations per source
nk   = size(M.x,3);                      % number of states per population
x    = reshape(x,ns,np,nk);              % hidden states 


% extrinsic connection strengths
%==========================================================================
 
% exponential transform to ensure positivity constraints
%--------------------------------------------------------------------------
A{1} = exp(P.A{1});                      % forward
A{2} = exp(P.A{2});                      % backward
C    = exp(P.C);                         % subcortical
 

% detect and reduce the strength of reciprocal (lateral) connections
%--------------------------------------------------------------------------
for i = 1:length(A)
    L    = (A{i} > exp(-8)) & (A{i}' > exp(-8));
    A{i} = A{i}./(1 + 8*L);
end

            
% intrinsic connection strengths
%==========================================================================

% condition specific effects
%--------------------------------------------------------------------------
G        = full(P.G);
G(2,2,:) = squeeze(G(2,2,:)) + P.H;
G        = exp(G);

% connectivity switches
%==========================================================================
% 1 - excitatory spiny stellate cells (granular input cells)
% 2 - superficial pyramidal cells     (forward  output cells)
% 3 - inhibitory interneurons         (intrisic interneuons)
% 4 - deep pyramidal cells            (backward output cells)

% extrinsic connections (F B) - from superficial and deep pyramidal cells
%--------------------------------------------------------------------------
SA   = [2   0  ;
        0   0  ;
        0   16  ;
        0   0];
 
% intrinsic connections (np x np) - excitatory
%--------------------------------------------------------------------------
GE   = [ 0     0     0     0
        96     0     0     0
         4     0     0   128
         0     2     0     0];
 
% intrinsic connections (np x np) - inhibitory
%--------------------------------------------------------------------------
GI   = [12     0     4     0
         0  2048     0     0
         0     0    64     0
         0     2    96     2];

% rate constants (ns x np) (excitatory 4ms, inhibitory 16ms)
%--------------------------------------------------------------------------
KE   = [1/2 1/2 1/2 1/2]*1000*4;                % excitatory rate constants
KI   = exp(-P.T)*[1/32 1/16 1/16 1/16]*1000;  % inhibitory rate constants
 
% Voltages
%--------------------------------------------------------------------------
VL   = -70;                               % reversal  potential leak (K)
VE   =  60;                               % reversal  potential excite (Na)
VI   = -90;                               % reversal  potential inhib (Cl)
VR   = -40;                               % threshold potential
 
CV   = exp(P.CV).*[24 60 42 10]/1000;     % membrane capacitance
GL   = 1;                                 % leak conductance
 
% mean-field effects:
%==========================================================================

% neural-mass approximation to covariance of states
%----------------------------------------------------------------------
Vx   = exp(P.S)*16;

% mean population firing and afferent extrinsic input
%--------------------------------------------------------------------------
m      = spm_Ncdf_jdw(x(:,:,1),VR,Vx);    % mean firing rate  
a(:,1) = A{1}*m(:,2);                     % forward afference
a(:,2) = A{2}*m(:,4);                     % backward afference

% Averge background activity and exogenous input
%==========================================================================
BE     = exp(P.E)*0.21;

% input
%--------------------------------------------------------------------------
if isfield(M,'u')
    
    % endogenous input
    %----------------------------------------------------------------------
    U = u(:)/32;
    
else
    
    % exogenous input
    %----------------------------------------------------------------------
    U = C*u(:)/256;
    
end

% flow over every (ns x np) subpopulation
%==========================================================================
f     = x;
for i = 1:ns
    for j = 1:np
 
        % intrinsic coupling
        %------------------------------------------------------------------
        E = (G(j,:,i).*GE(j,:))*m(i,:)';
        I = (G(j,:,i).*GI(j,:))*m(i,:)';
        
        % extrinsic coupling (excitatory only) and background activity
        %------------------------------------------------------------------
        E = E + BE + SA(j,:)*a(i,:)';

        % and exogenous input(U)
        %------------------------------------------------------------------
        if j == 1, E = E + U(i); end
        
        % Voltage
        %==================================================================
        f(i,j,1) =         (GL*(VL - x(i,j,1)) + ...
                      x(i,j,2)*(VE - x(i,j,1)) + ...
                      x(i,j,3)*(VI - x(i,j,1)) )/CV(j);
        
        % Conductance
        %==================================================================
        f(i,j,2) = (E - x(i,j,2))*KE(1,j);
        f(i,j,3) = (I - x(i,j,3))*KI(i,j);
        
    end
end
 
% vectorise equations of motion
%==========================================================================
f = spm_vec(f);
 
if nargout < 2, return, end

% Jacobian
%==========================================================================
J = spm_cat(spm_diff(M.f,x,u,P,M,1));

if nargout < 3, return, end

% Delays
%==========================================================================
% Delay differential equations can be integrated efficiently (but 
% approximately) by absorbing the delay operator into the Jacobian
%
%    dx(t)/dt     = f(x(t - d))
%                 = Q(d)f(x(t))
%
%    J(d)         = Q(d)df/dx
%--------------------------------------------------------------------------
% [specified] fixed parameters
%--------------------------------------------------------------------------
D  = [2 16];

d  = -D.*exp(P.D)/1000;
Sp = kron(ones(nk,nk),kron( eye(np,np),eye(ns,ns)));  % states: same pop.
Ss = kron(ones(nk,nk),kron(ones(np,np),eye(ns,ns)));  % states: same source

Dp = ~Ss;                            % states: different sources
Ds = ~Sp & Ss;                       % states: same source different pop.
D  = d(2)*Dp + d(1)*Ds;


% Implement: dx(t)/dt = f(x(t - d)) = inv(1 - D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%--------------------------------------------------------------------------
Q  = spm_inv(speye(length(J)) - D.*J);


