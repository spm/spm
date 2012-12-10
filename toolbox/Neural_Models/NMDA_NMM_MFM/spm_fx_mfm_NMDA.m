function [f,J,Q] = spm_fx_mfm_NMDA(x,u,P,M)
% state equations for neural-mass and mean-field models
% FORMAT [f,J,Q] = spm_fx_mfm_NMDA(x,u,P,M)
%
% x - states and covariances
%
% x{1}(i,j,k)   - k-th state of j-th population on i-th source
%                 i.e., running over sources, pop. and states
% x{2}(:,:,i,j) - covariance among k states
%                 i.e., running over states x states, sources and pop.
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
% $Id: spm_fx_mfm_NMDA.m 5106 2012-12-10 17:13:34Z rosalyn $
 
% get dimensions and configure state variables
%--------------------------------------------------------------------------

try, x = spm_unvec(x,M.x); end

xin = x;
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
C    = exp(P.C);                            % subcortical

% switches on extrinsic afferent connections (np x nc)
%--------------------------------------------------------------------------
try
    SA   = P.SA;
catch
    SA   = sparse([1 0 1;
                   0 1 1;
                   0 0 0]);
end
            
% intrinsic connection strengths
%==========================================================================
G    = exp(P.G);


% intrinsic connections (np x np) - excitatory

%----------------------------------------------------------------------

GE   = [0   0   1/2*exp(P.scale_geP1);
    0   0   1;
    1/2*exp(P.scale_geS) 0   exp(-5)*exp(P.scale_geP2)];

GE_NMDA   = [0   0   0;
    0   0   1;
    1*exp(P.scale_geS) 0   exp(-5)*exp(P.scale_geP2)];


% intrinsic connections (np x np) - inhibitory
%----------------------------------------------------------------------
GI   = [0   1/4 0;
    0   0    0;
    0   1*exp(P.scale_gi)  0];


    
% rate constants (ns x np) (excitatory 4ms, inhibitory 16ms)
%--------------------------------------------------------------------------
KE   = exp(-P.T)*1000/4;                     % excitatory time constants
KI   = exp(P.TI)*1000/16;                    % inhibitory time constants
KNMDA = exp(P.TN)*1000/100; 

% Voltages
%--------------------------------------------------------------------------
VL   = -70;                                  % reversal  potential leak (K)
VE   =  60;                                  % reversal  potential excite (Na)
VI   = -90;                                  % reversal  potential inhib (Cl)
VR   = -40;                                  % reversal Ca(NMDA)      
VN   = 60; 

CV   = 8/1000;                               % membrane capacitance
GL   = 1;                                    % leak conductance

%%% approx curvature got NMDA entry

 
% mean-field effects: the paramters of the sigmoid activation function
%==========================================================================
if mfm
    
    % covariance among states (mV^2)
    %----------------------------------------------------------------------    
    for i = 1:ns
        for j = 1:np
            Cx{i,j} = x{2}(:,:,i,j);
            Vx(i,j) = Cx{i,j}(1,1);          % population variance

        end
    end

    D   = sparse(diag([1/16 1 1 1]));          % diffusion
    D   = exp(P.S)*D;
    
else
    
    % neural-mass approximation to covariance of states
    %----------------------------------------------------------------------
    try
        Cx = M.Cx;
    catch
        Cx = [75    0.2     0.8 0;
              0.2   0.004  0     0;
              0.8   0      0.02  0;
              0     0      0     0];
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
 
% input
%==========================================================================
if isfield(M,'u')
    
    % endogenous input
    %----------------------------------------------------------------------
    U = u(:)/128;
    
else
    % exogenous input
    %----------------------------------------------------------------------
    U = C*u(:)*exp(P.scale_u);
end

% Exogenous input (to excitatory populations)
%--------------------------------------------------------------------------
try
    B = exp(P.X)/8;
catch
    B = 0;
end
 
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
        E_NMDA = G(i)*GE_NMDA(j,:)*m(i,:)';
        I =      GI(j,:)*m(i,:)';
        
        % extrinsic coupling (excitatory only) and background activity
        %------------------------------------------------------------------
        E =  E + SA(j,:)*a(i,:)' + B;
 
        % Voltage
        %------------------------------------------------------------------
        f{1}(i,j,1) =         (GL*(VL - x{1}(i,j,1)) + ...
                      x{1}(i,j,2)*(VE - x{1}(i,j,1)) + ...
                      x{1}(i,j,3)*(VI - x{1}(i,j,1)) )/CV;
                  
        % Exogenous input (U)
        %------------------------------------------------------------------
        if j == 1
            E = E + U(i);
        end
        
           if (j ==3 || j==2)  %% Pyamidal Cells & interneuron NMDA receptos

           mag_block = 1/(1 + 0.2*exp(-0.062*(exp(P.scale_NMDA))*x{1}(i,j,1))) ;
           f{1}(i,j,1) =  f{1}(i,j,1) + (x{1}(i,j,4)*(VN - x{1}(i,j,1))*mag_block)/CV;
           end
        
        % Conductances
        %------------------------------------------------------------------
        f{1}(i,j,2) = (E - x{1}(i,j,2))*KE(i);
        f{1}(i,j,3) = (I - x{1}(i,j,3))*KI;
        f{1}(i,j,4) = (E_NMDA - x{1}(i,j,4))*KNMDA ;%%and NMDA
     
        
        % 2nd moments - covariances
        %==================================================================
        if mfm
            
            % add curvature-dependent dispersion to flow
            %--------------------------------------------------------------
            %% approx curvature
                      
            if (j ==2 || j ==3)
            x_V = x{1}(i,j,1);
            x_G = x{1}(i,j,4);
            df_dvv = spm_diff('spm_fx_NMDA',x_V,x_G,P,M,[1 1]); 
            df_dvg = spm_diff('spm_fx_NMDA',x_V,x_G,P,M,[1 2]); 
            df_dgv = spm_diff('spm_fx_NMDA',x_V,x_G,P,M,[2 1]); 
      
                fxx  = [df_dvv{1}/CV -1/CV -1/CV  df_dvg{1}/CV;
                -1/CV  0    0   0;
                -1/CV  0    0   0;
                df_dgv{1}/CV   0    0   0];
            
          else
            fxx  = [0 -1/CV -1/CV  0;
                -1/CV  0    0   0;
                -1/CV  0    0   0;
                0  0    0   0];    %%%sparse([2 3 1 1],[1 1 2 3],-1/CV);    % curvature: df(V)/dxx
         end
            f{1}(i,j,1) = f{1}(i,j,1) + tr(Cx{i,j},fxx)/2;

            % df/dx
            %--------------------------------------------------------------
            if (j ==3 || j==2)  %%% NMDA dependant
            
            x_V = x{1}(i,j,1);
            x_G = x{1}(i,j,4);
            df_dv = spm_diff('spm_fx_NMDA',x_V,x_G,P,M,[1]); 
            
            %%% Product and Quotient Rules for NMDA component    
            Sg  = -GL - x{1}(i,j,2) - x{1}(i,j,3) + df_dv;
            fx  = [Sg/CV (VE - x{1}(i,j,1))/CV (VI - x{1}(i,j,1))/CV (mag_block*(VE - x{1}(i,j,1)))/CV;
                   0     -KE(i)              0                 0;
                   0      0                 -KI                0;
                   0      0                  0                 -KNMDA];
            else
            Sg  = -GL - x{1}(i,j,2) - x{1}(i,j,3) + 0;
            fx  = [Sg/CV (VE - x{1}(i,j,1))/CV (VI - x{1}(i,j,1))/CV 0;
                   0     -KE(i)              0                 0;
                   0      0                 -KI                0;
                   0      0                  0                 0];
            end
                
                
            % dCdt
            %--------------------------------------------------------------
           
            St            = fx*Cx{i,j} + D;
            f{2}(:,:,i,j) = St + St';
 

        else
            
            % fixed covariance (Cx)
            %--------------------------------------------------------------
               fxx  = [0 -1/CV -1/CV  0;
                -1/CV  0    0   0;
                -1/CV  0    0   0;
                0  0    0   0];  
            f{1}(i,j,1) = f{1}(i,j,1) + tr(Cx,fxx)/2;
            
        end
        
    end
end
 
% vectorise equations of motion
%==========================================================================
f = spm_vec(f);
 
if nargout < 2, return, end

% Jacobian
%==========================================================================
J = spm_cat(spm_diff('spm_fx_mfm',xin,u,P,M,1));

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
d  = -[2 16].*exp(P.D)/1000;
nk = 3;                                               % number of states
Sp = kron(ones(nk,nk),kron( eye(np,np),eye(ns,ns)));  % states: same pop.
Ss = kron(ones(nk,nk),kron(ones(np,np),eye(ns,ns)));  % states: same source

Dp = ~Ss;                            % states: different sources
Ds = ~Sp & Ss;                       % states: same source different pop.
D  = d(2)*Dp + d(1)*Ds;

% disable for mean field models (temporarily)
%--------------------------------------------------------------------------
if mfm

    Cp = kron(ones(nk,1),kron(kron(eye(np,np) ,eye(ns,ns)),ones(1,nk*nk)));
    Cs = kron(ones(nk,1),kron(kron(ones(np,np),eye(ns,ns)),ones(1,nk*nk)));
    Dp = kron(kron( eye(np,np),eye(ns,ns)),kron(ones(nk,nk),ones(nk,nk)));
    Ds = kron(kron(ones(np,np),eye(ns,ns)),kron(ones(nk,nk),ones(nk,nk)));

    Sp = spm_cat({Sp Cp; Cp' Dp});
    Ss = spm_cat({Ss Cs; Cs' Ds});

    Dp = ~Ss;                        % states: different sources
    Ds = ~Sp & Ss;                   % states: same source different pop.

    D  = d(2)*Dp + d(1)*Ds;
end


% Implement: dx(t)/dt = f(x(t - d)) = inv(1 - D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%--------------------------------------------------------------------------
Q  = inv(speye(length(J)) - D.*J);


% trace(a*b)
%--------------------------------------------------------------------------
function x = tr(a,b);
%__________________________________________________________________________

b   = b';
x   = a(:)'*b(:);
