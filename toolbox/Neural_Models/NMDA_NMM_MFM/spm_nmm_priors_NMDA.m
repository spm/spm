function [pE,pC] = spm_nmm_priors_NMDA(A,B,C)
% prior moments for a neural-mass model of ERPs
% FORMAT [pE,pC] = spm_nmm_priors_NMDA(A,B,C)
%
% A{3},B{m},C  - binary constraints on extrinsic connections
%
% pE - prior expectation - f(x,u,P,M)
%
% population variance
%--------------------------------------------------------------------------
%     E.S    - variance
%
% synaptic parameters
%--------------------------------------------------------------------------
%    pE.T    - synaptic time constants
%    pE.G    - intrinsic connectivity
%
% connectivity parameters
%--------------------------------------------------------------------------
%    pE.A    - extrinsic
%    pE.B    - trial-dependent
%    pE.C    - stimulus input
%
%    pE.SA - switches on extrinsic (excitatory)
%    pE.GE - switches on intrinsic (excitatory)
%    pE.GI - switches on intrinsic (inhibitory)
%
% stimulus and noise parameters
%--------------------------------------------------------------------------
%    pE.R    - onset and dispersion
%    pE.D    - delays
%    pE.X    - exogenous background activity
%
% pC - prior covariances: cov(spm_vec(pE))
%
% Because priors are specified under log normal assumptions, most
% parameters are simply scaling coefficients with a prior expectation
% and variance of one.  After log transform this renders pE = 0 and
% pC = 1;  The prior expectations of what they scale are specified in
% spm_erp_fx
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_nmm_priors_NMDA.m 4820 2012-08-01 12:20:00Z guillaume $
 
 
% disable log zero warning
%--------------------------------------------------------------------------
warning('off','MATLAB:log:logOfZero');
n     = size(C,1);                                % number of sources
u     = size(C,2);                                % number of inputs
 
% parameters for neural-mass forward model
%==========================================================================
 
% population variance
%--------------------------------------------------------------------------
pE.S   = 0;                 pC.S = 1/64;
pE.CV   = 0;                 pC.CV = 1/64;

% set intrinic [excitatory] time constants
%--------------------------------------------------------------------------
pE.T    = log(ones(n,1));    pC.T = ones(n,1)/4;    % excitatory constants
pE.TI   = log(ones(n,1));    pC.TI = ones(n,1)/4;    % inhibitory constants
pE.TN   = log(ones(n,1));    pC.TN = ones(n,1)/4;    % modulatory constants

pE.G   = log(ones(n,1));    pC.G = ones(n,1)/4;    % intrinsic connections



pE.scale_u = log(0.25) ;  pC.scale_u = 1/4; %%% %% MFM 2.25 %%was u= 0.25 with NMDA = 4 generates high theta  ~8Hz

pE.scale_NMDA =   log(4) ;  pC.scale_NMDA = 1/4;
pE.scale_geP1 =   log(1)  ;  pC.scale_geP1 = 1/4;
pE.scale_geP2 =   log(1)  ;  pC.scale_geP2 = 1/4;
pE.scale_geS =    log(1)  ;  pC.scale_geS = 1/4;
pE.scale_geS2 =   log(1)  ;  pC.scale_geS2 = 1/4;
pE.scale_gi =     log(1)  ;  pC.scale_gi = 1/4;

%%% ------- BELOW IS DOPAMINE SPECIFIC for 2 Trials

%% Dopamine should 
%% decrease AMPA
%% Increase GABA
%% Increase NMDA (i.e.alpha magnesiun switch NMDA)

pE.B_gi{1}   = 0 ; pC.B_gi{1}   =  1/4;
pE.B_geP1{1} = 0;  pC.B_geP1{1} =  1/4;
pE.B_geP2{1} = 0;  pC.B_geP2{1} =  1/4;
pE.B_geS{1}  = 0;  pC.B_geS{1}  =  1/4;
pE.B_geS2{1} = 0;  pC.B_geS2{1} =  1/4;
pE.B_NMDA{1} = 0;  pC.B_NMDA{1} =  1/4;
pE.B_u{1}    = 0;  pC.B_u{1}    =  1/4


% set extrinsic connectivity
%--------------------------------------------------------------------------
Q     = sparse(n,n);
for i = 1:length(A)
    A{i} = ~~A{i};
    pE.A{i} = A{i}*32 - 32;                       % forward
    pC.A{i} = A{i};                               % backward
    Q      = Q | A{i};                            % and lateral connections
end
 
for i = 1:length(B)
    B{i} = ~~B{i};
    pE.B{i} = 0*B{i};                             % input-dependent scaling
    pC.B{i} = B{i};
    Q      = Q | B{i};
end
C     = ~~C;
pE.C   = C*33 - 32;                               % where inputs enter
pC.C   = C;
 

% connectivity switches
%==========================================================================

% extrinsic connections
%--------------------------------------------------------------------------
pE.SA   = [1   0   1;  
           0   1   1;
           0   0   0];   pC.SA = sparse(3,3);
 
% intrinsic connections (np x np) - excitatory
%--------------------------------------------------------------------------
pE.GE   = [0   0   1/2;  
           0   0   1;
           1/2 0   0];   pC.GE = sparse(3,3);
 
% intrinsic connections (np x np) - inhibitory
%--------------------------------------------------------------------------
pE.GI   = [0   1/4 0;
           0   0   0;
           0   1   0];   pC.GI = sparse(3,3);


% set stimulus parameters: onset and dispersion
%--------------------------------------------------------------------------
pE.R   = sparse(u,2);    pC.R  = ones(u,1)*[1/16 1/16];

% and delays (intrinsic and extrinsic)
%--------------------------------------------------------------------------
pE.D   = [0 0];          pC.D  = [1 1]/128;

% Exogenous background activity
%--------------------------------------------------------------------------
pE.X   = 0;              pC.X  = 1/128;

warning('on','MATLAB:log:logOfZero');
 

% covariance of neural parameters
%==========================================================================
pC    = diag(sparse(spm_vec(pC)));
