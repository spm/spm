function [x,DEM] = DEM_demo_MMN_gen(P,G,U)
% generates hidden states for a MMN roving paradigm using spm_DEM
% FORMAT [x,DEM] = DEM_demo_MMN_gen(P,G,U);
%
% P   - parameters
% G   - generative (response) model
% U   - design
%
% x   - hidden neuronal states {(ns x (nr x 8)}; ... }
% DEM - stucture array of DEM structures
%
% see DEM_demo_MMN
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_MMN_gen.m 4146 2010-12-23 21:01:39Z karl $
 
 
 
 
% parameters (recognition) model
%==========================================================================
 
% model parameters; Jacobian and mixing parameters
%--------------------------------------------------------------------------
F       = [-2  4 ;
           -4 -1]/32;
A.g     = [ 0  1 ;                           % amplitude
            4  0];                           % frequency
A.f     = exp(P.f).*F;                       % dynamics
 
 
% stimulus parameters
%--------------------------------------------------------------------------
B.g     = [ 0       1 ;
          (P.g + 1) 0];
B.f     = A.f;
 
      
 
% model specification (which parameters can be learned)
%==========================================================================
ip      = [2];        
pC      = sparse(ip,ip,1,8,8);                % prior covariance  
 
% level 1
%--------------------------------------------------------------------------
M(1).m  = 1;                                 % 1 input or cause
M(1).n  = 2;                                 % 2 hidden states
M(1).l  = 2;                                 % 2 outputs (amplitude and Hz)
 
M(1).f  = inline('P.f*x + [v; 0]','x','v','P');
M(1).g  = inline('P.g*x + [0; 16]','x','v','P');
M(1).pE = A;                                 % The prior expectation
M(1).pC = pC;                                % The prior covariance
M(1).Q  = {eye(2)};                          % error precision (data)
M(1).hE = 5;                                 % error log-precision prior
M(1).hC = 1/32;                              % error log-precision prior
M(1).W  = exp(16);                           % error precision (dynamics)
 
% level 2
%--------------------------------------------------------------------------
M(2).l  = 1;                                 % 1 output 
M(2).V  = exp(16);                           % error precision (cause)
 
 
% number of steps
%--------------------------------------------------------------------------
M(1).E.nD = 1;
M(1).E.nE = 1;
M(1).E.nM = 8;
 
% parameterise input or cause
%==========================================================================
n       = size(U.X,1);                       % number of trials
ns      = G.ns;                              % samples per trial
t       = [1:ns]*U.dt;                       % time bin (sec)
c       = spm_erp_u(t,P,G)'/32;              % Gaussian cause
 
% DEM estimation
%==========================================================================
 
% Cycle through presentations
%--------------------------------------------------------------------------
DEM   = {};
randn('state',0);
for i = 1:n
 
    % Stimulus
    %----------------------------------------------------------------------
    DEM{i}   = spm_DEM_generate(M,c,B,{16 16},{16 []});
    DEM{i}.U = c;
    
    % Invert and pass posteriors to priors of M
    %----------------------------------------------------------------------
    DEM{i}   = spm_DEM(DEM{i});           % compute conditional densities
    M(1).pE  = DEM{i}.qP.P{1};            % update parameter estimates
    M(1).hE  = DEM{i}.qH.h{1};            % update hyperparameter estimates
 
end
 
% extract errors (cf. pyramidal activity)
%--------------------------------------------------------------------------
for i = 1:n
 
    % prediction error (precision weighted)
    %----------------------------------------------------------------------
    R      = spm_DEM_EEG(DEM{i},U.dt,[1 2]);
    
    % save trial (causal and hidden state error for each level)
    %----------------------------------------------------------------------
    x{i,1} = spm_cat(R(:))';
 
end
