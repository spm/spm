function DEM_demo_fMRI_PEB
% Demonstration of PEB for multisession spectral DCM studies
%__________________________________________________________________________
% This demonstration routine illustrates the analysis of a multisession
% fMRI study using spectral DCM. Crucially, the between session effects are
% characterised using empirical Bayes and Bayesian model reduction. This
% means that the original session data are only inverted once (at the
% within session level). The resulting posterior estimates and then used to
% make inferences about between session effects (e.g., time or drug
% effects). The basic question addressed in this sort of analysis is where
% between session effects are expressed in terms of connectivity or
% parameters of neuronal fluctuations. These sorts of effects are specified
% in a second level design matrix in the usual way and can be identified
% using Bayesian model reduction.
%
% in this example, we analyse three sessions with a monotonic change in the
% intrinsic (self) connectivity over three sessions. This involves
% decreases in diagonal A parameters at the first two levels of a simple
% three node hierarchy - and an increase at the highest (third) level.
% Physiologically, this corresponds to a decrease in self-inhibition (or
% increase in excitability) in the lower notes for regions, as time goes
% on.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_fMRI_PEB.m 7679 2019-10-24 15:54:07Z spm $
 
 
 
% Simulate timeseries
%==========================================================================
rng('default')
 
% assume we have data from eight subjects with a TR of two seconds
% -------------------------------------------------------------------------
T  = 256*8;                           % number of observations (scans)
TR = 2;                               % repetition time or timing
t  = (1:T)*TR;                        % observation times
n  = 3;                               % number of regions or nodes
 
% setup model
% -------------------------------------------------------------------------
options.nonlinear  = 0;
options.two_state  = 0;
options.stochastic = 0;
options.centre     = 1;
options.induced    = 1;
 
A   = ones(n,n);
B   = zeros(n,n,0);
C   = zeros(n,n);
D   = zeros(n,n,0);
pP  = spm_dcm_fmri_priors(A,B,C,D,options);
 
 
% true parameters - a simple hierarchy of three nodes
% -------------------------------------------------------------------------
pP.A = [  0  -.2    0;
         .2    0  -.1
          0   .3    0];
pP.C = eye(n,n);
pP.transit = randn(n,1)/16;
 
% spectral density of neuronal fluctuations and observation noise (eg)
% -------------------------------------------------------------------------
[Gu,Gn,Hz,dt] = spm_csd_fmri_gu(pP,TR);
Gu            = Gu(:,1,1)*ones(1,n);
Gn            = Gn(:,1,1)*ones(1,n);
 
 
% session specific parameters: here, intrinsic connectivity
% -------------------------------------------------------------------------
X    = exp(-(0:2))';                    % between session decay
B    = [  1    0    0;                  % of intrinsic connections
          0    1    0
          0    0   -1];
X    = X - mean(X);                     % detrended session-effects
 
 
% simulate response to endogenous fluctuations
%==========================================================================
for s = 1:numel(X)
    
    % parameters in this session
    % ---------------------------------------------------------------------
    P    = pP;
    P.A  = P.A + B*X(s);
    
    % integrate states with endogenous fluctuations (u)
    % ---------------------------------------------------------------------
    M.f  = 'spm_fx_fmri';
    M.x  = sparse(n,5);
    U.u  = spm_rand_power_law(Gu,Hz,dt,T);
    U.dt = TR;
    x    = spm_int_J(P,M,U);
    
    % haemodynamic observer
    % ---------------------------------------------------------------------
    for i = 1:T
        y(i,:) = spm_gx_fmri(spm_unvec(x(i,:),M.x),[],P)';
    end
    
    % response with observation noise (e)
    % ---------------------------------------------------------------------
    e       = spm_rand_power_law(Gn,Hz,dt,T);
    Y(s).y  = y + e;
    Y(s).dt = TR;
    PP(s)   = P;
    
end
 
 
% show simulated response (last session)
%--------------------------------------------------------------------------
i = 1:128;
spm_figure('Getwin','Figure 1'); clf
subplot(2,2,1)
plot(t(i),U.u(i,:))
title('Endogenous fluctuations','FontSize',16)
xlabel('Time (seconds)')
ylabel('Amplitude')
axis square
 
subplot(2,2,2), hold off
plot(t(i),x(i,n + 1:end),'c'), hold on
plot(t(i),x(i,1:n)), hold off
title('Hidden states','FontSize',16)
xlabel('Time (seconds)')
ylabel('Amplitude')
axis square
 
subplot(2,2,3)
plot(t(i),y(i,:),t(i),e(i,:),':')
title('Hemodynamic response and noise','FontSize',16)
xlabel('Time (seconds)')
ylabel('Amplitude')
axis square
 
 
% nonlinear system identification (DCM for CSD)
%==========================================================================
dcm.options = options;
 
dcm.a    = logical(pP.A);
dcm.b    = zeros(n,n,1);
dcm.c    = zeros(n,1);
dcm.d    = zeros(n,n,0);
 
% add responsse to each session
% -------------------------------------------------------------------------
for s = 1:numel(X)
    DCM{s,1}   = dcm;
    DCM{s,1}.Y = Y(s);
end
 
% first level inversion - spectral DCM
% =========================================================================
CSD = spm_dcm_peb_fit(DCM);
 
% show estimates for a single session
% -------------------------------------------------------------------------
spm_figure('Getwin','Figure 1');
 
j  = 1:n^2;
d  = 3;
subplot(2,2,4); hold off
spm_plot_ci(CSD{d}.Ep.A(:),CSD{d}.Cp(j,j)), hold on
bar(PP(d).A(:),1/4), hold off
title('True and MAP connections (Deterministic)','FontSize',16)
axis square
 
 
% Parametric empirical Bayes
%==========================================================================
% having inverted every session, we can construct a between session model
% at the second level. We are interested in identifying which parameters
% change according to the session specific effects encoded in the design
% matrix. This design matrix is supplemented with a constant term, creating
% two times the number of parameters at the second level. Using Bayesian
% model reduction, we can then examine models with and without each
% parameter for the  constant (first) and hypothesised (second) explanatory
% variables in the design matrix. Crucially, we want to explain as much  as
% possible using just between session differences. This can be implemented
% by setting the between session covariance of random effects to a
% relatively small value; here the original prior variance divided by 32.
%--------------------------------------------------------------------------
clear M;
beta  = 32;
 
M.X   = [X.^0 X];                  % between session explanatory variables
M.hE  = 0;                         % prior expectation of log precision
M.hC  = 1/16;                      % prior covariance of precision
M.bE  = CSD{1}.M.pE;               % prior expectations over sessions
M.bC  = CSD{1}.M.pC;               % prior covariance over sessions
M.pC  = CSD{1}.M.pC/beta;          % prior covariance between sessions
 
field = {'A','a'};                 % parameters of interest
PEB   = spm_dcm_peb(CSD,M,field);  % empirical Bayesian inversion
BMA   = spm_dcm_peb_bmc(PEB);      % Bayesian model reduction and averaging
 
% overlay true between session effects over BMA density over parameters
%--------------------------------------------------------------------------
j     = spm_find_pC(CSD{1},'A');
p     = spm_vec(B); p = p(j);
subplot(3,2,4), hold on, bar(p,1/4), hold off
 
return
 