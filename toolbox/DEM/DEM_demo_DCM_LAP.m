% Demo applying the Laplace scheme to DCM with hidden states
%__________________________________________________________________________
% This routine demonstrates Generalized filtering for a DCM (Dynamic Causal
% Model) of fMRI responses using simulated data. This is an endogenous 
% DCM in that there are no exogenous inputs. The demonstration specifies 
% and inverts a full connectivity model and then illustrates post-hoc model
% optimization to recover (discover) the true architecture. It concludes 
% with an automatic model optimization in terms of the prior variances over
% coupling parameters.
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_DCM_LAP.m 4339 2011-06-02 13:04:49Z karl $
 
% Specify a DCM to generate synthetic data
%==========================================================================
clear
 
% DEM Structure: create inputs
% -------------------------------------------------------------------------
T  = 256;
TR = 3.22;
t  = (1:T)*TR;
n  = 4;
U  = spm_conv(randn(n,T),0,2)/exp(4/2);
 
% set inversion parameters
% -------------------------------------------------------------------------
M(1).E.s  = 1/2;         % smoothness of random fluctuations
M(1).E.d  = 2;           % embedding dimension
M(1).E.n  = 6;           % embedding dimension
M(1).E.nE = 32;          % maximum number of DEM iterations
 
 
% priors
% -------------------------------------------------------------------------
A  = ones(n,n);
B  = zeros(n,n,0);
C  = zeros(n,n);
D  = zeros(n,n,0);
 
[pE,pC] = spm_dcm_fmri_priors(A,B,C,D);
 
 
% adjust M.f (GF works in time bins not seconds) and initialize M.P
% -------------------------------------------------------------------------
M(1).f  = inline(['spm_fx_fmri(x,v,P)*' num2str(TR)],'x','v','P');
M(1).g  = 'spm_gx_fmri';
M(1).x  = sparse(n,5);
M(1).pE = pE;
M(1).pC = pC;
 
M(2).v  = sparse(n,1);
 
 
% simulate endogenous dynamics
%==========================================================================
 
% true parameters (stochastic estimates)
% -------------------------------------------------------------------------
pP   = pE;
a    =  0.3;
b    = -0.3;
c    = -0.5;
pP.A = [c  a  0  0  0;
        a  c  b  0  0;
        0  b  c  a  0;
        0  0  a  c  b;
        0  0  0  a  c];
     
pP.A = pP.A(1:n,1:n);
pP.C = eye(n,n);
SIM  = spm_DEM_generate(M,U,{pP},{6,16},{16});
 
 
% Show simulated response
%--------------------------------------------------------------------------
spm_figure('Getwin','Figure 1');
spm_DEM_qU(SIM.pU)
 
 
% Specify generative model for inversion (DCM)
% =========================================================================
 
% set inversion parameters
% -------------------------------------------------------------------------
DCM.M       = M;
DCM.M(2).v  = 0;
 
% Specify hyper-priors on precisions
% -------------------------------------------------------------------------
W           = exp(spm_vec(sparse(1:n,1,(6 - 16),n,5) + 16));
DCM.M(1).xP = exp(6);
DCM.M(1).V  = exp(6);        % prior log precision (noise)
DCM.M(1).W  = diag(W);       % fixed precision (hidden-state)
DCM.M(2).V  = exp(16);       % fixed precision (hidden-cause)
 
% Add data
% -------------------------------------------------------------------------
DCM.Y       = SIM.Y;
 
% Full connectivity inversion
% =========================================================================
F  = ones(n,n);
B  = zeros(n,n,0);
C  = zeros(n,1);
D  = zeros(n,n,0);
 
[pE pC]     = spm_dcm_fmri_priors(F,B,C,D);
DCM.M(1).pE = pE;
DCM.M(1).pC = pC;
FULL        = spm_LAP(DCM);
 
 
% Search model space with Savage-Dickey density ratio
% =========================================================================
[A K Nk]  = spm_dcm_sparse_priors(n);
 
% find candidate models based on full-connectivity
% -------------------------------------------------------------------------
pE    = FULL.M(1).pE.A;
qE    = FULL.qP.P{1}.A;
qC    = FULL.qP.C(1:n*n,1:n*n);
pC    = FULL.M(1).pC(1:n*n,1:n*n);
for i = 1:length(A)
    k       = find(~A{i});       
    rE      = pE;
    rC      = pC;
    rE(k)   = 0;
    rC(k,k) = 0;
    P(i,1)  = spm_log_evidence(qE,qC,pE,pC,rE,rC);
end
 
% find true model
% -------------------------------------------------------------------------
for i = 1:length(A);
    if ~any(spm_vec(~~SIM.pP.P{1}.A - A{i}))
        tA = i; break
    end
end
 
% log-posterior (model)
% -------------------------------------------------------------------------
PP    = exp(P - max(P));
PP    = PP/sum(PP);
 
% Graphics (density on parameter and model space)
% -------------------------------------------------------------------------
spm_figure('Getwin','Figure 2'); clf
 
subplot(2,2,1)
spm_plot_ci(spm_vec(qE),qC),     hold on
bar(spm_vec(pP.A),1/2), hold off
title('true and MAP connections','FontSize',16)
axis square
 
subplot(2,2,2)
bar(P)
title('log-posterior','FontSize',16)
xlabel('model','FontSize',12)
ylabel('log-probability','FontSize',12)
axis square
 
subplot(2,2,3)
plot(Nk,    P,    '.k','MarkerSize',16), hold on
plot(Nk(tA),P(tA),'.r','MarkerSize',32), hold off
title('log-evidence','FontSize',16)
xlabel('graph size','FontSize',12)
ylabel('log-probability','FontSize',12)
axis square
 
subplot(2,2,4)
bar(PP)
title('posterior','FontSize',16)
xlabel('model','FontSize',12)
ylabel('probability','FontSize',12)
axis square
 
 
% Compare true and AMS adjacency
% -------------------------------------------------------------------------
spm_figure('Getwin','Figure 3'); clf
 
% Automatic model selection (optimizing prior variances on parameters)
% =========================================================================
[rE,rC]   = spm_dcm_optimise(qE,qC,pE,pC);
rA        = spm_unvec(diag(rC),pE);
 
subplot(2,2,1)
imagesc(full(A{tA} - diag(diag(A{tA}))))
title('true adjacency','FontSize',16)
xlabel('source','FontSize',12)
ylabel('target','FontSize',12)
axis square
 
subplot(2,2,2)
imagesc(full(rA))
title('optmised priors','FontSize',16)
xlabel('source','FontSize',12)
ylabel('target','FontSize',12)
axis square
