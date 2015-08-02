function DEM_demo_ontology
% This demonstration routine illustrates how a generative model can be used
% to furnish a computational nosology. In brief, it generates symptoms and
% diagnostic profiles from hidden or latent exogenous causes (e.g.,
% therapeutic interventions) that are mediated by latent (pathophysiological
% and psychopathological) states.  Pathophysiological trajectories  are
% modelled with a Lorenz attractor that (with a linear mapping)
% produces (two-dimensional) psychopathology. In turn, the
% psychopathological states generate symptoms (with a non-linear function
% of linear mixtures) and diagnostic outcomes (with a softmax function of
% diagnostic potential). The psychopathological state of a subject is
% associated with a diagnostic potential in terms of its Euclidean distance
% from disease categories (locations in the associated state space).
%
% We start by simulating a relapsing-remitting disease process and then
% infer the latent states and parameters of a particular subject.
% This is then repeated in the setting of a therapeutic intervention.
% The demonstration then briefly considers model identification and
% selection by focusing on the mapping between pathophysiology and
% psychopathology. Finally, We consider, prognosis and prediction by
% estimating subject-specific parameters prior to therapy and then
% predicting putative response in the future, based upon a posterior
% predictive density.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_ontology.m 6511 2015-08-02 15:05:41Z karl $
 
 
% Set up the generative model
%==========================================================================
rng('default')
 
% length of trajectory
%--------------------------------------------------------------------------
N         = 64;                               % number of assessments
 
% serial correlations and other inversion parameters
%--------------------------------------------------------------------------
M(1).E.s  = 1/2;
M(1).E.n  = 2;
M(1).E.K  = 1/16;
M(1).E.nE = 32;
 
% therapeutic intervention; ranging between zero and one, starting halfway
% through the patient assessment
%--------------------------------------------------------------------------
U        = spm_phi(((1:N) - N/2));
 
 
% level 1: the level that generates diagnostic and symptom profiles g(v) from
% latent (psychopathological) causes (v)
%--------------------------------------------------------------------------
P.A      = randn(6,2)/32;
P.B      = [0 0 1 1; 0 1 0 1];
 
M(1).g   = @(x,v,P)[spm_softmax(-sum((P.B - v*ones(1,4)).^2)'); tanh(P.A*exp(v))];
M(1).pE  = P;
M(1).V   = exp(12);
 
 
% level 2: the level that generates latent causes v = g(x) from (pathophysiological)
% states (x) that are subject to interventions (U)
%--------------------------------------------------------------------------
P.A      = [10 24 1];
P.B      = [2 0;1/2 1];
 
M(2).f   = @(x,v,P)[-P.A(1) P.A(1) 0; ((1 - v*P.A(3))*P.A(2) - x(3)) -1 0; x(2) 0 -8/3]*x/32;
M(2).g   = @(x,v,P) P.B*x([2 3])/16;
M(2).x   = [2; 4; 32];
M(2).pE  = P;
M(2).pC  = diag([1 1 1 0 0 0 0]);
M(2).V   = exp(8);
M(2).W   = exp(4);
 
% level 3: the exogenous (therapeutic) interventions
%--------------------------------------------------------------------------
M(3).v   = U(:,1);
M(3).V   = exp(16);
 
 
% illustrate trajectories with and without therapy
%==========================================================================
 
% set subject-specific parameters: parameters of the attractor in P{2}.A
%--------------------------------------------------------------------------
P         = {M.pE};
P{2}.A(2) = 32;
 
% natural progression without therapy: DEM.U = U*0
%--------------------------------------------------------------------------
DEM.U  = U*0;
DEM    = spm_DEM_generate(M,DEM.U,P);
DEM    = spm_DEM(DEM);
 
spm_figure('GetWin','Figure 1'); clf
spm_DEM_plot(DEM)
 
% repeat with therapy: DEM.U = U
%--------------------------------------------------------------------------
DEM.U  = U;
DEM    = spm_DEM_generate(M,DEM.U,P);
DEM    = spm_DEM(DEM);
 
% plot true and inferred trajectories
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_DEM_plot(DEM)
 
 
% prognosis and prediction
%==========================================================================
 
% Estimate subject-specific parameters at presentation (1:N/2 assessments)
%--------------------------------------------------------------------------
PEM         = DEM;
PEM.U       = PEM.U(:,1:N/2);
PEM.Y       = PEM.Y(:,1:N/2);
PEM         = spm_DEM(PEM);
 
% plot true and inferred parameters
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
spm_DEM_qP(PEM.qP,PEM.pP)
 
% update priors and predict with zero precision data
%--------------------------------------------------------------------------
PEM         = spm_ADEM_update(PEM,0);
PEM.M(1).V  = 0;
PEM.M(2).W  = exp(8);
PEM.M(2).pC = 0;
 
PEM.U       = U(N/4:end);
PEM.Y       = zeros(size(PEM.Y,1),size(PEM.U,2));
PEM         = spm_DEM(PEM);
 
% plot predicted trajectories
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4'); clf
spm_DEM_plot(PEM)
 
% now repeat without therapeutic intervention
%==========================================================================
PEM.U = U(N/4:end)*0;
PEM   = spm_DEM(PEM);
 
spm_figure('GetWin','Figure 5'); clf
spm_DEM_plot(PEM)
 
% overlay trajectories with increasing levels of therapy (u)
%--------------------------------------------------------------------------
u     = linspace(0,2,8);
for i = 1:length(u)
    
    % predict a response to treatment level
    %----------------------------------------------------------------------
    PEM.U   = U(N/4:end)*u(i);
    PEM     = spm_DEM(PEM);
    
    % plot response trajectory
    %----------------------------------------------------------------------
    spm_figure('GetWin','Figure 5'); subplot(3,2,6), hold on
    plot(PEM.qU.v{2}(1,:),PEM.qU.v{2}(2,:),'r:')
    
    % record predictive efficacy of treatment
    %----------------------------------------------------------------------
    R(i)    = PEM.qU.v{1}(1,end);
    
end
 
%  show implicit dose-response relationships
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 5'); subplot(3,2,5), hold off
bar(u,R)
title('dose response curve','fontsize',16)
xlabel('level of therapeutic intervention')
xlabel('probability of remission')
 
% return
 
 
% model identification and selection
%==========================================================================
% in this section, we assume a cohort of subjects with known
% pathophysiology have been identified and estimate the generic parameters
% mapping pathophysiology to psychopathology.  These estimates are then
% entered into any empirical Bayesian analysis to perform Bayesian model
% comparison to see which of these parameters are necessary. By
% construction, the parameter linking the first pathological state to the
% second psychopathology is redundant.
%--------------------------------------------------------------------------
n      = 8;                                        % number of subjects
for i  = 1:n
    DCM{i,1}           = spm_DEM_generate(M,U);
    DCM{i,1}.U         = U;
    DCM{i,1}.M(2).pE.B = eye(2);
    DCM{i,1}.M(2).pC   = diag([0 0 0 1 1 1 1]);
end
 
%  Bayesian model inversion and hierarchical (empirical Bayesian) modelling
%--------------------------------------------------------------------------
DCM    = spm_dcm_fit(DCM);
PEB    = spm_dcm_peb(DCM);
 
% Bayesian model comparison using Bayesian model reduction
%--------------------------------------------------------------------------
spm_dcm_bmr_all(PEB)
subplot(3,2,3), hold on, bar(1:4,M(2).pE.B(:),1/3), hold off
 
 
 
 
% plotting sub function
%==========================================================================
function spm_DEM_plot(DEM)
 
% plot hidden states and causes
%--------------------------------------------------------------------------
if isfield(DEM,'pU')
    spm_DEM_qU(DEM.qU,DEM.pU)
else
    spm_DEM_qU(DEM.qU)
end
 
% supplement with trajectories
%--------------------------------------------------------------------------
subplot(6,2,2), imagesc(DEM.qU.v{1}(1:4,:))
title('diagnostic & symptom profile  ','fontsize',16)
subplot(6,2,4), imagesc(DEM.qU.v{1}(5:end,:))
 
% supplement with trajectories
%--------------------------------------------------------------------------
a   = [-1 3 -1 3];
vi  = linspace(a(1),a(2),64);
vj  = linspace(a(3),a(4),64);
for i = 1:length(vi)
    for j = 1:length(vj)
        x      = DEM.M(1).x;
        v      = [vi(i); vj(j)];
        p      = DEM.M(1).g(x,v,DEM.qP.P{1});
        p      = p(1:4);
        s(j,i) = p'*log(p);
        [p,q]  = max(p);
        d(j,i) = q;
    end
end
d    = d.*(min(s(:)) - s);
 
subplot(3,2,6), imagesc(a(1:2),a(3:4),d), axis xy, hold on
plot(DEM.qU.v{2}(1,:),DEM.qU.v{2}(2,:),'r'), hold off
title('latent psychopathology','fontsize',16)
 
 
subplot(3,2,4),title('latent pathophysiology','fontsize',16)
subplot(3,2,3),title('latent psychopathology','fontsize',16)
subplot(3,2,1),title('predicted symptoms and error','fontsize',16)
subplot(3,2,5),title('pathology and therapy','fontsize',16)
set(gca,'YLim',[-.2 1.2])
