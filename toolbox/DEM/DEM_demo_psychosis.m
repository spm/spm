function DEM_demo_psychosis
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
 
 
%% Set up the generative model
%==========================================================================
rng('default')
 
% length of trajectory
%--------------------------------------------------------------------------
T         = 128;                               % number of assessments

% dynamics and parameters of a Lorentz system (with Jacobian)
%==========================================================================
% dxdt = f(x)
% y    = h(g(x)) + e
%--------------------------------------------------------------------------

% level 1: the level that generates symptom profiles g(v) from
% latent (psychopathological) causes (v)
%--------------------------------------------------------------------------
M.h  = @(x,P)spm_phi([x(:,1) - P.u(1), x(:,2) - P.u(2)]*P.s);
 
% level 2: the level that generates latent causes v = g(x) from (pathophysiological)
% states (x) that are subject to interventions (U)
%--------------------------------------------------------------------------
M.f  = @(x,v,P,M)[-P.A(1) P.A(1) 0; ((1 - v*P.A(3))*P.A(2) - x(3)) -1 0; x(2) 0 -8/3]*x/P.t;
M.d  = @(x,P) x(:,[1 3])*P.B/16;

% level 3: the level that generates latent causes v = g(x) from (pathophysiological)
% states (x) that are subject to interventions (U)
%--------------------------------------------------------------------------
M.v  = @(T,P) spm_dctmtx(T,numel(P.v))*P.v(:);

% y    = h(g(x)) + e
%--------------------------------------------------------------------------
pE.A  = [10 16 1];
pE.B  = [2 0;0 1];
pE.C  = 8;
pE.x  = [2; 4; 32];
pE.u  = [1,1,1.5];
pE.v  = [0 0];
pE.s  = 16;
pE.t  = 64;
pC    = diag([1 1 1  0 0 0 0  1/64  0 0 0  0 0 0 1/32 1/32  1 1]);

% model specification
%==========================================================================
M.Nmax = 32;                   % maximum number of iterations
M.G    = @spm_psychosis_gen;   % generative function
M.FS   = @(Y)spm_conv(Y,16,0); % feature selection  (link function)
M.pE   = pE;                   % prior expectations (parameters)
M.pC   = pC;                   % prior covariances  (parameters)
M.hE   = 4;                    % prior expectation  (log-precision)
M.hC   = 1/512;                % prior covariances  (log-precision)
M.T    = T;                    % number of time points
U      = zeros(T,1);

% model inversion with Variational Laplace (Gauss Newton)
%==========================================================================

% generate synthetic symptom scores with a particular Rayleigh parameter
%--------------------------------------------------------------------------
P      = pE;
P.A(2) = 21;
P.v    = randn(size(pE.v))/4;
Y      = spm_psychosis_gen(P,M,U);

% multi-start inversion with Variational Laplace
%--------------------------------------------------------------------------
pA    = 8:4:24;                % initial values of Rayleigh parameter
F     = -Inf;
for i = 1:numel(pA)
    
    % initialise parameters (i.e., prior expectations) and invert
    %----------------------------------------------------------------------
    M.pE.A(2)    = pA(i);
    [ep,cp,eh,f] = spm_nlsi_GN(M,U,Y);
    
    % retain posteriors if model evidence has increased
    %----------------------------------------------------------------------
    if f > F
        F      = f;
        DCM.Ep = ep;
        DCM.Cp = cp;
        DCM.Eh = eh;
        DCM.F  = F;
    end
end

% illustrate results
%--------------------------------------------------------------------------
DCM.M  = M;
DCM.Y  = Y;
DCM.U  = U;
DCM.P  = P;

spm_DCM_plot(DCM)



%% repeat for a group subjects with random parametric effects
%--------------------------------------------------------------------------
n     = 16;                           % number of subjects in each group
pP    = zeros(n,1);                   % value of Rayleigh parameter
qP    = zeros(n,1);                   % posterior estimates
pY    = [];                           % true symptom scores
qY    = [];                           % predicted scores
pA    = 8:2:24;                       % initial of Rayleigh parameters
for j = 1:n
    
    % generate synthetic symptom scores with a particular Rayleigh parameter
    %----------------------------------------------------------------------
    P      = pE;
    P.A(2) = rand*16 + 8;
    P.v    = randn(size(pE.v))/4;
    Y      = spm_psychosis_gen(P,M,U);
    
    % model inversion with Variational Laplace
    %----------------------------------------------------------------------
    F     = -Inf;
    for i = 1:numel(pA)
        
        % initialise parameters (i.e., prior expectations) and invert
        %------------------------------------------------------------------
        M.pE.A(2)    = pA(i);
        [ep,cp,eh,f] = spm_nlsi_GN(M,U,Y);
        
        % retain posteriors if model evidence has increased
        %------------------------------------------------------------------
        if f > F
            F  = f;
            Ep = ep;
        end
    end

    y     = spm_psychosis_gen(Ep,M,U);
    pY    = [pY; Y];
    qY    = [qY; y];
    pP(j) =  P.A(2)
    qP(j) = Ep.A(2);
    
end


% report true and inferred symptom scores and parameters
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
subplot(3,2,2)
plot(pP,qP,'o')

 

function [y,s,x] = spm_psychosis_gen(P,M,U)
% FORMAT [y,s,x] = spm_psychosis_gen(P,M,U)

% generate symptom scores from underlying physiological dynamics
%--------------------------------------------------------------------------
T    = size(U,1);
M.x  = P.x;                          % initial physiological states
U    = M.v(T,P);                     % fluctuations
x    = spm_int_J(P,M,U);             % integrate trajectory over time
s    = M.d(x,P);                     % generate psychological states
y    = M.h(s,P);                     % and threshold symptom scores

return

% plotting sub function
%==========================================================================
function spm_DCM_plot(DCM)
 
% generate true and posterior expectations
%--------------------------------------------------------------------------
[Y,S,X]  = spm_psychosis_gen(DCM.P ,DCM.M,DCM.U);
[y,s,x]  = spm_psychosis_gen(DCM.Ep,DCM.M,DCM.U);

% plot symptoms, psychopathology and pathophysiology
%--------------------------------------------------------------------------
subplot(3,2,1), plot(Y,':'); hold on, set(gca,'ColorOrderIndex',1)
plot(y), hold off
title('Symptoms scores','fontsize',16), xlabel('time')
subplot(3,2,3), plot(S,':'); hold on, set(gca,'ColorOrderIndex',1)
plot(s), hold off
title('Latent psychopathology','fontsize',16), xlabel('time')
subplot(3,2,4), plot(X,':'); hold on, set(gca,'ColorOrderIndex',1), 
plot(x), hold off
title('Latent pathophysiology','fontsize',16), xlabel('time')

% supplement with trajectories in phase space
%--------------------------------------------------------------------------
subplot(3,2,5), cla
plot(S(:,1),S(:,2)), hold on;
plot(s(:,1),s(:,2))
title('Psychopathology','fontsize',16), box off
legend({'true','inferred'})

subplot(3,2,6), cla
plot3(X(:,1),X(:,2),X(:,3)), hold on;
plot3(x(:,1),x(:,2),x(:,3))
title('Pathophysiology','fontsize',16), box off

% supplement with image format of symptom scores
%--------------------------------------------------------------------------
subplot(6,2,2), imagesc(Y')
title('Observed and predicted scores','fontsize',16)
subplot(6,2,4), imagesc(y'), xlabel('time')



