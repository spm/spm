function DEMO_model_reduction_ERP
% Illustration of (post hoc)the neuronal mass model optimisation
%__________________________________________________________________________
% This demonstration routine illustrates the post-hoc optimisation of
% dynamic causal models for event related responses. To assess performance
% in relation to ground truth, it uses simulated data. We will simulate a
% simple two source model with exogenous input to the first source and
% reciprocal (extrinsic) connections between the two sources. the ERPs are
% simulated and two conditions, where the second condition induces a change
% in the intrinsic coupling of the first source and the forward extrinsic
% coupling. We then explore a simple model space; created by increasing the
% precision of shrinkage priors on the intrinsic condition specific effect.
% Because this effect was responsible for generating the data, we expect
% the free energy (log evidence) to fall as the shrinkage covariance falls
% to 0). Crucially, we compare and contrast the estimates of the free
% energy (and parameter estimates) using an explicit inversion of the
% reduced models (with tighter shrinkage priors) and a post-hoc model
% reduction procedure - that is computationally more efficient and
% robust to local minima.
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEMO_model_reduction_ERP.m 7679 2019-10-24 15:54:07Z spm $

% model specification - a simple two source model with two electrodes
% =========================================================================
rng('default')

Nc    = 2;                                        % number of channels
Ns    = 2;                                        % number of sources

options.spatial  = 'LFP';
options.model    = 'ERP';
options.analysis = 'ERP';
M.dipfit.model   = options.model;
M.dipfit.type    = options.spatial;
M.dipfit.Nc      = Nc;
M.dipfit.Ns      = Ns;

% sspecify connectivity - reciprocal connections with condition specific
% changes in intrinsic and extrinsic connectivity
%--------------------------------------------------------------------------
A{1}    = [0 0; 1 0];
A{2}    = [0 1; 0 0];
A{3}    = [0 0; 0 0];
B{1}    = [1 0; 1 0];
C       = [1; 0];

[pE,pC] = spm_dcm_neural_priors(A,B,C,options.model);
[gE,gC] = spm_L_priors(M.dipfit);
[x,f]   = spm_dcm_x_neural(pE,options.model);

% hyperpriors (assuming a high signal to noise)
%--------------------------------------------------------------------------
hE      = 6;
hC      = 1/128;

% create model
%--------------------------------------------------------------------------
M.IS   = 'spm_gen_erp';
M.G    = 'spm_lx_erp';
M.f    = f;
M.x    = x;
M.pE   = pE;
M.pC   = pC;
M.gE   = gE;
M.gC   = gC;
M.hE   = hE;
M.hC   = hC;
M.m    = length(B);
M.n    = length(spm_vec(M.x));
M.l    = Nc;
M.ns   = 64;

% create input structure
%--------------------------------------------------------------------------
dt     = 4/1000;
pst    = (1:M.ns)*dt;
M.ons  = 64;
M.dur  = 16;
U.dt   = dt;
U.X    = [0; 1];

% specified true connectivity (P) and spatial parameters (G) - with
% condition specific effects on the intrinsic connectivity of the first
% source and its forward extrinsic connection
%--------------------------------------------------------------------------
P      = pE;
G      = gE;
P.B{1} = [-1/4 0; 1/2 0];


% generate neuronal response and data
%--------------------------------------------------------------------------
x     = spm_gen_erp(P,M,U);                 % neuronal response
L     = spm_lx_erp(G,M.dipfit);             % lead field
V     = spm_sqrtm(spm_Q(1/2,M.ns));         % square root of noise covariance
for i = 1:length(x)
    n    = exp(-hE/2)*V*randn(M.ns,Nc);     % noise
    s{i} = x{i}*L';                         % signal
    y{i} = s{i} + n;                        % data (signal plus noise)
end

% data structure specification
%--------------------------------------------------------------------------
Y.y   = y;
Y.Q   = {spm_Q(1/2,M.ns,1)};
Y.dt  = dt;
Y.pst = pst;


% display
%--------------------------------------------------------------------------
spm_figure('Getwin','Figure 1');
subplot(2,1,1)
plot(pst,x{1},'r',pst,x{2},'b')
xlabel('time');ylabel('amplitude');
title('Hidden neuronal states','FontSize',16)

subplot(2,1,2)
plot(pst,s{1},':r',pst,s{2},':b',pst,y{1},'r',pst,y{2},'b')
xlabel('time');ylabel('amplitude');
title('Observed response','FontSize',16)



% Invert model under increasing shrinkage priors on the condition specific
% change in the intrinsic (B) parameter of the first source. This range
% specified by alpha. Because the true value is non-zero, we expect  free-
% energy to decrease when the prior covariance falls to 0 and this
% parameter is effectively eliminated. The ensuing estimates are obtained
% while optimising the parameters of the spatial model and the  (neuronal)
% extrinsic parameter..
% =========================================================================


% fix all (neuronal) parameters (except those of interest)
% -------------------------------------------------------------------------
M.pC           = spm_unvec(spm_vec(pC)*0,pC);
M.pC.B{1}(1,1) = 1/8;
M.pC.B{1}(2,1) = 1/8;

% vary the prior covariance
% -------------------------------------------------------------------------
alpha = exp(-8:2);
for i = 1:length(alpha)   

    % reset shrinkage prior
    % ---------------------------------------------------------------------
    M.pC.B{1}(1,1) = alpha(i);

    % full inversion
    % ---------------------------------------------------------------------
    [Ep,Eg,Cp,Cg,S,F] = spm_nlsi_N(M,U,Y);
    Ep_all{i} = spm_vec(Ep);
    Cp_all{i} = diag(Cp);
    F_all(i)  = F;
    
end

% Post-hoc reduction for the full (complex) model: this evaluates the free
% energy and posterior distribution over the parameters, given just the
% full prior and posterior - for any required reduced prior.
% -------------------------------------------------------------------------
for i = 1:length(alpha)
    
    rC           = M.pC;
    rC.B{1}(1,1) = alpha(i);

    [F,rEp,rCp]  = spm_log_evidence_reduce(Ep,Cp,M.pE,M.pC,M.pE,rC);
    rEp_all{i}   = spm_vec(rEp);
    rCp_all{i}   = diag(rCp);
    rF_all(i)    = F;
    
end


% compare full inversion and model reduction
% -------------------------------------------------------------------------
spm_figure('Getwin','Figure 2');

subplot(2,1,1)
semilogx(alpha,F_all - F_all(end),alpha,rF_all - rF_all(end))
xlabel('prior covariance');ylabel('relative log evidence');
legend('full','posthoc');
title('Full and reduced log evidence','FontSize',16)


% compare parameter estimates
% -------------------------------------------------------------------------
j     = find(spm_vec(M.pC));
pP    = spm_vec(P);
for i = 1:length(alpha)   

    fQp(i,:) = Ep_all{i}(j);
    fCp(i,:) = Cp_all{i}(j);
    rQp(i,:) = rEp_all{i}(j);
    vCp(i,:) = rCp_all{i}(j);
    
    Pp(i,:)  = pP(j);
end

subplot(2,1,2)
spm_plot_ci(rQp',vCp'),  hold on
plot(fQp,'--'),          hold on
plot(Pp,'-.'),           hold off

xlabel('prior log-covariance');ylabel('difference in log evidence');
title('MAP estimates (full -- reduced  - true -.)','FontSize',16)
set(gca,'XTickLabel',log(alpha))


return


% free energy landscape:
% Here, we evaluate the free energy, which is a functional of the data and
% conditional or posterior expectations (noting that posterior precisions
% are functions of the expectations). The free energy can be computed in a
% simple way by inverting the model using a single iteration. In what
% follows, we evaluate the free energy over a range of the two (intrinsic
% and extrinsic) coupling parameters (at the true values of the remaining
% parameters)
% =========================================================================
beta      = linspace(-1,1,32);                       % range of parameters
M.nograph = 1;
M.Nmax    = 1;
M.Gmax    = 1;
M.Hmax    = 1;

M.P   = P;
M.Q   = G;
for i = 1:length(beta)
    for j = 1:length(beta)
        
        % specify posterior expectations (and implicitly precisions)
        % -----------------------------------------------------------------
        M.P           = P;
        M.P.B{1}(1,1) = beta(i);
        M.P.B{1}(2,1) = beta(j);
        
        % evaluate free energy (without updating expectations)
        % -----------------------------------------------------------------
        [Ep,Eg,Cp,Cg,S,F] = spm_nlsi_N(M,U,Y);
        FF(i,j)           = F;
        
        % record conditional covariance at maximum
        % -----------------------------------------------------------------
        if F >= max(FF(:)), CP  = Cp; end

    end
end

% apply Occam's window
% -------------------------------------------------------------------------
FF    = FF - max(FF(:));
F     = max(FF,-64)

% Free energy landscape and associated conditional covariance: the
% conditional covariance  is the inverse precision - which is proportional
% to the curvature of the variational energy (equivalent to the curvature
% of the free energy functional of posterior or conditional expectations)
% -------------------------------------------------------------------------
spm_figure('Getwin','Figure 3');

j     = find(spm_vec(M.pC));
subplot(2,1,1)
imagesc(beta,beta,F),       hold on
contour(beta,beta,F,8,'m'), hold on
plot(pP(j(2)),pP(j(1)),'.r','MarkerSize',32), hold off
xlabel('extrinsic');ylabel('intrinsic');
title('Free-energy landscape','FontSize',16)
axis square

subplot(2,1,2)
imagesc(CP(j,j))
xlabel('parameters');ylabel('parameters');
title('Posterior covariance','FontSize',16)
axis square

