function DEM_demo_connectivity_fMRI
% Demonstration of DCM for fMRI–CSD with hierarchical constraints
%__________________________________________________________________________
% This demonstration routine illustrates the inversion of resting state
% fMRI timeseries using a generative model of the adjacency matrix. This
% model is based upon an embedding space of dimensions the in which the
% (log) connectivity among nodes is a (radial basis) function of their
% metric separation. This generative model of connectivity requires a
% hierarchical constraints on the edges and therefore uses the expectation
% aand maximisation step of dynamic expectationmaximisation. Here, the
% hidden causes at the first level are the effective connectivity and the
% hidden causes at the second level are the locations in embedding states.
%
% simulated timeseries are generated and inverted under typical priors.
% This routine that performs a model space search over precisions on the
% hierarchical constraints and the dimensionality of the embedding space.
% This illustrates: (i) the increase in model evidence afforded by
% hierarchical constraints (when they are true) and (ii) the optimal
% prior precision that reflects the amplitude of random variations in
% connectivity about the constraints. (iii) Finally,the search over moral
% dimension illustrates how Bayesian model comparison can identify the
% dimensionality of the metric space generating hierarchical connectivity.
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_demo_connectivity_fMRI.m 5708 2013-10-22 09:20:59Z karl $

% Simulate timeseries
%==========================================================================
rng('default')

% DEM Structure: create random inputs
% -------------------------------------------------------------------------
D   = 3;                               % embedding dimensional
T   = 512;                             % number of observations (scans)
TR  = 2;                               % repetition time or timing
n   = 5;                               % number of regions or nodes
t   = (1:T)*TR;                        % observation times
v.a = 3;                               % log amplitude
v.x = randn(D,n);                      % location

% priors
% -------------------------------------------------------------------------
options.nmax       = 8;                % effective number of notes

options.two_state  = 1;
options.induced    = 1;
options.stochastic = 0;
options.nonlinear  = 0;
options.embedding  = 3;

A   = ones(n,n);
B   = zeros(n,n,0);
C   = zeros(n,n);
D   = zeros(n,n,0);


% true parameters (reciprocal connectivity)
% -------------------------------------------------------------------------
[pE,pC,x]  = spm_dcm_fmri_priors(A,B,C,D,options);
pP         = spm_dcm_fmri_graph_gen([],v,pE);
pP.A       = pP.A + randn(n,n)*exp(-2);
pP.C       = eye(n,n);
pP.transit = randn(n,1)/16;

% simulate response to endogenous fluctuations
%==========================================================================

% integrate states
% -------------------------------------------------------------------------
U.u  = spm_rand_mar(T,n,1/2)/2;      % endogenous fluctuations
U.dt = TR;
M.f  = 'spm_fx_fmri';
M.x  = x;
x    = spm_int_J(pP,M,U);

% haemodynamic observer
% -------------------------------------------------------------------------
for i = 1:T
    y(i,:) = spm_gx_fmri(spm_unvec(x(i,:),M.x),[],pP)';
end

% observation noise process
% -------------------------------------------------------------------------
e    = spm_rand_mar(T,n,1/2)/8;

% show simulated response
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


% nonlinear system identification (DCM for CSD) over subjects
%==========================================================================
DCM.options = options;

DCM.a    = ones(n,n);
DCM.b    = zeros(n,n,0);
DCM.c    = zeros(n,1);
DCM.d    = zeros(n,n,0);

DCM.Y.y  = y + e;
DCM.Y.dt = TR;
DCM.U.u  = zeros(T,1);
DCM.U.dt = TR;


% nonlinear system identification (Variational Laplace)
% =========================================================================

% classical
% -------------------------------------------------------------------------
DCM  = spm_dcm_fmri_csd(DCM);

% hierarchical
% -------------------------------------------------------------------------
DEM  = spm_dcm_fmri_csd_DEM(DCM);


% summary
% -------------------------------------------------------------------------
spm_figure('Getwin','Figure 2'); clf

subplot(2,1,1); hold off
j   = 1:numel(DCM.Ep.A);
Ep  = exp(DEM.Ep.A(:))/8;
Qp  = exp(DCM.Ep.A(:))/8;
Cp  = exp(DCM.Cp(j,j)/64)/8;
Pp  = exp(pP.A(:))/8;
spm_plot_ci(Ep,Cp), hold on
bar(Pp,1/2), hold off
title('True and MAP connections','FontSize',16)
axis square

subplot(2,2,3); cla
plot(Pp,Ep,'.','MarkerSize',16), hold off
title('MAP vs. true','FontSize',16)
xlabel('true')
ylabel('estimate')
axis square

subplot(2,2,4); cla
plot(Pp,Qp,'.','MarkerSize',16), hold off
title('classical vs. true ','FontSize',16)
xlabel('true')
ylabel('estimate')
axis square

% load empirical DCM for search over precision and embedding dimension
%==========================================================================
empirical = 1;
if empirical
    load DCM_ATT
    
    n     = DCM.n;
    DCM.a = ones(n,n);
    DCM.b = zeros(n,n,0);
    DCM.d = zeros(n,n,0);
    
    pP.A  = zeros(n,n);
end


% search over precision of hidden causes
%==========================================================================
DCM.options = options;
DCM.options.v = [4, -4];

V     = linspace(1,8,8);
F     = [];
RMS   = [];
for i = 1:length(V)
    
    DCM.options.v(1) = V(i);
    
    % invert
    %======================================================================
    DCM.options.v(2) = -4;
    DEM = spm_dcm_fmri_csd_DEM(DCM);
    
    % root mean square error
    % ---------------------------------------------------------------------
    dp             = (exp(DEM.Ep.A) - exp(pP.A))/8;
    RMS(end + 1,1) = sqrt(mean(dp(~~dp).^2));
    
    % free energy
    % ---------------------------------------------------------------------
    F(end + 1,1)   = DEM.F;
    
    
    % repeat with precise full priors
    %======================================================================
    DCM.options.v(2) = 16;
    DEM = spm_dcm_fmri_csd_DEM(DCM);
    
    % root mean square error
    % ---------------------------------------------------------------------
    dp             = (exp(DEM.Ep.A) - exp(pP.A))/8;
    RMS(end,2)     = sqrt(mean(dp(~~dp).^2));
    
    % free energy
    % ---------------------------------------------------------------------
    F(end,2)       = DEM.F;
    
end


% summary
% -----------------------------------------------------------------
spm_figure('Getwin','Figure 3'); clf

subplot(2,2,1); cla
bar(F), hold on
plot([4 4],[min(F(:)) 1.2*max(F(:))],'r:','LineWidth',4), hold off
title('log-evidence and precision','FontSize',16)
xlabel('prior precision')
ylabel('free energy')
axis square

subplot(2,2,2);
bar(RMS), hold on
plot([1 length(V)],[0.05 0.05],'r:','LineWidth',4), hold off
title('root mean square error','FontSize',16)
xlabel('prior precision')
ylabel('RMS')
axis square


% search over embedding dimension
%==========================================================================
DCM.options.v = [4, -4];

D     = 1:6;
DF    = [];
DRMS  = [];
for i = 1:length(D)

    % invert
    %======================================================================
    DCM.options.embedding = D(i);
    DEM = spm_dcm_fmri_csd_DEM(DCM);
        
    % root mean square error
    % ---------------------------------------------------------------------
    dp              = (exp(DEM.Ep.A) - exp(pP.A))/8;
    DRMS(end + 1,1) = sqrt(mean(dp(~~dp).^2));
    
    % free energy
    % ---------------------------------------------------------------------
    DF(end + 1,1)   = DEM.F;
    
end
    

% summary
% -----------------------------------------------------------------
spm_figure('Getwin','Figure 3');

subplot(2,2,3); cla
bar(DF), hold on
title('log-evidence and embedding','FontSize',16)
xlabel('embedding dimension')
ylabel('free energy')
set(gca,'YLim',[min(DF - 16) max(DF) + 16])
axis square

subplot(2,2,4);
bar(DRMS), hold on
title('root mean square error','FontSize',16)
xlabel('embedding dimension')
ylabel('RMS')
set(gca,'YLim',[min(DRMS - std(DRMS)) max(DRMS) + std(DRMS)])
axis square

% and save matlab file
% -----------------------------------------------------------------
if empirical
    save empirical
else
    save simulations
end


return

% NOTES: illustrate the ill-posed nature of the problem
%==========================================================================
M     = DCM.M;
U     = DCM.U;
M.x   = zeros(n,5); 

nA    = 32;
pA    = linspace(-.4,.4,nA);
Y     = [];
P     = [];
for i = 1:nA
    for j = 1:nA
        
        % map from parameter space to data space
        %------------------------------------------------------------------
        pp           = pP;
        pp.A(1,2)    = pA(i);
        pp.A(2,1)    = pA(j);
        Y(:,end + 1) = spm_vec(spm_csd_fmri_mtf(pp,M,U));
        P(:,end + 1) = spm_vec(pp.A);
        
    end
end

% distance measures
%--------------------------------------------------------------------------
Up      = P([2 (n + 1)],:)';
[Uy Sy] = spm_svd(spm_detrend(Y'));
Uy      = real(Uy);

Cp    = Up;
for i = 1:2
    Cp(:,i) = Up(:,i) - min(Up(:,i));
    Cp(:,i) = 0.001 + Cp(:,i)./(max(Cp(:,i))*1.1);
end

spm_figure('Getwin','Figure 4'); clf
subplot(2,1,1), cla
for  i = 1:nA*nA
    plot(Up(i,1),Up(i,2),'.','Markersize',32,'Color',[1/2 Cp(i,1) Cp(i,2)]), hold on
end
axis square
title('Parameter space','FontSize',16)
xlabel('Forward connection')
ylabel('Backward connection')
axis square

subplot(2,1,2), cla
for  i = 1:nA*nA
    plot3(Uy(i,1),Uy(i,2),Uy(i,3),'.','Markersize',32,'Color',[1/2 Cp(i,1) Cp(i,2)]), hold on
end
axis square
title('Data space','FontSize',16)
xlabel('1st PC')
ylabel('2nd PC')
ylabel('3rd PC')
axis square

return


