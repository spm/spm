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
% $Id: DEM_demo_connectivity_fMRI.m 5736 2013-11-10 13:17:10Z karl $

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
v.a = [2; -1];                         % log amplitude
v.x = randn(D,n)/2;                    % location

% priors
% -------------------------------------------------------------------------
options.nmax       = 8;                % effective number of notes

options.two_state  = 0;
options.induced    = 1;
options.stochastic = 0;
options.nonlinear  = 0;
options.embedding  = 3;
options.backwards  = 0;
options.v          = 4;

A   = ones(n,n);
B   = zeros(n,n,0);
C   = zeros(n,n);
D   = zeros(n,n,0);


% true parameters (reciprocal connectivity)
% -------------------------------------------------------------------------
[pE,pC,x]  = spm_dcm_fmri_priors(A,B,C,D,options);
pP         = spm_dcm_fmri_graph_gen([],v,pE);
pP.A       = pP.A + randn(size(pP.A))*exp(-4); disp(pP.A)
pP.C       = eye(n,n);
pP.transit = randn(n,1)*exp(-4);

% simulate response to endogenous fluctuations
%==========================================================================

% integrate states
% -------------------------------------------------------------------------
U.u  = spm_rand_mar(T,n,1/2)/8;      % endogenous fluctuations
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
e    = spm_rand_mar(T,n,1/2)/4;

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
j   = 1:numel(pP.A);
Ep  = DEM.Ep.A(j);
Qp  = DCM.Ep.A(j);
Cp  = DCM.Cp(j,j);
Pp  = pP.A(:);
spm_plot_ci(Ep(:),Cp), hold on
bar(Pp,1/2,'b'), hold off
title('True and MAP connections','FontSize',16)
axis square


subplot(2,2,3); cla
plot(Pp,Ep,'r.','MarkerSize',16), hold off
title('MAP vs. true','FontSize',16)
xlabel('true')
ylabel('estimate')
axis square

subplot(2,2,4); cla
plot(Pp,Qp,'r.','MarkerSize',16), hold off
title('classical vs. true ','FontSize',16)
xlabel('true')
ylabel('estimate')
axis square


% search over precision of hidden causes
%==========================================================================
V     = 4:8;
F     = [];
R     = [];
for i = 1:length(V)
    
    % invert
    %======================================================================
    DCM.options.v         = V(i);
    DCM.options.embedding = 3;
    DEM = spm_dcm_fmri_csd_DEM(DCM);
    
    % correlation
    % ---------------------------------------------------------------------
    R(end + 1,1) = corr(DEM.Ep.A(:),pP.A(:));
    
    % free energy
    % ---------------------------------------------------------------------
    F(end + 1,1) = DEM.F;
    
    
    % repeat with precise full priors
    %======================================================================
    DCM.options.embedding = 0;
    DEM = spm_dcm_fmri_csd_DEM(DCM);
    
    % correlation
    % ---------------------------------------------------------------------
    R(end,2)     = corr(DEM.Ep.A(:),pP.A(:));
    
    % free energy
    % ---------------------------------------------------------------------
    F(end,2)     = DEM.F;
    
end


% summary
% -----------------------------------------------------------------
spm_figure('Getwin','Figure 3'); clf

subplot(2,2,1); cla
bar(V,F), hold on
plot([4 4],[min(F(:)) 1.2*max(F(:))],'r:','LineWidth',4), hold off
title('log-evidence and precision','FontSize',16)
xlabel('prior precision')
ylabel('free energy')
axis square

subplot(2,2,3);
bar(V,R), hold on
plot([1 length(V)],[0.05 0.05],'r:','LineWidth',4), hold off
title('correlation','FontSize',16)
xlabel('prior precision')
ylabel('rho')
axis square


% search over embedding dimension
%==========================================================================
D     = 0:4;
DF    = [];
DR    = [];
for i = 1:length(D)

    % invert
    %======================================================================
    DCM.options.v         = 6;
    DCM.options.embedding = D(i);
    DEM = spm_dcm_fmri_csd_DEM(DCM);
        
    % correlation
    % ---------------------------------------------------------------------
    DR(end + 1,1)   = corr(DEM.Ep.A(:),pP.A(:));
    
    % free energy
    % ---------------------------------------------------------------------
    DF(end + 1,1)   = DEM.F;
    
end
    

% summary
% -----------------------------------------------------------------
spm_figure('Getwin','Figure 3');

subplot(2,2,2); cla
bar(D,DF), hold on
title('log-evidence and embedding','FontSize',16)
xlabel('embedding dimension')
ylabel('free energy')
set(gca,'YLim',[min(DF - 16) max(DF) + 16])
axis square

subplot(2,2,4);
bar(D,DR), hold on
title('correlation','FontSize',16)
xlabel('embedding dimension')
ylabel('rho')
set(gca,'YLim',[min(DR - std(DR)) max(DR) + std(DR)])
axis square


% load empirical DCM for search over precision and embedding dimension
%==========================================================================
load DCM_stochastic

n     = DCM.n;
DCM.a = ones(n,n);
DCM.b = zeros(n,n,0);
DCM.d = zeros(n,n,0);

% search over precision of hidden causes
%==========================================================================
DCM.options   = options;
eF    = [];
for i = 1:length(V)
    
    % invert
    %======================================================================
    DCM.options.v         = V(i);
    DCM.options.embedding = 1;
    DEM = spm_dcm_fmri_csd_DEM(DCM);
    
    % free energy
    % ---------------------------------------------------------------------
    eF(end + 1,1) = DEM.F;
    
    
    % repeat with precise full priors
    %======================================================================
    DCM.options.embedding = 0;
    DEM = spm_dcm_fmri_csd_DEM(DCM);
    
    % free energy
    % ---------------------------------------------------------------------
    eF(end,2) = DEM.F;
    
end


% summary
% -----------------------------------------------------------------
spm_figure('Getwin','Figure 4'); clf

subplot(2,2,1); cla
bar(V,eF), hold on
title('precision (empirical)','FontSize',16)
xlabel('prior precision ')
ylabel('free energy')
axis square


% search over embedding dimension
%==========================================================================
eDF   = [];
for i = 1:length(D)

    % invert
    %======================================================================
    DCM.options.v         = 6;
    DCM.options.embedding = D(i);
    DEM = spm_dcm_fmri_csd_DEM(DCM);
            
    % free energy
    % ---------------------------------------------------------------------
    eDF(end + 1,1)   = DEM.F;
    
end
    

% summary
% -------------------------------------------------------------------------
spm_figure('Getwin','Figure 4');

subplot(2,2,2); cla
bar(D,eDF), hold on
title('embedding (empirical)','FontSize',16)
xlabel('embedding dimension')
ylabel('free energy')
set(gca,'YLim',[min(eDF - 16) max(eDF) + 16])
axis square

    
% and save matlab file
% -------------------------------------------------------------------------
save paper



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


% graphics
%--------------------------------------------------------------------------
spm_figure('Getwin','Figure 5'); clf

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


