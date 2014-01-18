function DEM_demo_modes_fMRI
% Demonstration of spectral DCM for fMRI with eigenvector constraints
%__________________________________________________________________________
% This demonstration routine illustrates the inversion of resting state
% fMRI timeseries using a generative model of the adjacency matrix. This
% model is based upon the eigenmodes of the functional connectivity matrix,
% which are the eigenvectors of the effective connectivity matrix or
% Jacobian - assuming the effective connectivity is symmetrical. This means
% it is only necessary to estimate the eigenvalues; in other words, one
% unknown parameter per node.
%
% Simulated timeseries are generated and inverted under typical priors.
% This routine then performs a model space search over the decay rates of
% stable (dissipative) modes and the number of unstable modes.
% This illustrates: (i) the increase in model evidence afforded by
% hierarchical constraints (when they are true) and (ii) the identification
% of the principle modes underlying connectivity.
%
% NB: To see the model optimisation delete the 'return' at about line 200
% 
% see also: DEM_demo_connectivity_fMRI
%           spm_dcm_fmri_mode_gen
%           spm_dcm_fmri_mode
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_demo_modes_fMRI.m 5838 2014-01-18 18:40:37Z karl $


% Simulate timeseries
%==========================================================================
rng('default')

% DEM Structure: create random inputs
% -------------------------------------------------------------------------
ED    = 3;                               % embedding dimensional
DE    = 1;                               % log decay of stable modes (secs)
T     = 1024;                            % number of observations (scans)
TR    = 2;                               % repetition time or timing
n     = 6;                               % number of regions or nodes
t     = (1:T)*TR;                        % observation times


% empirical covariance matrix
% -------------------------------------------------------------------------
Cy    = [
    1.0073   -0.4669   -0.3647    0.5984    0.4378    0.6177
   -0.4669    0.8979    0.2469   -0.4208   -0.2764   -0.4382
   -0.3647    0.2469    0.4027   -0.2821   -0.1795   -0.2702
    0.5984   -0.4208   -0.2821    0.7553    0.4419    0.6125
    0.4378   -0.2764   -0.1795    0.4419    0.4634    0.3947
    0.6177   -0.4382   -0.2702    0.6125    0.3947    0.9421];


% modes and exponents
% -------------------------------------------------------------------------
u     = spm_svd(Cy);


% priors
% -------------------------------------------------------------------------
options.nmax       = 8;                % effective number of nodes
options.Nmax       = 64;               % maximum number of iterationsno
options.two_state  = 0;
options.induced    = 1;
options.stochastic = 0;
options.nonlinear  = 0;
options.decay      = DE;
options.precision  = 0;


A   = full(sparse(1:ED,1,1,n,1));
B   = zeros(n,n,0);
C   = zeros(n,n);
D   = zeros(n,n,0);


% true parameters (reciprocal connectivity)
% -------------------------------------------------------------------------
[pP,pC,x]  = spm_dcm_fmri_priors(A,B,C,D,options);
pP.A(1:ED) = linspace(2,0,ED); 
pA         = spm_dcm_fmri_mode_gen(pP.A,u); 
pP.C       = eye(n,n);
pP.transit = randn(n,1)*exp(-4);
disp(pA)

% simulate response to endogenous fluctuations
%==========================================================================

% integrate states
% -------------------------------------------------------------------------
U.u     = spm_rand_mar(T,n,1/2)/4;      % endogenous fluctuations
U.dt    = TR;
M.f     = 'spm_fx_fmri';
M.x     = x;
M.modes = u;
x       = spm_int_J(pP,M,U);

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
i = 1:256;
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

DCM.a    = sparse(1:ED,1,1,n,1);
DCM.b    = zeros(n,n,0);
DCM.c    = zeros(n,1);
DCM.d    = zeros(n,n,0);

DCM.Y.y  = y + e;
DCM.Y.dt = TR;
DCM.U.u  = zeros(T,1);
DCM.U.dt = TR;

% nonlinear system identification (Variational Laplace)
% =========================================================================

% classical (with symmetry constraints on effective connectivity)
% -------------------------------------------------------------------------
DCM.M.symmetry = 1;

DCM.a = ones(n,n);
DCM   = spm_dcm_fmri_csd(DCM);
FC    = DCM.F;

% with explicit (modal) stability constraints
% -------------------------------------------------------------------------
DCM.a = sparse(1:ED,1,1,n,1);
DEM   = spm_dcm_fmri_csd(DCM);


% summary
% -------------------------------------------------------------------------
spm_figure('Getwin','Figure 2'); clf

[Ep Cp] = spm_dcm_fmri_mode_gen(DEM.Ep.A,u,DEM.Cp(1:n,1:n)); 
Pp      = spm_dcm_fmri_mode_gen(pP.A,u);
Qp      = DCM.Ep.A; Qp = Qp - diag(diag(Qp) + exp(diag(Qp))/2);

subplot(2,1,1); hold off
spm_plot_ci(Ep(:),Cp), hold on
bar(Pp(:),1/4,'k'), hold off
title('True and MAP connections','FontSize',16)
axis square

j   = find(eye(n,n));

subplot(2,1,2); cla
plot(Pp(:),Ep(:),'b.','MarkerSize',32),    hold on
plot(Pp(:),Qp(:),'c.','MarkerSize',16),    hold on
plot(Pp(j),Ep(j),'r.','MarkerSize',32), hold on
plot(Pp(j),Qp(j),'m.','MarkerSize',16), hold on
plot([-2 1],[-2 1],'k--'), hold on
title('MAP vs. true','FontSize',16)
xlabel('True')
ylabel('Estimate')
axis square
legend ('hierarchical','conventional')


% proximity space
%==========================================================================
spm_figure('Getwin','Figure (MAP)'); clf

Ep   = DEM.Ep;
Cp   = spm_unvec(spm_vec(diag(DEM.Cp)),Ep);

subplot(2,2,1)
spm_dcm_graph_functional(u')
title('True','FontSize',16)

subplot(2,2,2)
bar(exp(pP.A)),axis square
xlabel('Mode')
ylabel('Time Constants (secs)')
title('True','FontSize',16)

qu = DEM.M.modes;
qu = qu*diag(diag(sign(qu'*u)));

subplot(2,2,3)
spm_dcm_graph_functional(qu')
title('Estimated','FontSize',16)

subplot(2,2,4)
spm_plot_ci(Ep.A,Cp.A,[],[],'exp'),axis square
xlabel('Mode')
ylabel('Time Constants (secs)')
title('Estimated','FontSize',16)



% return if called as a demo - otherwise perform model search
%--------------------------------------------------------------------------
if ~exist('DCM_stochastic.mat','file')
    return
end


% search over decay of stable causes
%==========================================================================
DCM.a = sparse(1:ED,1,1,n,1);
V     = 0:1/4:2;
F     = [];
R     = [];
for i = 1:length(V)
    
    % invert
    %======================================================================
    DCM.options.decay    = V(i);
    DEM = spm_dcm_fmri_csd(DCM);
    
    % RMS
    % ---------------------------------------------------------------------
    R(i) = paper_rms(DEM.Ep.A,pP.A,u);
    
    % free energy
    % ---------------------------------------------------------------------
    F(i) = DEM.F;
    
end


% summary
% -----------------------------------------------------------------
spm_figure('Getwin','Figure 3'); clf

subplot(2,2,1);
bar(V,F - min(F(:)) + 16)
title('Log-evidence and embedding','FontSize',16)
xlabel('Log-decay')
ylabel('Free energy')
axis square

subplot(2,2,3);
bar(V,R), hold on
plot(V,V*0 + 0.05,'r:'), hold off
title('RMS error','FontSize',16)
xlabel('Log-decay')
ylabel('Root mean square error')
axis square



% search over embedding dimension
%==========================================================================
DCM.options.decay = DE;
D     = 1:5;
DF    = [];
DR    = [];
for i = 1:length(D)

    % invert
    %======================================================================
    DCM.a = sparse(1:D(i),1,1,n,1);
    DEM   = spm_dcm_fmri_csd(DCM);
        
    % RMS
    % ---------------------------------------------------------------------
    DR(i)   = paper_rms(DEM.Ep.A,pP.A,u);
    
    % free energy
    % ---------------------------------------------------------------------
    DF(i)   = DEM.F;
    
end
    

% summary
% -----------------------------------------------------------------
spm_figure('Getwin','Figure 3');

RF  = DF - min(DF(:)) + 16;

subplot(2,2,2); cla
bar(D,RF),                                              hold on
plot(D,D - D + max(RF),'r',D,D - D + max(RF) - 3,'r:'), hold off
title('Log-evidence and embedding','FontSize',16)
xlabel('Embedding dimension')
ylabel('Free energy')
axis square

subplot(2,2,4);
bar(D,DR),               hold on
plot(D,D*0 + 0.05,'r:'), hold off
title('RMS error','FontSize',16)
xlabel('Embedding dimension')
ylabel('Root mean square error')
axis square


% load empirical DCM for search over precision and embedding dimension
%==========================================================================
try
    load DCM_stochastic
catch
    return
end

n     = DCM.n;

DCM.b = zeros(n,n,0);
DCM.d = zeros(n,n,0);
DCM.options = options;

% search over embedding dimension
%==========================================================================
DCM.options.decay = DE;
eDF   = [];
for i = 1:length(D)

    % invert
    %======================================================================
    DCM.a = sparse(1:D(i),1,1,n,1);
    DEM   = spm_dcm_fmri_csd(DCM);
            
    % free energy
    % ---------------------------------------------------------------------
    eDF(i)   = DEM.F;
    
end

% save parameter estimates and get dimension that maximises free energy
% -------------------------------------------------------------------------
Ep    = DEM.Ep;
Cp    = spm_unvec(spm_vec(diag(DEM.Cp)),Ep);
[f j] = max(eDF);

% search over precision of hidden causes
%==========================================================================
DCM.a = sparse(1:D(j),1,1,n,1);
eF    = [];
for i = 1:length(V)
    
    % invert
    %======================================================================
    DCM.options.decay    = V(i);
    DEM = spm_dcm_fmri_csd(DCM);
    
    % free energy
    % ---------------------------------------------------------------------
    eF(i) = DEM.F;
    
end


% summary
% -------------------------------------------------------------------------
spm_figure('Getwin','Figure 4'); clf

subplot(2,2,1);
bar(V,eF - min(eF(:)) + 16)
title('Decay (empirical)','FontSize',16)
xlabel('Log-decay ')
ylabel('Free energy')
axis square


% summary
% -------------------------------------------------------------------------
spm_figure('Getwin','Figure 4');

eRF  = eDF - min(eDF(:)) + 16;

subplot(2,2,2); cla
bar(D,eRF),                                               hold on
plot(D,D - D + max(eRF),'r',D,D - D + max(eRF) - 3,'r:'), hold off
title('Embedding (empirical)','FontSize',16)
xlabel('Embedding dimension')
ylabel('Free energy')
axis square


% and save matlab file
% -------------------------------------------------------------------------
save paper

% get functional space
%==========================================================================
spm_figure('Getwin','Figure 5'); clf

subplot(2,2,1)
spm_dcm_graph_functional(DEM.M.modes')
title('Scaling space','FontSize',16)

subplot(2,2,2)
spm_plot_ci(Ep.A,Cp.A,[],[],'exp'),axis square
xlabel('Mode')
ylabel('Time Constants (secs)')
title('Time Constants','FontSize',16)


% and compare functional and effective connectivity matrices
%--------------------------------------------------------------------------
spm_figure('Getwin','Figure 6'); clf

[qu,E,F] = spm_dcm_fmri_mode(Ep.A,u);
F        = spm_cov2corr(F);
EC       = spm_vec(E - diag(diag(E)));
FC       = spm_vec(F - diag(diag(F)));
EC       = EC(~~EC);
FC       = FC(~~FC);

subplot(3,2,1), imagesc(E)
xlabel('Node'), ylabel('Node'), axis square
title('Effective connectivity','FontSize',16)

subplot(3,2,2), imagesc(F)
xlabel('Node'), ylabel('Node'), axis square
title('Functional connectivity','FontSize',16)

subplot(3,2,3), hist(abs(EC),8)
xlabel('Connection strength (Hz)'), ylabel('Frequency')
title('Distribution','FontSize',16), axis square

subplot(3,2,4), hist(abs(FC),8)
xlabel('Correlation strength (Hz)'), ylabel('Frequency')
title('Distribution','FontSize',16), axis square

subplot(3,2,5), imagesc(abs(E) < 0.3)
xlabel('Node'), ylabel('Node'), axis square
title('Strong connections','FontSize',16)

subplot(3,2,6), imagesc(abs(F) < 0.3)
xlabel('Node'), ylabel('Node'), axis square
title('Strong correlations','FontSize',16)


return



function rms = paper_rms(A,B,u)
% Root mean square difference in (% extrinsic connectivity)
% -------------------------------------------------------------------------
A   = spm_dcm_fmri_mode_gen(A,u);
B   = spm_dcm_fmri_mode_gen(B,u);
D   = A - B;
rms = sqrt(mean(D(:).^2));

return


% notes for power law scaling
%==========================================================================

spm_figure('Getwin','Figure 7'); clf

% frequencies
%--------------------------------------------------------------------------
dw   = 1/32;
w    = (-32:32)*dw;

% time constants
%--------------------------------------------------------------------------
t    = [128 64 32 16]/64;

% plot Lorentzian functions
%--------------------------------------------------------------------------
G    = 0;
for i = 1:length(t)
    g(:,i) = 1./((2*pi*w).^2 + t(i)^(-2));
end

subplot(2,1,1)
plot(w,g);
xlabel('Frequencies')
ylabel('Spectral density')
title('Lorentzian functions','FontSize',16);

% time constants
%--------------------------------------------------------------------------
a    = 2;
t    = (1:1024)/256;
p    = t.^-a;
p    = p/sum(p);

% plot Lorentzian functions
%--------------------------------------------------------------------------
w    = (1:256)*dw;
G    = 0;
for i = 1:length(t)
    G      = G + p(i)*1./((2*pi*w).^2 + t(i)^(-2));
end

subplot(2,2,3)
plot(w,G);
xlabel('Frequencies')
ylabel('Spectral density')
title('Mixture of Lorentzians','FontSize',16);

LG   =  log(G);
PL   = -log(w);
PL   = PL - PL(end) + LG(end);

subplot(2,2,4)
plot(log(w),LG,log(w),PL);
xlabel('Frequencies')
ylabel('Spectral density')
title('Power law scaling','FontSize',16);



