function DEM_demo_connectivity_fMRI
% Demonstration of DCM for CSD (fMRI) with simulated responses
%__________________________________________________________________________
% This demonstration compares generalised filtering and deterministic DCM
% (generating complex cross spectra) in the context of a nonlinear
% convolution (fMRI) model using simulated data. Here, the dynamic
% convolution model for fMRI responses is converted into a static
% non-linear model by generating not the timeseries per se but their
% second-order statistics – in the form of cross spectra and covariance
% functions. This enables model parameters to the estimated using the
% second order data features through minimisation of variational free
% energy. For comparison, the same data are inverted (in timeseries form)
% using generalised filtering. This example uses a particularly difficult
% problem – with limited data - to emphasise the differences. The
% generalised filtering uses the posteriors from the deterministic scheme
% as priors. This is equivalent to Bayesian parameter averaging using
% (orthogonal) second and first order data features.
%
% NB - the generalised filtering trakes much longer than the deterministic
% scheme
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_demo_connectivity_fMRI.m 5696 2013-10-15 19:10:26Z karl $

% Simulate timeseries
%==========================================================================
rng('default')
clear DCM
global DCM

% DEM Structure: create random inputs
% -------------------------------------------------------------------------
T   = 512;                             % number of observations (scans)
TR  = 2;                               % repetition time or timing
n   = 4;                               % number of regions or nodes
t   = (1:T)*TR;                        % observation times
v   = randn(3,n);                      % location
xyz = v + randn(3,n)/8;                % location

% priors
% -------------------------------------------------------------------------
options.nmax       = 8;               % effective number of notes

options.nonlinear  = 0;
options.two_state  = 0;
options.stochastic = 0;
options.centre     = 1;
options.induced    = 1;
options.backwards  = 1;

A   = ones(n,n);
B   = zeros(n,n,0);
C   = zeros(n,n);
D   = zeros(n,n,0);
pP  = spm_dcm_fmri_priors(A,B,C,D,options);
DCM.M.pE  = pP;


% true parameters (reciprocal connectivity)
% -------------------------------------------------------------------------
pP         = spm_dcm_fmri_graph_gen([],v,[]);
pP.C       = eye(n,n);
pP.transit = randn(n,1)/16;

% simulate response to endogenous fluctuations
%==========================================================================

% integrate states
% -------------------------------------------------------------------------
U.u  = spm_rand_mar(T,n,1/2)/2;      % endogenous fluctuations
U.dt = TR;
M.f  = 'spm_fx_fmri';
M.x  = sparse(n,5);
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
DCM.xY.xyz  = xyz;

DCM.a    = ones(n,n);
DCM.b    = zeros(n,n,0);
DCM.c    = zeros(n,1);
DCM.d    = zeros(n,n,0);

DCM.Y.y  = y + e;
DCM.Y.dt = TR;
DCM.U.u  = zeros(T,1);
DCM.U.dt = TR;



% provisional inversion
%--------------------------------------------------------------------------
DCM   = spm_dcm_fmri_csd(DCM);


% replace original connectivty with posterior expectations
%--------------------------------------------------------------------------
pP.A  = DCM.Ep.A;
M.g   = 'spm_gx_fmri';

% re-integrate states and observation noise process
%--------------------------------------------------------------------------
U.u   = spm_rand_mar(T,n,1/2)/2;      % endogenous fluctuations
y     = spm_int_J(pP,M,U);            % integrate with observer
e     = spm_rand_mar(T,n,1/2)/8;

% response
% -----------------------------------------------------------------
DCM.Y.y  = y + e;

% nonlinear system identification (Variational Laplace)
% =================================================================
DCM   = spm_dcm_fmri_csd_DEM(DCM);


% summary
% -----------------------------------------------------------------
spm_figure('Getwin','Figure 2'); clf

subplot(2,1,1); hold off
j   = 1:numel(DCM.Ep.A);
spm_plot_ci(DCM.Ep.A(:),DCM.Cp(j,j)), hold on
bar(pP.A(:),1/4), hold off
title('True and MAP connections','FontSize',16)
axis square

subplot(2,2,3); cla
plot(spm_vec(pP.A),spm_vec(DCM.Ep.A),'.','MarkerSize',16), hold off
title('True and MAP connections (Extrinsic)','FontSize',16)
xlabel('True')
ylabel('Estimated')
axis square



return

% NOTES: illustrate the ill-posed nature of the problem
%==========================================================================
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
        Y(:,end + 1) = spm_vec(spm_csd_fmri_mtf(pp,DCM.M,DCM.U));
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

spm_figure('Getwin','Figure 2'); clf
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


