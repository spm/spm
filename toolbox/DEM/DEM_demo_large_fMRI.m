function DEM_demo_large_fMRI
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
% $Id: DEM_demo_large_fMRI.m 5691 2013-10-11 16:53:00Z karl $

% Simulate timeseries
%==========================================================================
rng('default')

% DEM Structure: create random inputs
% -------------------------------------------------------------------------
N  = 8;                               % number of runs
T  = 512;                             % number of observations (scans)
TR = 2;                               % repetition time or timing
n  = 8;                               % number of regions or nodes
t  = (1:T)*TR;                        % observation times

% priors
% -------------------------------------------------------------------------
options.nmax       = 8;               % effective number of notes

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


% true parameters (reciprocal connectivity)
% -------------------------------------------------------------------------
pP.A = randn(n,n)/8;
pP.A = pP.A - diag(diag(pP.A));
pP.C = eye(n,n);
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
CSD   = {};
RMS   = [];
Qp    = [];
Pp    = pP.A - diag(diag(pP.A));
for i = 1:N
    try

        % integrate states
        % -----------------------------------------------------------------
        U.u  = spm_rand_mar(T,n,1/2)/2;      % endogenous fluctuations
        y    = spm_int_J(pP,M,U);            % integrate with observer
        
        % observation noise process
        % -----------------------------------------------------------------
        e    = spm_rand_mar(T,n,1/2)/8;
        
        % response
        % -----------------------------------------------------------------
        DCM.Y.y  = y + e;
        
        % nonlinear system identification (Variational Laplace)
        % =================================================================
        CSD{end + 1} = spm_dcm_fmri_csd(DCM);
        BPA    = spm_dcm_average(CSD,'simulation',1);
        
        % MAP estimates
        % -----------------------------------------------------------------
        qp     = CSD{end}.Ep.A;
        Qp(:,end + 1) = spm_vec(qp - diag(diag(qp)));
                
        % root mean square error
        % -----------------------------------------------------------------
        dp     = BPA.Ep.A - pP.A;
        dp     = dp - diag(diag(dp));
        RMS(end + 1) = sqrt(mean(dp(~~dp).^2));
        
        
        % summary
        % -----------------------------------------------------------------
        spm_figure('Getwin','Figure 2'); clf
        
        subplot(2,1,1); hold off
        spm_plot_ci(BPA.Ep.A(:),BPA.Cp(1:n*n,1:n*n)), hold on
        bar(pP.A(:),1/4), hold off
        title('True and MAP connections','FontSize',16)
        axis square
        
        subplot(2,2,3); cla
        plot(Pp(:),Qp,'cd','MarkerSize',8),hold on
        plot(Pp,BPA.Ep.A - diag(diag(BPA.Ep.A)),'.','MarkerSize',16), hold off
        title('True and MAP connections (Extrinsic)','FontSize',16)
        xlabel('True')
        ylabel('Estimated')
        axis square
        
        subplot(2,2,4);
        plot(RMS)
        title('root mean square error','FontSize',16)
        xlabel('number of sessions')
        ylabel('RMS')
        axis square
        drawnow
        
    end
end
