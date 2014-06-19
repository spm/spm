function spm_dcm_Granger_demo
% Demo routine for induced responses
%==========================================================================
%
% This routine illustrates the relationship between Geweke Granger 
% causality (GC) in frequency space and modulation transfer functions 
% (MTF).  We first compare and contrast analytic results for GC with 
% estimates based on a simulated time series. These synthetic data are 
% chosen to show that (analytic) GC can, in principle, detect sparsity 
% structure in terms of missing causal connections (however, GC estimates 
% are not so efficient). We then demonstrate the behaviour of (analytic) 
% GC by varying the strength of forward connections, backward connections 
% and intrinsic gain.  There is reasonable behaviour under these 
% manipulations. However, when we introduce realistic levels of (power law) 
% measurement noise, GC fails. The simulations conclude by showing that DCM 
% recovery of the underlying model parameters can furnish  (analytic) GC 
% among sources (in the absence of measurement noise). [delete the 'return'
% below to see these simulations].
% 
% See also:
%  spm_ccf2csd.m, spm_ccf2mar, spm_csd2ccf.m, spm_csd2mar.m, spm_mar2csd.m, 
%  spm_csd2coh.m, spm_ccf2gew, spm_dcm_mtf.m, spm_Q.m, spm_mar.m and 
%  spm_mar_spectral.m
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_Granger_demo.m 6059 2014-06-19 11:57:31Z vladimir $
 
 
% Model specification
%==========================================================================
rng('default')
 
% number of regions
%--------------------------------------------------------------------------
Nc    = 2;                                       % number of channels
Ns    = 2;                                       % number of sources
ns    = 2*128;                                   % sampling frequency
dt    = 1/ns;                                    % time bins
Hz    = 1:128;                                   % frequency
p     = 16;                                      % autoregression order
options.spatial  = 'LFP';
options.model    = 'CMC';
options.analysis = 'CSD';
M.dipfit.model = options.model;
M.dipfit.type  = options.spatial;
M.dipfit.Nc    = Nc;
M.dipfit.Ns    = Ns;
M.pF.D         = [1 4];
 
% extrinsic connections (forward an backward)
%--------------------------------------------------------------------------
A{1} = [0 0; 0 0];
A{2} = [0 0; 0 0];
A{3} = [0 0; 0 0];
B    = {};
C    = sparse(2,0);
 
% get priors
%--------------------------------------------------------------------------
pE    = spm_dcm_neural_priors(A,B,C,options.model);
pE    = spm_L_priors(M.dipfit,pE);
pE    = spm_ssr_priors(pE);
x     = spm_dcm_x_neural(pE,options.model);

% (log) connectivity parameters
%--------------------------------------------------------------------------
pE.A{1}(2,1) = 2;
pE.S         = 1/8;

% (log) amplitude of fluctations and noise
%--------------------------------------------------------------------------
pE.a(1,:) = -2;
pE.b(1,:) = -8;
pE.c(1,:) = -8;

 
% orders and model
%==========================================================================
nx    = length(spm_vec(x));
 
% create forward model
%--------------------------------------------------------------------------
M.f   = 'spm_fx_cmc';
M.g   = 'spm_gx_erp';
M.x   = x;
M.n   = nx;
M.pE  = pE;
M.m   = Ns;
M.l   = Nc;
M.Hz  = Hz;
M.Rft = 4;


% specify M.u - endogenous input (fluctuations) and intial states
%--------------------------------------------------------------------------
M.u   = sparse(Ns,1);
 
% solve for steady state
%--------------------------------------------------------------------------
M.x   = spm_dcm_neural_x(pE,M);


% Analytic spectral chararacterisation
%==========================================================================
spm_figure('GetWin','Figure 1'); clf

[csd,Hz,mtf] = spm_csd_mtf(pE,M);
csd          = csd{1};
mtf          = mtf{1};
ccf          = spm_csd2ccf(csd,Hz,dt);
mar          = spm_ccf2mar(ccf,p);
mar          = spm_mar_spectra(mar,Hz,ns);

spm_figure('GetWin','Figure 1'); clf
spm_spectral_plot(Hz,csd,  'b', 'frequency','density')
spm_spectral_plot(Hz,mar.P,'r', 'frequency','density')

legend('cross spectral density',...
       'autoregressive model')


%  comparison of expected results
%==========================================================================
spm_figure('GetWin','Figure 2'); clf

dtf  = mar.dtf;
gew  = mar.gew;
new  = spm_csd2gew(csd,Hz);

spm_spectral_plot(Hz,mtf,'b', 'frequency','density')
spm_spectral_plot(Hz,dtf,'g', 'frequency','density')
spm_spectral_plot(Hz,gew,'r', 'frequency','density')
spm_spectral_plot(Hz,new,'r:','frequency','density')


subplot(2,2,3), a = axis; subplot(2,2,2), axis(a);

legend('modulation transfer function',...
       'directed transfer function',...
       'Granger causality (parametric)',...
       'Granger causality (non-parametric)')


% effects of changing various model parameters
%==========================================================================

% (log) scaling, and parameters
%--------------------------------------------------------------------------
logs  = [ ((1:4)/1 - 2);
          ((1:4)/1 - 2);
          ((1:4)/6 + 0);
          ((1:4)/1 - 6)];

param = {'A{1}(2,1)','A{3}(1,2)','S','c(1,2)'};
str   = {     'forward connectivity',
              'backward connectivity',
              'intrinsic gain',
              'amplitude of noise'};


% expected transfer function and Gramger causality
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
a     = [0 Hz(end) 0 .35];
ca    = 0;
for i = 1:size(logs,1)
    for j = 1:size(logs,2)
        
        P     = pE;
        eval(sprintf('P.%s = %i;',param{i},logs(i,j)));
        
        % create forward model and solve for steady state
        %------------------------------------------------------------------
        M.x   = spm_dcm_neural_x(P,M);
        
        % Analytic spectral chararacterisation
        %==================================================================
        [csd,Hz,mtf,qew] = spm_csd_mtf(P,M);
        
        ccf   = spm_csd2ccf(csd{1},Hz,dt);
        mar   = spm_ccf2mar(ccf,p);
        mar   = spm_mar_spectra(mar,Hz,ns);
        
        
        % plot forwards and backwards functions
        %------------------------------------------------------------------
        subplot(4,2,ca + 1)
        plot(Hz,abs(mar.gew(:,2,1)),Hz,abs(qew{1}(:,2,1)))
        xlabel('frequency')
        ylabel('absolute value')
        title(sprintf('%s',str{i}),'FontSize',16)
        axis square, hold on, set(gca,'XLim',[0 Hz(end)])
        axis(a);
        
        subplot(4,2,ca + 2)
        plot(Hz,abs(mar.gew(:,1,2)),Hz,abs(qew{1}(:,1,2)))
        xlabel('frequency')
        ylabel('absolute value')
        title(sprintf('backward'),'FontSize',16)
        axis square, hold on, set(gca,'XLim',[0 Hz(end)])
        axis(a);

    end
    ca  = ca + 2;

end


% Lyapunov functions
%==========================================================================

% (log) scaling, and parameters
%--------------------------------------------------------------------------
param = {'A{1}(2,1)','A{3}(1,2)','S','D(1,1)'};
str   = {     'forward connectivity',
              'backward connectivity',
              'intrinsic gain',
              'intrinsic delays'};

% changing parameters and Lyapunov functions
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4'); clf
k     = linspace(0,1,32);
for i = 1:length(param)
    for j = 1:length(k)
        
        P     = pE;
        eval(sprintf('P.%s = %i;',param{i},k(j)));
        
        % create forward model and solve for steady state
        %------------------------------------------------------------------
        M.x   = spm_dcm_neural_x(P,M);
        
        % Analytic spectral chararacterisation
        %==================================================================
        [f,dfdx,D] = feval(M.f,M.x,M.u,P,M);
        dfdx       = D*dfdx;
        S(i,j)     = max(real(eig(full(dfdx),'nobalance')));
        
        % recursive check on MAR approximation
        %------------------------------------------------------------------
        csd        = spm_csd_mtf(P,M);
        ccf        = spm_csd2ccf(csd{1},Hz);
        [mar,c]    = spm_ccf2mar(ccf,p);
        pcond(i,j) = c;
        
    end
end

% plot
%--------------------------------------------------------------------------
subplot(2,2,1)
plot(k,S),hold on
xlabel('parameter (log) scaling')
ylabel('principal eigenvalue')
title('Lyapunov exponents','FontSize',16)
axis square; spm_axis tight
plot(k,k*0,'k:',k,k*0 - 1/(dt*p),'k-.'),hold off

subplot(2,2,2)
semilogy(k,pcond),hold on
xlabel('parameter (log) scaling')
ylabel('condition number')
title('condition of cross covariance','FontSize',16)
axis square;
legend(str)

plot(k,k*0 + 1e4,'k:'),hold off

% MAR approximation
%==========================================================================
param = {'S'};
for i = 1:length(param)
    for j = 1:length(k)
        
        P     = pE;
        eval(sprintf('P.%s = %i;',param{i},k(j)));
        
        % create forward model and solve for steady state
        %------------------------------------------------------------------
        M.x   = spm_dcm_neural_x(P,M);
        
        % check on MAR approximation
        %------------------------------------------------------------------
        csd   = spm_csd_mtf(P,M);
        mar   = spm_csd2mar(csd{1},Hz,p);
        mar   = spm_mar_spectra(mar,Hz);
        
        CSD(:,j) = csd{1}(:,2,2);
        PSD(:,j) = mar.P(:,2,2);
        
    end
end

subplot(2,2,3)
imagesc(k,Hz,log(abs(CSD)))
xlabel('parameter (log) scaling')
ylabel('frequency')
title('(log) spectral density','FontSize',16)
axis square

subplot(2,2,4)
imagesc(k,Hz,abs(PSD - CSD).^(1/4))
xlabel('parameter (log) scaling')
ylabel('frequency')
title('MAR approximation error','FontSize',16)
axis square

% a more careful examination of [in]stability on Granger causality
%==========================================================================
spm_figure('GetWin','Figure 5'); clf
k     = linspace(0,2/3,8);
for j = 1:length(k)
    
    
    % amplitude of observation noise
    %----------------------------------------------------------------------
    P        = pE;
    P.S      = k(j);
       
    % create forward model and solve for steady state
    %----------------------------------------------------------------------
    M.x      = spm_dcm_neural_x(P,M);
    
    % Analytic spectral chararacterisation (parametric)
    %======================================================================
    [csd,Hz,mtf,qew] = spm_csd_mtf(P,M);
    
    ccf      = spm_csd2ccf(csd{1},Hz,dt);
    mar      = spm_ccf2mar(ccf,p);
    mar      = spm_mar_spectra(mar,Hz,ns);
    
    % and non-parametric
    %======================================================================
    gew      = spm_csd2gew(csd{1},Hz);
    
    % save forwards and backwards functions
    %----------------------------------------------------------------------
    GCF(:,j) = abs(mar.gew(:,2,1));
    GCB(:,j) = abs(mar.gew(:,1,2));
    
    % plot forwards and backwards functions (parametric)
    %----------------------------------------------------------------------
    spm_figure('GetWin','Figure 5');
    
    subplot(3,2,1)
    plot(Hz,abs(mar.gew(:,2,1)),Hz,abs(qew{1}(:,2,1)))
    xlabel('frequency')
    ylabel('absolute value')
    title('forward','FontSize',16)
    axis square, hold on, set(gca,'XLim',[0 Hz(end)])
    a  = axis;
    
    subplot(3,2,2)
    plot(Hz,abs(mar.gew(:,1,2)),Hz,abs(qew{1}(:,1,2)))
    xlabel('frequency')
    ylabel('absolute value')
    title('backward','FontSize',16)
    axis square, hold on, set(gca,'XLim',[0 Hz(end)])
    axis(a);

    
    % plot forwards and backwards functions (nonparametric)
    %----------------------------------------------------------------------
    subplot(3,2,3)
    plot(Hz,abs(gew(:,2,1)),Hz,abs(qew{1}(:,2,1)))
    xlabel('frequency')
    ylabel('absolute value')
    title('forward','FontSize',16)
    axis square, hold on, set(gca,'XLim',[0 Hz(end)])
    axis(a);
    
    subplot(3,2,4)
    plot(Hz,abs(gew(:,1,2)),Hz,abs(qew{1}(:,1,2)))
    xlabel('frequency')
    ylabel('absolute value')
    title('backward','FontSize',16)
    axis square, hold on, set(gca,'XLim',[0 Hz(end)])
    axis(a);
    
end

a = .2;
subplot(3,2,5)
image(Hz,k,GCF'*64/a)
xlabel('frequency')
ylabel('log(exponent)')
title('forward','FontSize',16)
axis square

subplot(3,2,6)
image(Hz,k,GCB'*64/a)
xlabel('frequency')
ylabel('log(exponent)')
title('backward','FontSize',16)
axis square

% a more careful examination of measurement noise
%==========================================================================
spm_figure('GetWin','Figure 6' ); clf
spm_figure('GetWin','Figure 6a'); clf
k     = linspace(-8,-0,8);
for j = 1:length(k)
    
    
    % amplitude of observation noise
    %----------------------------------------------------------------------
    P        = pE;
    P.b(1)   = k(j);
    P.c(1,:) = k(j) + 1;
       
    % create forward model and solve for steady state
    %----------------------------------------------------------------------
    M.x      = spm_dcm_neural_x(P,M);
    
    % Analytic spectral chararacterisation (parametric)
    %======================================================================
    [csd,Hz,mtf,qew] = spm_csd_mtf(P,M);
    
    ccf      = spm_csd2ccf(csd{1},Hz,dt);
    mar      = spm_ccf2mar(ccf,p);
    mar      = spm_mar_spectra(mar,Hz,ns);
    
    % and non-parametric
    %======================================================================
    gew      = spm_csd2gew(csd{1},Hz);
    qew      = qew{1};
    
    % save forwards and backwards functions
    %----------------------------------------------------------------------
    GCF(:,j) = abs(gew(:,2,1));
    GCB(:,j) = abs(gew(:,1,2));
    
    % plot forwards and backwards functions  (parametric)
    %----------------------------------------------------------------------
    spm_figure('GetWin','Figure 6');
    
    subplot(3,2,1)
    plot(Hz,abs(mar.gew(:,2,1)),Hz,abs(qew(:,2,1)))
    xlabel('frequency')
    ylabel('absolute value')
    title('forward','FontSize',16)
    axis square, hold on, set(gca,'XLim',[0 Hz(end)])
    a  = axis;
    
    subplot(3,2,2)
    plot(Hz,abs(mar.gew(:,1,2)),Hz,abs(qew(:,1,2)))
    xlabel('frequency')
    ylabel('absolute value')
    title('backward','FontSize',16)
    axis square, hold on, set(gca,'XLim',[0 Hz(end)])
    axis(a);

    
    % plot forwards and backwards functions  (nonparametric)
    %----------------------------------------------------------------------
    subplot(3,2,3)
    plot(Hz,abs(gew(:,2,1)),Hz,abs(qew(:,2,1)))
    xlabel('frequency')
    ylabel('absolute value')
    title('forward','FontSize',16)
    axis square, hold on, set(gca,'XLim',[0 Hz(end)])
    axis(a);
    
    subplot(3,2,4)
    plot(Hz,abs(gew(:,1,2)),Hz,abs(qew(:,1,2)))
    xlabel('frequency')
    ylabel('absolute value')
    title('backward','FontSize',16)
    axis square, hold on, set(gca,'XLim',[0 Hz(end)])
    axis(a);
    
    % plot associated coherence
    %----------------------------------------------------------------------
    spm_figure('GetWin','Figure 6a' );
    
    if j > 4, str = ':'; else, str = '-'; end
    coh = spm_csd2coh(csd{1},Hz);
    spm_spectral_plot(Hz,coh,str,'frequency','coherence')
    
end

spm_figure('GetWin','Figure 6' );

subplot(3,2,5)
imagesc(Hz,k,GCF')
xlabel('frequency')
ylabel('log(exponent)')
title('forward','FontSize',16)
axis square

subplot(3,2,6)
imagesc(Hz,k,GCB')
xlabel('frequency')
ylabel('log(exponent)')
title('backward','FontSize',16)
axis square



% return if in demonstration mode
%--------------------------------------------------------------------------
DEMO = 1;
if DEMO, return, end


% DCM estimates of coupling
%==========================================================================
rng('default')

% get priors and generate data
%--------------------------------------------------------------------------
pE    = spm_dcm_neural_priors(A,B,C,options.model);
pE    = spm_L_priors(M.dipfit,pE);
pE    = spm_ssr_priors(pE);

% (log) connectivity parameters (forward connection only)
%--------------------------------------------------------------------------
pE.A{1}(2,1) = 2;
pE.S         = 1/8;

% (log) amplitude of fluctations and noise
%--------------------------------------------------------------------------
pE.a(1,:) = -2;
pE.b(1,:) = -8;
pE.c(1,:) = -8;


% expected cross spectral density
%--------------------------------------------------------------------------
csd       = spm_csd_mtf(pE,M);

% Get spectral profile of fluctuations and noise
%--------------------------------------------------------------------------
[Gu,Gs,Gn] = spm_csd_mtf_gu(pE,Hz);

% Integrate with power law process (simulate multiple trials)
%--------------------------------------------------------------------------
PSD   = 0;
CSD   = 0;
N     = 1024;
U.dt  = dt;
for t = 1:16
    
    % neuronal fluctuations
    %----------------------------------------------------------------------
    U.u      = spm_rand_power_law(Gu,Hz,dt,N);
    LFP      = spm_int_L(pE,M,U);
    
    % and measurement noise
    %----------------------------------------------------------------------
    En       = spm_rand_power_law(Gn,Hz,dt,N);
    Es       = spm_rand_power_law(Gs,Hz,dt,N);
    E        = Es + En*ones(1,Ns);
    
    % and estimate spectral features under a MAR model
    %----------------------------------------------------------------------
    MAR      = spm_mar(LFP + E,p);
    MAR      = spm_mar_spectra(MAR,Hz,ns);
    CSD      = CSD + MAR.P;
    
    % and using Welch's method
    %----------------------------------------------------------------------
    PSD      = PSD + spm_csd(LFP + E,Hz,ns);
    
    CCD(:,t) = abs(CSD(:,1,2)/t);
    PCD(:,t) = abs(CSD(:,1,1)/t);
   
    % plot
    %----------------------------------------------------------------------
    spm_figure('GetWin','Figure 7'); clf
    spm_spectral_plot(Hz,csd{1},'r',  'frequency','density')
    spm_spectral_plot(Hz,CSD/t, 'b',  'frequency','density')
    spm_spectral_plot(Hz,PSD/t, 'g',  'frequency','density')
    legend('real','estimated (AR)','estimated (CSD)')
    drawnow
    
end

%  show convergence of spectral estimators
%--------------------------------------------------------------------------
subplot(2,2,3), hold off
imagesc(Hz,1:t,log(PCD'))
xlabel('frequency')
ylabel('trial number')
title('log auto spectra','FontSize',16)
axis square

subplot(2,2,4), hold off
imagesc(Hz,1:t,log(CCD'))
xlabel('frequency')
ylabel('trial number')
title('log cross spectra','FontSize',16)
axis square

% DCM set up (allow for both forward and backward connections)
%==========================================================================

% (log) connectivity parameters (forward connection only)
%--------------------------------------------------------------------------
pE.A{1}(2,1) = 2;
pE.S         = 1/8;

% (log) amplitude of fluctations and noise
%--------------------------------------------------------------------------
pE.a(1,:) = -2;
pE.b(1,:) = -4;
pE.c(1,:) = -4;

% expected cross spectral density
%--------------------------------------------------------------------------
[csd,Hz,mtf] = spm_csd_mtf(pE,M);

DCM.options.model   = 'CMC';
DCM.options.spatial = 'LFP';
DCM.options.DATA    = 0;

DCM.M.dipfit.Nc  = Nc;
DCM.M.dipfit.Ns  = Ns;
DCM.M.Nmax = 32;
DCM.M.U    = eye(Nc,Nc);

DCM.A      = {[0 0; 1 0],[0 1; 0 0],[0 0; 0 0]};
DCM.B      = {};
DCM.C      = sparse(Ns,0);

% place in data structure
%--------------------------------------------------------------------------
DCM.xY.y  = csd;
DCM.xY.dt = dt;
DCM.xY.Hz = Hz;

% estimate
%--------------------------------------------------------------------------
DCM       = spm_dcm_csd(DCM);

% show results in terms of transfer functions and Granger causality
%==========================================================================
spm_figure('GetWin','Figure 8'); clf

% transfer functions in the absence of measurement noise
%--------------------------------------------------------------------------
Ep    = DCM.Ep; 
Ep.b  = Ep.b - 32;              % and suppress non-specific and
Ep.c  = Ep.c - 32;              % specific channel noise

Gu           = spm_csd_mtf_gu(Ep,Hz);
[psd Hz dtf] = spm_csd_mtf(Ep,DCM.M);

mtf   = mtf{1};
csd   = csd{1};
dtf   = dtf{1};

% transfer functions and Granger causality among sources and channels
%--------------------------------------------------------------------------
gwe  = spm_dtf2gew(mtf,Gu);
qew  = spm_dtf2gew(dtf,Gu);
gew  = spm_csd2gew(csd,Hz);

spm_spectral_plot(Hz,mtf,'b:', 'frequency','density')
spm_spectral_plot(Hz,gwe,'b',  'frequency','density')
spm_spectral_plot(Hz,qew,'r',  'frequency','density')
spm_spectral_plot(Hz,gew,'g',  'frequency','density')
legend('modulation transfer function',...
       'Granger causality (true)',...
       'Granger causality (source)',...
       'Granger causality (channel)')

subplot(2,2,3), a = axis; subplot(2,2,2), axis(a);

return


   
% NOTES: a more careful examination of delays
%==========================================================================
spm_figure('GetWin','Figure 9'); clf

k     = linspace(8,11,8);
for j = 1:length(k)
    
    
    % keep total power of fluctuations constant
    %----------------------------------------------------------------------
    M.pF.D(2,1)  = k(j);
       
    % create forward model and solve for steady state
    %----------------------------------------------------------------------
    M.x          = spm_dcm_neural_x(pE,M);
    
    % Analytic spectral chararacterisation (parametric)
    %======================================================================
    [csd,Hz,mtf] = spm_csd_mtf(pE,M);
    ccf          = spm_csd2ccf(csd{1},Hz,dt);
    mar          = spm_ccf2mar(ccf,p);
    mar          = spm_mar_spectra(mar,Hz,ns);
    
    % and non-parametric)
    %======================================================================
    gew          = spm_csd2gew(csd{1},Hz);
    
    
    % plot forwards and backwards functions
    %----------------------------------------------------------------------
    subplot(3,2,1)
    plot(Hz,abs(mar.gew(:,2,1)),Hz,abs(mtf{1}(:,2,1)))
    xlabel('frequency')
    ylabel('absolute value')
    title('forward','FontSize',16)
    axis square, hold on, set(gca,'XLim',[0 Hz(end)])
    a  = axis;
    
    subplot(3,2,2)
    plot(Hz,abs(mar.gew(:,1,2)),Hz,abs(mtf{1}(:,1,2)))
    xlabel('frequency')
    ylabel('absolute value')
    title('backward','FontSize',16)
    axis square, hold on, set(gca,'XLim',[0 Hz(end)])
    axis(a);

    
    % plot forwards and backwards functions
    %----------------------------------------------------------------------
    subplot(3,2,3)
    plot(Hz,abs(gew(:,2,1)),Hz,abs(mtf{1}(:,2,1)))
    xlabel('frequency')
    ylabel('absolute value')
    title('forward','FontSize',16)
    axis square, hold on, set(gca,'XLim',[0 Hz(end)])
    axis(a);
    
    subplot(3,2,4)
    plot(Hz,abs(gew(:,1,2)),Hz,abs(mtf{1}(:,1,2)))
    xlabel('frequency')
    ylabel('absolute value')
    title('backward','FontSize',16)
    axis square, hold on, set(gca,'XLim',[0 Hz(end)])
    axis(a);
    
end



