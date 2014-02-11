function spm_Granger_demo
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
% $Id: spm_Granger_demo.m 5877 2014-02-11 20:03:34Z karl $
 
 
% Model specification
%==========================================================================
rng('default')
 
% number of regions
%--------------------------------------------------------------------------
N     = 2048;                                    % unmber of time bins
Nc    = 2;                                       % number of channels
Ns    = 2;                                       % number of sources
Hz    = 1:128;                                   % frequency
dt    = 4/1000;                                  % time bins
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
[x,f] = spm_dcm_x_neural(pE,options.model);

% (log) connectivity parameters
%--------------------------------------------------------------------------
pE.A{1}(2,1) = 2;
pE.S         = 1/8;

% (log) amplitude of fluctations and noise
%--------------------------------------------------------------------------
pE.a(1,:) =  2;
pE.b(1,:) = -4;
pE.c(1,:) = -4;

 
% orders and model
%==========================================================================
nx    = length(spm_vec(x));
 
% create forward model
%--------------------------------------------------------------------------
M.f   = f;
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
[csd,Hz,mtf] = spm_csd_mtf(pE,M);
csd          = csd{1};
mtf          = mtf{1};
ccf          = spm_csd2ccf(csd,Hz,dt);
mar          = spm_ccf2mar(ccf,p);
mar          = spm_mar_spectra(mar,Hz,1/dt);


% Numerical spectral chararacterisation
%==========================================================================

% Get spectral profile of fluctuations and noise
%--------------------------------------------------------------------------
[Gu,Gs,Gn]   = spm_csd_mtf_gu(pE,M.Hz);


% Integrate with pink noise process
%--------------------------------------------------------------------------
U.dt  = dt;
U.u   = spm_rand_power_law(Gu,Hz,dt,N);
LFP   = spm_int_L(pE,M,U);

% and measurement noise
%--------------------------------------------------------------------------
En    = spm_rand_power_law(Gn,Hz,dt,N);
Es    = spm_rand_power_law(Gs,Hz,dt,N);
E     = Es + En*ones(1,Ns);

% and estimate spectral features under a MAR model
%--------------------------------------------------------------------------
MAR           = spm_mar(LFP + E,p);
MAR.noise_cov = mar.noise_cov;
MAR           = spm_mar_spectra(MAR,Hz,1/dt);

CSD = MAR.P;
CSP = mar.P;


% Graphics
%==========================================================================
spm_figure('GetWin','Figure 1'); clf

% compare analytic and numerical spctral densities
%--------------------------------------------------------------------------
subplot(2,2,1)
plot((1:512)*dt,LFP(1:512,:))
xlabel('time')
title('LFP','FontSize',16)
axis square, spm_axis tight

subplot(2,2,2)
plot([spm_vec(mar.lag),spm_vec(MAR.lag)])
xlabel('expected')
ylabel('estimated')
title('auto-regression coefficients','FontSize',16)
axis square, spm_axis tight

subplot(2,2,3)
plot(Hz,abs(csd(:,1,1)),Hz,abs(CSD(:,1,1)),'--',Hz,abs(CSP(:,1,1)),':')
xlabel('frequency')
ylabel('absolute value')
title('auto-spectral density','FontSize',16)
axis square, set(gca,'XLim',[0 Hz(end)])

subplot(2,2,4)
plot(Hz,abs(csd(:,2,1)),Hz,abs(CSD(:,2,1)),'--',Hz,abs(CSP(:,2,1)),':')
xlabel('frequency')
ylabel('absolute value')
title('cross-spectral density','FontSize',16)
axis square, set(gca,'XLim',[0 Hz(end)])
legend('expected','estimated','expected (MAR)')


%  comparison of expected and numerical results
%==========================================================================
spm_figure('GetWin','Figure 2'); clf

% compare analytic and numerical spctral densities
%--------------------------------------------------------------------------
dtf  = mar.dtf;
gew  = mar.gew;
DTF  = MAR.dtf;
GEW  = MAR.gew;
spm_spectral_plot(Hz,mtf,'b',  'frequency','density')
spm_spectral_plot(Hz,dtf,'g',  'frequency','density')
spm_spectral_plot(Hz,DTF,'g-.','frequency','density')
spm_spectral_plot(Hz,gew,'r',  'frequency','density')
spm_spectral_plot(Hz,GEW,'r-.','frequency','density')

subplot(2,2,3), a = axis; subplot(2,2,2), axis(a);

legend('modulation transfer function','directed transfer function','estimate','Granger causality','estimate')


% effects of changing various model parameters
%==========================================================================

% (log) connectivity parameters
%--------------------------------------------------------------------------
pE.A{1}(2,1) = 2;
pE.S         = 1/8;

% (log) amplitude of fluctations and noise
%--------------------------------------------------------------------------
pE.a(1,:) =  0;
pE.b(1,:) = -4;
pE.c(1,:) = -4;


% (log) scaling, and parameters
%--------------------------------------------------------------------------
logs  = [ ((1:4)/1 - 2);
          ((1:4)/1 - 2);
          ((1:4)/8 + 0);
          ((1:4)*2 - 8)];


param = {'A{1}(2,1)','A{3}(1,2)','S','b(1,:)'};
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
        [csd,Hz,mtf] = spm_csd_mtf(P,M);
        ccf          = spm_csd2ccf(csd{1},Hz,dt);
        mar          = spm_ccf2mar(ccf,p);
        mar          = spm_mar_spectra(mar,Hz,1/dt);
        
        
        % plot forwards and backwards functions
        %------------------------------------------------------------------
        subplot(4,2,ca + 1)
        plot(Hz,abs(mar.gew(:,2,1)),Hz,abs(mtf{1}(:,2,1)),'--')
        xlabel('frequency')
        ylabel('absolute value')
        title(sprintf('%s',str{i}),'FontSize',16)
        axis square, hold on, set(gca,'XLim',[0 Hz(end)])
        axis(a);
        
        subplot(4,2,ca + 2)
        plot(Hz,abs(mar.gew(:,1,2)),Hz,abs(mtf{1}(:,1,2)),'--')
        xlabel('frequency')
        ylabel('absolute value')
        title(sprintf('backward'),'FontSize',16)
        axis square, hold on, set(gca,'XLim',[0 Hz(end)])
        axis(a);

    end
    ca  = ca + 2;

end

% a more careful examination of fluctuations
%==========================================================================
spm_figure('GetWin','Figure 4'); clf
k     = linspace(-4,-1,8);
for j = 1:length(k)
    
    
    % amplitude of observation noise
    %----------------------------------------------------------------------
    P        = pE;
    P.b(1,:) = k(j);
       
    % create forward model and solve for steady state
    %----------------------------------------------------------------------
    M.x          = spm_dcm_neural_x(P,M);
    
    % Analytic spectral chararacterisation (parametric)
    %======================================================================
    [csd,Hz,mtf] = spm_csd_mtf(P,M);
    ccf          = spm_csd2ccf(csd{1},Hz,dt);
    mar          = spm_ccf2mar(ccf,p);
    mar          = spm_mar_spectra(mar,Hz,1/dt);
    
    % and non-parametric)
    %======================================================================
    gew          = spm_csd2gwe(csd{1},Hz);
    
    % save forwards and backwards functions
    %----------------------------------------------------------------------
    GCF(:,j)     = abs(gew(:,2,1));
    GCB(:,j)     = abs(gew(:,1,2));
    
    % plot forwards and backwards functions
    %----------------------------------------------------------------------
    subplot(3,2,1)
    plot(Hz,abs(mar.gew(:,2,1)),Hz,abs(mtf{1}(:,2,1)),'--')
    xlabel('frequency')
    ylabel('absolute value')
    title('forward','FontSize',16)
    axis square, hold on, set(gca,'XLim',[0 Hz(end)])
    a  = axis;
    
    subplot(3,2,2)
    plot(Hz,abs(mar.gew(:,1,2)),Hz,abs(mtf{1}(:,1,2)),'--')
    xlabel('frequency')
    ylabel('absolute value')
    title('backward','FontSize',16)
    axis square, hold on, set(gca,'XLim',[0 Hz(end)])
    axis(a);

    
    % plot forwards and backwards functions
    %----------------------------------------------------------------------
    subplot(3,2,3)
    plot(Hz,abs(gew(:,2,1)),Hz,abs(mtf{1}(:,2,1)),'--')
    xlabel('frequency')
    ylabel('absolute value')
    title('forward','FontSize',16)
    axis square, hold on, set(gca,'XLim',[0 Hz(end)])
    axis(a);
    
    subplot(3,2,4)
    plot(Hz,abs(gew(:,1,2)),Hz,abs(mtf{1}(:,1,2)),'--')
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


% a more careful examination of delays
%==========================================================================
spm_figure('GetWin','Figure 5'); clf

% (log) amplitude of fluctations and noise
%--------------------------------------------------------------------------
pE.a(1,:) =  0;
pE.b(1,:) = -4;
pE.c(1,:) = -4;


k     = linspace(8,12,8);
for j = 1:length(k)
    
    
    % keep total power of fluctuations constant
    %----------------------------------------------------------------------
    M.pF.D         = [1 k(j)];
       
    % create forward model and solve for steady state
    %----------------------------------------------------------------------
    M.x          = spm_dcm_neural_x(pE,M);
    
    % Analytic spectral chararacterisation (parametric)
    %======================================================================
    [csd,Hz,mtf] = spm_csd_mtf(pE,M);
    ccf          = spm_csd2ccf(csd{1},Hz,dt);
    mar          = spm_ccf2mar(ccf,p);
    mar          = spm_mar_spectra(mar,Hz,1/dt);
    
    % and non-parametric)
    %======================================================================
    gew          = spm_csd2gwe(csd{1},Hz);
    
    
    % plot forwards and backwards functions
    %----------------------------------------------------------------------
    subplot(3,2,1)
    plot(Hz,abs(mar.gew(:,2,1)),Hz,abs(mtf{1}(:,2,1)),'--')
    xlabel('frequency')
    ylabel('absolute value')
    title('forward','FontSize',16)
    axis square, hold on, set(gca,'XLim',[0 Hz(end)])
    a  = axis;
    
    subplot(3,2,2)
    plot(Hz,abs(mar.gew(:,1,2)),Hz,abs(mtf{1}(:,1,2)),'--')
    xlabel('frequency')
    ylabel('absolute value')
    title('backward','FontSize',16)
    axis square, hold on, set(gca,'XLim',[0 Hz(end)])
    axis(a);

    
    % plot forwards and backwards functions
    %----------------------------------------------------------------------
    subplot(3,2,3)
    plot(Hz,abs(gew(:,2,1)),Hz,abs(mtf{1}(:,2,1)),'--')
    xlabel('frequency')
    ylabel('absolute value')
    title('forward','FontSize',16)
    axis square, hold on, set(gca,'XLim',[0 Hz(end)])
    axis(a);
    
    subplot(3,2,4)
    plot(Hz,abs(gew(:,1,2)),Hz,abs(mtf{1}(:,1,2)),'--')
    xlabel('frequency')
    ylabel('absolute value')
    title('backward','FontSize',16)
    axis square, hold on, set(gca,'XLim',[0 Hz(end)])
    axis(a);
    
end



return


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

% Get spectral profile of fluctuations and noise
%--------------------------------------------------------------------------
[Gu,Gs,Gn] = spm_csd_mtf_gu(pE,M.Hz);

% Integrate with power law process
%--------------------------------------------------------------------------
U.dt  = dt;
U.u   = spm_rand_power_law(Gu,Hz,dt,N);
LFP   = spm_int_L(pE,M,U);

% and measurement noise
%--------------------------------------------------------------------------
En    = spm_rand_power_law(Gn,Hz,dt,N);
Es    = spm_rand_power_law(Gs,Hz,dt,N);
E     = Es + En*ones(1,Ns);

% and estimate spectral features under a MAR model
%--------------------------------------------------------------------------
MAR   = spm_mar(LFP + E,p);
MAR   = spm_mar_spectra(MAR,Hz,1/dt);
CSD   = MAR.P;


% DCM set up (allow for both forward and backward connections)
%==========================================================================
DCM.A     = {[0 0; 1 0],[0 1; 0 0],[0 0; 0 0]};
DCM.B     = {};
DCM.C     = sparse(Ns,0);
DCM.xU    = U;
DCM.xY.y  = {CSD};
DCM.xY.dt = dt;
DCM.xY.Hz = Hz;

DCM.options.model   = 'CMC';
DCM.options.spatial = 'LFP';

DCM.M.dipfit.Nc     = Nc;
DCM.M.dipfit.Ns     = Ns;

% estimate
%--------------------------------------------------------------------------
DCM.options.DATA = 0;
DCM       = spm_dcm_csd(DCM);

% show results in terms of transfer functions and Granger causality
%==========================================================================
spm_figure('GetWin','Figure 6'); clf

% transfer functions and Granger causality among sources
%--------------------------------------------------------------------------
mtf = DCM.dtf{1};
ccf = DCM.ccf{1};
gew = spm_ccf2gwe(ccf,Hz,dt,p);

spm_spectral_plot(Hz,mtf,'b',  'frequency','density')
spm_spectral_plot(Hz,gew,'r',  'frequency','density')

% Granger causality among channels
%--------------------------------------------------------------------------
P   = DCM.Ep;
csd = spm_csd_mtf(P,DCM.M);
ccf = spm_csd2ccf(csd{1},Hz,dt);
gew = spm_ccf2gwe(ccf,Hz,dt,p);

spm_spectral_plot(Hz,gew,'g--',  'frequency','density')
legend('modulation transfer function','Granger causality (source)','Granger causality (channel)')


