function spm_erp2csd_demo
% Demo routine for local field potential models
%==========================================================================
%
% This routine illustrates the use of empirical Bayes, for dynamic causal
% modelling, in identifying the causes of paroxysmal seizure activity; as
% expressed in terms of spectral density responses. We first simulate data
% by generating (endogenous neuronal) inputs under a scale free or power
% law assumption (the priors used for DCM for CSD). The inputs are used to
% generate responses over two seconds, whose spectral density is then used
% to estimate the neural mass model parameters. This is repeated for
% several different values of a particular intrinsic connection strength.
% Empirical Bayes is then used to compare competing models of between
% epoch changes in intrinsic connections. The  posterior distributions
% are then compared with the true values, under the selected model.
%
% The key aspects of this demonstration are to show that cross spectral
% density data features can be used to summarise evoked responses – and
% that trial to trial (or condition to condition) variations in model
% parameters can be identified using model selection, under a parametric
% random effect or empirical Bayesian model, which furnishes posterior
% densities over parameters at the first or within trial Level.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_erp2csd_demo.m 6759 2016-03-27 19:45:17Z karl $


% Model specification
%==========================================================================
rng('default')

% number of regions in simulated seizure activity and model specification
%--------------------------------------------------------------------------
rG  = [0,2];                               % range of parameter
Nc  = 1;                                   % number of channels
Ns  = 1;                                   % number of sources
options.spatial  = 'LFP';
options.model    = 'CMC';
options.analysis = 'CSD';
M.dipfit.model   = options.model;
M.dipfit.type    = options.spatial;
M.dipfit.Nc      = Nc;
M.dipfit.Ns      = Ns;
M.Hz             = 1:128;

% get associated priors
%--------------------------------------------------------------------------
A      = {0 0 0};
B      = {};
C      = 1;
pE     = spm_dcm_neural_priors(A,B,C,options.model);
pE     = spm_L_priors(M.dipfit,pE);
pE     = spm_ssr_priors(pE);
[x,f]  = spm_dcm_x_neural(pE,options.model);

% suppress channel noise (assuming many trials would be averaged)
%--------------------------------------------------------------------------
pE.a   = [ 0; 0];                  % log amplitude and f^(-a) exponent
pE.b   = [-2; 0];                  % log amplitude and f^(-a) exponent
pE.c   = [-2; 0];                  % log amplitude and f^(-a) exponent

% number of hidden states and endogenous inputs
%--------------------------------------------------------------------------
nx     = length(spm_vec(x));
nu     = size(pE.C,2);

% create LFP model
%==========================================================================
M.f    = f;
M.g    = 'spm_gx_erp';
M.x    = x;
M.n    = nx;
M.pE   = pE;
M.m    = nu;
M.l    = Nc;

% Volterra Kernels and transfer functions
%==========================================================================
%  illustrate response characteristics in terms of Volterra Kernels and
%  transfer functions. The idea here is to illustrate the dependency on a
%  particular connection parameter; here the inhibitory connection between
%  superficial pyramidal cells and spiny stellate cells.
%--------------------------------------------------------------------------
spm_figure('GetWin','Volterra kernels and transfer functions');

% augment and bi-linearise (with delays)
%--------------------------------------------------------------------------
[f,J,D]       = spm_fx_cmc(x,0,pE,M);
M.u           = sparse(Ns,1);
M.D           = D;
[M0,M1,L1,L2] = spm_bireduce(M,pE);

% compute kernels (over 64 ms)
%--------------------------------------------------------------------------
N          = 64;
dt         = 1/1000;
pst        = (1:N)*dt*1000;
[K0,K1,K2] = spm_kernels(M0,M1,L1,L2,N,dt);

subplot(2,2,1)
plot(pst,K1(:,:,1))
title('1st-order Volterra kernel','FontSize',16)
axis square
xlabel('time (ms)')

subplot(2,2,2)
imagesc(pst,pst,K2(1:64,1:64,1,1,1))
title('2nd-order Volterra kernel','FontSize',16)
axis square
xlabel('time (ms)')

% compute transfer functions for different inhibitory connections
%--------------------------------------------------------------------------
p     = linspace(-1,3,64);
for i = 1:length(p)
    P       = pE;
    P.G(3)  = p(i);
    [G,w]   = spm_csd_mtf(P,M);
    GW(:,i) = abs(G{1});
end

subplot(2,2,3)
plot(w,GW)
xlabel('frequency {Hz}')
title('transfer function','FontSize',16)
axis square

subplot(2,2,4)
imagesc(p,w,log(GW)), hold on
plot([rG(1),rG(1)],[w(1),w(end)],':w'), hold on
plot([rG(2),rG(2)],[w(1),w(end)],':w'), hold off
title('transfer functions','FontSize',16)
ylabel('Frequency')
xlabel('Inhibitory connection','FontSize',16)
axis xy
axis square


% illustrate responses to random fluctuations
%==========================================================================
spm_figure('GetWin','spontaneous fluctuations');clf

U.dt  = 4/1000;
N     = 2/U.dt;
pst   = (1:N)'*U.dt;


% spectral densities  of inputs and noise (under prior expectations)
%--------------------------------------------------------------------------
[Gu,Gs] = spm_csd_mtf_gu(pE,M.Hz);
[n,Kn]  = spm_rand_power_law(Gs,M.Hz,U.dt,N);
[u,Ku]  = spm_rand_power_law(Gu,M.Hz,U.dt,N);

% generate data under increasing values (G) of the inhibitory connection:
%--------------------------------------------------------------------------
nt    = 8;
M.p   = 8;
M.dt  = U.dt;
G     = linspace(rG(1),rG(2),4);
for i = 1:length(G)
    
    % change intrisic connectivity
    %----------------------------------------------------------------------
    P      = pE;
    P.G(3) = G(i);
    
    % predicted spectral density
    %----------------------------------------------------------------------
    psd(i) = spm_csd_mtf(P,M,U);
    
    % sample density over (nt) trials
    %----------------------------------------------------------------------
    csd{i} = 0;
    for j = 1:nt
        
        % generate random neuronal fluctuations - u (and noise - n)
        %------------------------------------------------------------------
        n    = randn(N,1);
        n    = Kn*n;
        u    = randn(N,1); u(fix(N/8)) = 16;
        u    = Ku*u;
        
        
        % generated evoked response and estimate induced response (with noise)
        %------------------------------------------------------------------
        U.u    = u;
        y      = spm_int_L(P,M,U);
        mar    = spm_mar(y + n,M.p);
        mar    = spm_mar_spectra(mar,M.Hz,1/U.dt);
        csd{i} = csd{i} + mar.P/nt;
        
    end
    
    % plot
    %----------------------------------------------------------------------
    col = [1/2 1 1]*(1 - i/length(G));
    t   = find(pst < .5);
    
    if i == 1
        
        subplot(3,2,1), plot(pst(t),u(t),'Color',col), spm_axis tight
        xlabel('time (s)'), title(' neuronal input','FontSize',16)
        
        % LFP – random fluctuations
        %----------------------------------------------------------------------
        subplot(3,2,2), plot(pst(t),y(t),'Color',col), hold on
        plot(pst(t),n(t),':','Color',col), hold off, spm_axis tight
        xlabel('time (s)'), title('LFP response (low connectivity)','FontSize',16)
        
        
    elseif i == length(G)
        
        subplot(3,2,3), plot(pst(t),u(t),'Color',col), spm_axis tight
        xlabel('time (s)'), title(' neuronal input','FontSize',16)
        
        % LFP – random fluctuations
        %----------------------------------------------------------------------
        subplot(3,2,4), plot(pst(t),y(t),'Color',col), hold on
        plot(pst(t),n(t),':','Color',col), hold on, spm_axis tight
        xlabel('time (s)'), title('(high connectivity)','FontSize',16)
        
    end
    
    
    % plot spectral density with changes in connectivity
    %----------------------------------------------------------------------
    subplot(3,2,5), plot(M.Hz,abs(csd{i}),'Color',col), hold on
    xlabel('time (s)'), title('spectral density','FontSize',16)
    spm_axis tight
    
    subplot(3,2,6), plot(M.Hz,abs(psd{i}),'Color',col), hold on
    xlabel('time (s)'), title(' predicted','FontSize',16)
    spm_axis tight
    drawnow
    
    
end


%  inverse spectral responses for multiple levels of the parameter
%==========================================================================

%  setup model (DCM structure
%--------------------------------------------------------------------------
options.DATA = 0;
DCM.A        = A;
DCM.B        = B;
DCM.C        = C;
DCM.options  = options;
DCM.M.dipfit = M.dipfit;
DCM.M.p      = M.p;
DCM.M.dt     = U.dt;

DCM.xU       = U;

DCM.xY.Hz    = M.Hz;
DCM.xY.dt    = U.dt;
DCM.xY.pst   = pst*1000;

% invert
%--------------------------------------------------------------------------
for i = 1:length(csd);
    DCM.xY.y  = csd(i);
    CSD{i,1}  = DCM;
end

CSD = spm_dcm_fit(CSD);

% Parametric empirical Bayes
%==========================================================================
%  first identify the most likely parameter explaining changes. This
%  involves  comparing models defined in terms of the relative within and
%  between trial variance. In other words, isolate a single parameter by
%  setting the precision of between trial variance of the remaining
%  parameters to a relatively small value. The free energy of the
%  empirical Bayes model can then be used to evaluate which parameter
%  estimates are likely to be responsible for between trial variations.
%--------------------------------------------------------------------------
clear M;
M.hE  = 0;
M.hC  = 1/16;
M.bE  =      spm_vec(CSD{1}.M.pE);
M.bC  = diag(spm_vec(CSD{1}.M.pC));


%  Define model space in terms of parameters allowed to vary
%--------------------------------------------------------------------------
field = {'T(1)','T(2)','T(3)','G(1)','G(2)','G(3)'};
beta  = 32;
for i = 1:length(field)
    
    % restrict between trial precision (beta)
    %----------------------------------------------------------------------
    j         = spm_find_pC(CSD{1},field{i});
    M.pC      = M.bC/beta;
    M.pC(j,j) = M.bC(j,j);
    
    % invert and record free energy
    %----------------------------------------------------------------------
    PEB    = spm_dcm_peb(CSD,M,'all');
    F(i,1) = PEB.F;
    
end

% results in terms of model comparison
%--------------------------------------------------------------------------
spm_figure('GetWin','empirical Bayesian results');clf

subplot(2,1,1)
bar(spm_softmax(F));set(gca,'XtickLabel',field)
xlabel('model'), ylabel('probability')
title('model comparison','FontSize',16), axis square


% recover posterior density and the selected model
%--------------------------------------------------------------------------
j         = spm_find_pC(CSD{1},field{end});
M.pC      = M.bC/beta;
M.pC(j,j) = M.bC(j,j);
[PEB,DCM] = spm_dcm_peb(CSD,M,'all');

for i = 1:length(csd);
    q  = spm_vec(DCM{i}.Ep); qE(i,1)  = q(j);
    q  = diag(   DCM{i}.Cp); qC(i,1)  = q(j);
end

subplot(2,1,2)
spm_plot_ci(qE,qC), hold on, bar(G,1/4), hold off
xlabel('trial'), ylabel('parameter estimate')
title('true and estimated changes','FontSize',16), axis square


return


