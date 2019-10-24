function spm_dcm_prior_responses(Ep)
% Demo routine that computes transfer functions for free parameters
%==========================================================================
%
% This routine provides a survey of responses under stationarity
% assumptions for the suite of neural mass and mean field models used in
% DCM. It characterises the steady-state responses - under prior
% expectations - using spectral density and autocovariance functions
% with and with out channel noise. it then proceeds to evaluate evoked
% responses to a canonical input.
%
% This function is used primarily to check the prior expectations to ensure
% the expected responses within a comparable and appropriate range for
% scale empirical data. The amplitude of the responses are set by the
% scaling of U in the equations of motion for each model.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_prior_responses.m 7679 2019-10-24 15:54:07Z spm $


% Model specification
%==========================================================================


% models
%--------------------------------------------------------------------------
model = {'ERP','SEP','LFP','CMC','TFM','CMM','NMM','MFM'};
model = model(:);
Nm    = length(model);

% spatial model and sources
%--------------------------------------------------------------------------
M.dipfit.type  = 'LFP';
M.dipfit.Nc    = 1;
M.dipfit.Ns    = 1;


% within-trial effects: adjust onset relative to pst
%----------------------------------------------------------------------
M.ns   = 64;
M.ons  = 64;
M.dur  = 16;
U.dt   = 0.004;
U.X    = [];

ifig  = 1;
for i = 1:Nm
    
    
    % dear figure and plot
    %======================================================================
    iplot = rem(i - 1,3);
    if iplot == 0
        spm_figure('GetWin',sprintf('stationary responses: %i',ifig));
        ifig = ifig + 1;
    end
    
    
    % set model option and priors
    %----------------------------------------------------------------------
    options.spatial = 'LFP';
    options.model   = model{i};
    M.dipfit.model  = options.model;
    
    [pE pC] = spm_dcm_neural_priors({0 0 0},{},1,options.model);
    [pE pC] = spm_L_priors(M.dipfit,pE,pC);
    [pE pC] = spm_ssr_priors(pE,pC);
    
    
    % hidden neuronal states of interest
    %----------------------------------------------------------------------
    [x,f]   = spm_dcm_x_neural(pE,options.model);
    
    % orders and model
    %======================================================================
    nx      = length(spm_vec(x));
    nu      = size(pE.C,2);
    u       = sparse(1,nu);
    
    % create LFP model
    %----------------------------------------------------------------------
    M.f     = f;
    M.g     = 'spm_gx_erp';
    M.x     = x;
    M.n     = nx;
    M.pE    = pE;
    M.pC    = pC;
    M.m     = nu;
    M.l     = 1;
    
    % solve for steady state
    %----------------------------------------------------------------------
    M.x     = spm_dcm_neural_x(pE,M);
    
    % evaluate stationary responses
    %======================================================================
    M.u     = u;
    M.Hz    = 4:64;
    
    
    % get spectral responses, auto covariance function and ERP
    %----------------------------------------------------------------------
    [csd, Hz ] = spm_csd_mtf(pE,M,[]);
    [ccf, lag] = spm_csd2ccf(csd,Hz);
    [erp, pst] = spm_gen_erp(pE,M,U);  
    
    
    % plot
    %----------------------------------------------------------------------
    lag = 1000*lag;
    pst = 1000*pst;
    
    subplot(3,3,3*iplot + 2)
    plot(Hz,abs(csd{1})),  hold on
    
    subplot(3,3,3*iplot + 3)
    plot(lag,ccf{1}), hold on
    
    
    % now repeat without channel noise
    %----------------------------------------------------------------------
    pE.b = pE.b - 32;
    pE.c = pE.c - 32;
    [csd, Hz]  = spm_csd_mtf(pE,M,[]);
    [ccf, lag] = spm_csd2ccf(csd,Hz);
    
    % plot
    %----------------------------------------------------------------------
    subplot(3,3,3*iplot + 2)
    plot(Hz,abs(csd{1}),'b--'), hold off
    xlabel('frequency {Hz}')
    title('autospectrum','FontSize',16)
    
    subplot(3,3,3*iplot + 3)
    plot(lag,ccf{1},'b--'), hold off
    title('autocovariance','FontSize',16)
    
    subplot(3,3,3*iplot + 1), hold off
    plot(pst,erp{1},'b'), hold on
    title(sprintf('ERP: %s',model{i}),'FontSize',16)
    
    % hidden states
    %----------------------------------------------------------------------
    M.g  = @(x,u,P,M) x;
    erp  = spm_gen_erp(pE,M,U);
    plot(pst,erp{1}/32,'r:'), hold off
    
end

