function [HMM,csd] = spm_dcm_HMM(GCM,N,b)
% PEB Inversion of a DCM under a hidden Markov model of state transitions
% FORMAT [HMM,CSD] = spm_dcm_HMM(DCM,N,b)
% FORMAT [HMM]     = spm_dcm_HMM(CSD,b)
% -------------------------------------------------------------------------
% DCM{p} - DCMs for p sessons: DCM.b encodes state-dependent connections
% N      - number of windows within which to evaluate states
% b{s}   - Cell array state transition priors (Dirichlet parameters) 
%
% returns HMM(s):
%     HMM(s).X  - posterior expectation of hidden states
%     HMM(s).qB - posterior expectation of HMM parameters
%     HMM(s).qb - and Dirichlet concentration parameters
%     HMM(s).qP - posterior expectation of PEB parameters
%     HMM(s).qC - posterior covariances of PEB parameters
%     HMM(s).iP - indices of DCM parameters
%     HMM(s).Ep - posterior expectation of DCM parameters
%     HMM(s).Cp - posterior covariances of DCM parameters
%     HMM(s).L  - free energy components
%     HMM(s).F  - total free energy (model evidence)
% s  -  index of HMM structure (prior model of state transitions)
%
% CSD{N,P} - inverted DCM of each window; with window functions in CSD{n}.W
%__________________________________________________________________________
%
%  This routine characterises a single timeseries in terms of latent or
%  hidden  brain states manifest in terms of state dependent connectivity.
%  It first inverts  the complex cross spectral density of the observed
%  timeseries and  then estimates epoch specific fluctuations in state
%  dependent connectivity in a subset of connections (specified in the
%  logical field DCM.b). The ensuing  sequence of posterior densities are
%  then subject to Bayesian model reduction to provide evidence for
%  sequences of state transitions under a hidden Markov model. Effectively,
%  this involves supplying the evidence that the brain is in a particular
%  connectivity state at each epoch – using the reduced free energy -  to a
%  variational message passing scheme based upon a Markov decision process.
%  The higher (discrete state space for hidden Markov model) level that
%  returns the Bayesian model average for iterative optimisation of the
%  state dependent connection (PEB) parameters, and the epoch specific
%  connectivity (DCM) parameters. The products of this inversion are
%  posteriors at the DCM (epoch specific), PEB, (state specific)and HMM
%  (transition) level. These  posterior densities fully characterise a
%  given time series in terms of discrete  state transitions, where each
%  brain state is associated with a location in (connectivity) parameter
%  space; in other words, a discrete characterisation of dynamic or
%  fluctuating effective connectivity.
%__________________________________________________________________________
% Copyright (C) 2015-2016 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_HMM.m 7279 2018-03-10 21:22:44Z karl $


%  get windowed cross spectra if necessary
%==========================================================================
if nargin == 3
    csd = spm_dcm_window(GCM,N);
else
    csd = GCM;                     % cross spectra of windows
    b   = N;                       % state transition priors
end


%% inversion of hierarchical (empirical) Bayesian HMM model
%==========================================================================
[N,P] = size(csd);                 % number of windows and sessions
a     = spm_find_pC(csd{1},{'A'}); % indices of state-dependent parameters
for k = 1:numel(b)
    
    % initialise this model
    % ---------------------------------------------------------------------
    clear F L M O MDP
    b0      = b{k};                % prior contraints on transitions
    s       = size(b0,1);          % number of hidden states
    D       = sparse(1,1,1,s,1);   % intial state
    
    % initialise succession of hidden states (assume an orbit)
    % ---------------------------------------------------------------------
    X       = kron(ones(1,N),eye(s,s));
    X       = spm_softmax(X(1:s,1:N));
    X       = kron(ones(1,P),X);
    
    % initialise priors on hidden Markov model
    % ---------------------------------------------------------------------
    MDP.A   = {eye(s,s)};          % likelihood of outcomes given hidden states
    MDP.D   = {full(D)};           % initial states
    MDP.b   = {full(b0)};          % transitions among states
    
    % prepare
    % ---------------------------------------------------------------------
    for i = 1:16
        
        % update state-dependent parameters; using PEB
        %==================================================================
        M.X       = X';                          % expected hidden states
        [PEB,CSD] = spm_dcm_peb(csd(:),M,{'A'}); % expected connectivity
        
        % evaluate the likelihood of each state (at each epoch)
        %==================================================================
        for t = 1:numel(CSD)
            qE    = spm_vec(CSD{t}.Ep);
            qP    = spm_inv(CSD{t}.Cp);
            for o = 1:s
                rE     = PEB.Ep(:,o) - qE(a);
                F(o,t) = -rE'*qP(a,a)*rE/2;
            end
        end

        % update expected hidden states, using MDP (HMM)
        %==================================================================
        MDP.O{1} = spm_softmax(F);
        mdp      = spm_MDP_VB_X(MDP);
        
        % update transition probabilities and inital states
        %==================================================================
        MDP.B    = mdp.b;
        X        = mdp.X{1};
        
        % record free energy terms
        % -----------------------------------------------------------------
        L(1,i)   = mdp.F(end);      % HMM: complexity from hidden states
        L(2,i)   = mdp.Fb;          % HMM: complexity from HMM parameters
        L(3,i)   = PEB.F;           % PEB: complexity from PEB parameters
        
        
        % test for convergence
        % -----------------------------------------------------------------
        if i > 4 && norm(L(:,i) - L(:,i - 1)) < 1/128
            break
        end
    end
    
    % store posteriors for this hidden Markov model
    % ---------------------------------------------------------------------
    for t = 1:numel(CSD)
        Ep{t} = CSD{t}.Ep;      % epoch specific (PEB) expectation
        Cp{t} = CSD{t}.Cp;      % epoch specific (PEB) covariances
    end
    
    HMM(k).X  = X;              % posterior expectation of hidden states
    HMM(k).qB = mdp.B{1};       % posterior expectation of HMM parameters
    HMM(k).qb = mdp.b{1};       % and Dirichlet concentration parameters
    HMM(k).qP = PEB.Ep;         % posterior expectation of PEB parameters
    HMM(k).qC = PEB.Cp;         % posterior covariances of PEB parameters
    HMM(k).iP = PEB.Pind;       % indices of DCM parameters
    HMM(k).Ep = Ep;             % posterior expectation of DCM parameters
    HMM(k).Cp = Cp;             % posterior covariances of DCM parameters
    HMM(k).L  = L;              % free energy components
    HMM(k).F  = sum(L(:,end));  % total free energy (model evidence)
    
end


if nargout, return, end

% report analysis
%==========================================================================

%  get figure
%--------------------------------------------------------------------------
spm_figure('Getwin','HMM'); clf
spm_dcm_HMM_plot(HMM)

return


function [csd] = spm_dcm_window(GCM,N)
% Auxiliary function: cross spectra of windowed segments
% FORMAT [HMM,CSD] = spm_dcm_HMM(DCM,N)

% -------------------------------------------------------------------------
% check DCM is a cell array and number of sessions: e.g., participants (P)
%--------------------------------------------------------------------------
if isstruct(GCM), GCM = {GCM}; end
P     = numel(GCM);

% nonlinear system identification (DCM for CSD)
%==========================================================================
% First, invert the entire timeseries to estimate DCM parameters that
% do not change over time. These parameters will be fixed using precise
% shrinkage priors to estimate fluctuating (state-dependent) parameters
% using windowed data to evaluate complex cross spectra
%--------------------------------------------------------------------------

%  invert sessions and take parameter average
%--------------------------------------------------------------------------
for p =1:P
    CSD{p} = spm_dcm_fmri_csd(GCM{p});
end
CSD   = spm_dcm_bpa(CSD(:));

% for each session (e.g., subject)
%--------------------------------------------------------------------------
for p = 1:P
    
    % preliminaries: extract state-dependent parameter specification (b)
    %----------------------------------------------------------------------
    DCM   = GCM{p};
    n     = size(DCM.b,1);
    B     = DCM.b;
    DCM.b = zeros(n,n,0);
    
    % (overlapping) Hanning windows (epochs) W
    %----------------------------------------------------------------------
    nS    = size(DCM.Y.y,1);           % number of scans
    nw    = nS/N;                      % window length
    iw    = (1:N)*nw - nw/2;           % indices of window centres
    jw    = -nw:nw;
    hw    = spm_hanning(length(jw));
    W     = zeros(nS,N);
    for i = 1:N
        k         = iw(i) + jw;
        j         = find(k > 0 & k < nS);
        W(k(j),i) = hw(j);
    end
    
    % use posterior expectations as precise priors that relax hyperpriors
    %----------------------------------------------------------------------
    pC       = spm_zeros(CSD.Ep);      % infinitely precise shrinkage priors
    pC.A     = B/16;                   % state dependent parameters
    DCM.M.pE = CSD.Ep;                 % prior expectation of parameters
    DCM.M.pC = diag(spm_vec(pC));      % prior covariances of parameters
    DCM.M.hE = 0;                      % expected log degrees of freedom
    DCM.M.hC = 1/64;                   % intermediate covariance
    
    % invert each window
    %----------------------------------------------------------------------
    DCM.options.order = 6;
    Y.y   = DCM.Y.y;
    for t = 1:N
        w          = diag(W(:,t));           % get windowing matrix
        w          = w(W(:,t) > 0,:);        % and truncate timeseries
        DCM.Y.y    = w*Y.y;                  % windowed timeseries
        csd{t,p}   = spm_dcm_fmri_csd(DCM);  % invert window
        csd{t,p}.W = sum(w)';                % store window
    end
    
end

