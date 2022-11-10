function [DCM] = DEM_CLIMATE_India
% FORMAT DCM = DEM_CLIMATE_India
%
% Demonstration of climate and sustainability modelling
%__________________________________________________________________________
%
% This demonstration routine illustrates the use of dynamic causal
% modelling to test hypotheses about the causal architecture that couples
% economic activity to meteorological (and climate) variables. This example
% uses various timeseries from different regions (states) in India
% quantified in terms of temperature, rainfall, drought, irrigation,
% economic activity, measures of malnutrition et cetera.
%
% The dynamic causal model in this instance is a nonlinear deterministic
% (mean field) approximation to the expected (i.e., average) evolution of
% various latent states that generate observable data. Crucially, these
% data can be sparse and discontinuous in time. This means that the unknown
% variables of the implicit forward or generative model are the model
% parameters; namely, the rate or time constants coupling latent
% (unobservable) states and the parameters of a likelihood mapping from
% latent states to (observable) outcomes. In other words, because the model
% is deterministic, the latent states are fully determined by the
% parameters of the generative model, where these parameters can include
% the initial states.
%
% The structure and functional form of the model is described in the
% annotated routines that generates timeseries from latent states, where
% the latent states are the solution to ordinary differential equations
% that describe the influence of one latent state on another. For example,
% there is a latent state called 'anthropomorphic activity' (c.f.,
% population size) that drives a slow meteorological variable, which
% increases the amplitude of annual fluctuations in two (fast)
% meteorological variables. The meteorological variables determine the
% natural yield of agricultural activity, which in turn influences the use
% of irrigation and fertilisers. This influences crop production that
% contributes to food production and, ultimately, anthropomorphic activity,
% via a latent state called 'malnutrition'. In short, we have a nonlinear
% dynamical system in which anthropomorphic activity is coupled directly to
% meteorological (i.e., climate-like) states that are vicariously coupled
% back to anthropomorphic states. The implicit separation of timescales
% within the meteorological states results in itinerant dynamics at a fast
% (seasonal) and slow (decades) timescale.
%
% This routine first illustrates the way in which data are assembled and
% sorted. Here, we average away random fluctuations by averaging over
% regions and then log transform any non-negative data. This allows one to
% use a simple likelihood model with additive Gaussian noise. We next
% assemble the prior density over the parameters and hyperparameters
% controlling the amplitude of observation noise. finally, we invert the
% model to estimate the parameters in terms of a posterior density. The
% posterior density over model parameters can then be used in a variety of
% ways:
%
% First, one can ask whether any parameters are redundant. In other words,
% is the model too expressive or over-parameterised. In the example
% provided, we ask a slightly subtler question: did this parameter need to
% be a free parameter or could we have fixed it to its prior expectation.
% This question can be asked by comparing models with uninformative and
% very precise shrinkage priors over each parameter or combinations of
% parameters. The same kind of comparison can also be used to test
% hypotheses by comparing the log evidence of a model with and without a
% particular link or parameter. In Bayesian model reduction, different
% models correspond to different shrinkage priors (i.e., a model with out a
% particular parameter is specified with priors that shrink it towards a
% small value).
%
% The posterior density can also be used to assess the role of different
% parameters, in generating outcomes, using a straightforward sensitivity
% analysis. This is based on the change in an outcome produced by a change
% in the parameter, where the outcome is a function of time. Finally, one
% can integrate (i.e., solve) the model using the posterior estimates of
% the model parameters to predict what might happen in the future, under
% different scenarios or plausible interventions.
%
% Details about the model structure and assumptions (i.e., priors) can
% be found in the key routines that return the priors (i.e.,
% spm_CLIMATE_priors) and the routine that generates predicted outcomes,
% as a function of the parameters (i.e., spm_CLIMATE_gen). The remaining
% routines are part of the standard SPM software; including the model
% inversion or deconvolution scheme which, in this deterministic setting,
% rests on something called variational Laplace
%
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% Load illustrative data
%==========================================================================
D      = readtable('Indian_data.csv');

% specify (12) outcome variables to explain
%--------------------------------------------------------------------------
vnames = D.Properties.VariableNames;
sj     = find(ismember(vnames,'State_name'));
state  = unique(D(:,sj));
vars   = vnames(4:end);
vars   = vars([1 3 5 6 9 10 11 12 15 19 21 7]);

%% organise via region or state and type of outcome
%--------------------------------------------------------------------------
Y     = struct([]);
for i = 1:numel(state)
    
    % get state
    %----------------------------------------------------------------------
    s   = find(ismember(D(:,sj),state(i,1)));
    disp(state(i,1))
    
    d     = NaN(size(s));
    in    = zeros(size(s));
    is    = zeros(size(s));
    for j = 1:numel(s)
        year     = table2array(D(s(j),find(ismember(vnames,'Year'))));
        month    = table2array(D(s(j),find(ismember(vnames,'season'))));
        year     = num2str(year);
        
        % test whether data are specified in seasons or years
        %------------------------------------------------------------------
        if any(ismember({'NaN','All'},month))
            d(j) = datenum(['06-' year],'mm-yyyy');
            in(j) = 1;
        end
        if any(ismember({'Kharif'},month))
            d(j) = datenum(['03-' year],'mm-yyyy');
            is(j) = 1;
        end
        if any(ismember({'Rabi'},month))
            d(j) = datenum(['09-' year],'mm-yyyy');
            is(j) = 1;
        end
    end
    
    in  = find(in);
    is  = find(is);
    
    % create data structure for this region and outcome variable
    %----------------------------------------------------------------------
    for v = 1:numel(vars)
        
        y     = table2array(D(s,find(ismember(vnames,vars{v}))));
        ns    = sum(isfinite(y(is)));
        
        % retain seasonal data if supplied
        %------------------------------------------------------------------
        if ns > 8
            j = is;
            Y(v,i).period = 6;               % months
        else
            j = in;
            Y(v,i).period = 12;              % months
        end
        Y(v,i).state = table2array(state(i,1));
        Y(v,i).type  = vars{v};
        Y(v,i).unit  = 'Normalised';
        Y(v,i).U     = v;
        Y(v,i).Y     = y(j);
        Y(v,i).date  = d(j);            
        
    end
end

%% Select states to averge
%==========================================================================

% retain regions with at least eight outcome variables
%--------------------------------------------------------------------------
d0    = datestr(min(spm_vec(Y.date)),'dd-mm-yyyy');
for i = 1:size(Y,2)
    ny(i) = numel(spm_COVID_Y(Y(:,i),d0,0));
end
Y(:,ny < 8) = [];

% retain regions with seasonal data
%--------------------------------------------------------------------------
for i = 1:size(Y,2)
    ns(i) = sum([Y(:,i).period] == 6);
end
Y(:,ns < 12) = [];


%% Create grand average (A) over regions for each outcome variable
%==========================================================================
for i = 1:size(Y,1)
    
    dates      = unique(spm_vec({Y(i,:).date}));
    
    A(i).state = 'Average';                    % region or state
    A(i).type  = Y(i).type;                    % type of variable
    A(i).U     = Y(i).U;                       % index in spm_CLIMATE_gen
    A(i).unit  = Y(i).unit;                    % units
    A(i).date  = dates;                        % dates of outcomes
    A(i).Y     = NaN(size(dates));             % outcome variables
    
    % average definitive data
    %----------------------------------------------------------------------
    for d = 1:numel(dates)
        y = [];
        for s = 1:size(Y,2)
            t = find(ismember(Y(i,s).date,dates(d)),1);
            y = [y Y(i,s).Y(t)];
        end
        s = ~(isnan(y) | y == 0);
        if any(s)
            A(i).Y(d,1) = mean(y(s));
        end
    end
end

%% average crop production over seasons (to simplify the modelling)
%==========================================================================
for i = 5
    dv    = datevec(A(i).date);
    year  = unique(dv(:,1));
    for j = 1:numel(year)
        k = find(dv == year(j));
        A(i).Y(k) = mean(A(i).Y(k));
    end
end


%% log transform data (so latent states model proportional changes)
%==========================================================================
A     = spm_COVID_Y(A(:),min(dates),0);
for i = 1:size(A,1)
    if min(A(i).Y > 0)
        A(i).Y = log(A(i).Y);
    end
end

% organise and sort timeseries, removing NaNs
%--------------------------------------------------------------------------
A     = spm_COVID_Y(A(:),min(dates),0);

%% graphics: show the selected outcome variables prior to fitting
%==========================================================================
spm_figure('GetWin','Time series (average)'); clf;

d1 = min(spm_vec({A.date}));             % start states
dT = max(spm_vec({A.date}));             % end date
for i = 1:12
    subplot(4,3,i)
    plot(A(i).date,A(i).Y,'.')
    set(gca,'XLim',[d1,dT])
    datetick('x','yyyy','keeplimits')
    title(A(i).type)
end

% This concludes the assembly of data features that will inform the model.
% In the next section, we specify the priors over the model parameters (and
% precision parameters) and assemble the requisite information into a model
% structure (M). An accompanying data structure (xY) is then passed to a
% standard variational scheme (spm_nsli_GN) that returns the posteriors
% over the unknown parameters. This variational Laplace scheme assumes the
% posteriors and priors are Gaussian. This means that they are sufficiently
% specified in terms of means and covariances.


%% fit averaged data
%==========================================================================

% data structure with vectorised data and expected likelihood precision
%--------------------------------------------------------------------------
xY.y  = spm_vec(A.Y);
xY.Q  = spm_Ce([A.n]);
hE    = zeros(numel(A),1);
for i = 1:numel(A)
    hE(i) = 4 - log(var(A(i).Y));
end

% explain crop production accurately by increasing its log precision
%--------------------------------------------------------------------------
hE(5) = hE(5) + 4;


%% get and set priors over model parameters
%----------------------------------------------------------------------
[pE,pC,str] = spm_CLIMATE_priors;

% model specification
%======================================================================
M.Nmax  = 32;                   % maximum number of iterations
M.G     = @spm_CLIMATE_gen;     % generative function
M.pE    = pE;                   % prior expectations (parameters)
M.pC    = pC;                   % prior covariances  (parameters)
M.hE    = hE;                   % prior expectation  (log-precision)
M.hC    = 1/512;                % prior covariances  (log-precision)
M.T     = A;                    % original data structure and times
U       = [A.U];                % outputs to model

% model inversion with Variational Laplace (Gauss Newton)
%----------------------------------------------------------------------
[Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,xY);

% save prior and posterior estimates (and log evidence (F))
%----------------------------------------------------------------------
DCM.M  = M;
DCM.Ep = Ep;
DCM.Eh = Eh;
DCM.Cp = Cp;
DCM.Y  = A;
DCM.xY = xY;
DCM.F  = F;

%% illustrate posterior predictions of hidden states
%==========================================================================
spm_figure('GetWin','States'); clf;
%--------------------------------------------------------------------------
M.T     = {'01-01-1960', '01-01-2020'};       % time period to plot
u       = 5;                                  % output variable to plot
[Y,X,T] = spm_CLIMATE_gen(Ep,M,u);            % generate predictions
spm_CLIMATE_plot(Y,X,u,T,A);                  % plot with empirical data

% and predicted outcomes (with confidence intervals)
%--------------------------------------------------------------------------
spm_figure('GetWin','outcomes (1)'); clf
global CIPLOT, CIPLOT = false;

j     = 0;
for i = 1:numel(U)
    j = j + 1;
    subplot(4,3,j)
    spm_CLIMATE_ci(Ep,Cp,A,U(i),M);
    title(str.outcome{U(i)}), ylabel(A(i).unit)
end

% This concludes the illustration of model inversion or fitting. We next
% turn to various applications of the ensuing posterior estimates:
% illustrating the use of Bayesian model comparison, sensitivity analyses
% and the use of posterior predictive densities over outcomes in the
% future, for forecasting and scenario modelling.


%% An example of Bayesian model reduction to find redundant parameters
%==========================================================================
% This section illustrates the use of Bayesian model reduction to compare
% the evidence for models with and without certain parameters. This can be
% useful for identifying redundant parameters – or ensuring that all the
% model parameters are necessary to explain the data. In this illustration,
% we evaluate the evidence for all models, where different models are
% defined by switching off different combinations of model parameters. One
% then gathers together the marginal likelihood or evidence for all models
% with and without a particular parameter to establish the posterior
% probability that the parameter was necessary to explain the data.
% Switching parameters on and off, in this case, involves reducing the
% prior covariance to a small value but leaving the prior expectation
% unchanged. In effect, this converts a free parameter into a fixed
% parameter.

% Bayesian model reduction, testing all combinations of model parameters
%--------------------------------------------------------------------------
% FORMAT [DCM,BMR,BMA] = spm_dcm_bmr_all(DCM,field,OPT)
%
% DCM      - A single estimated DCM (or PEB) structure:
%
%  DCM.M.pE  - prior expectation
%  DCM.M.pC  - prior covariance
%  DCM.Ep    - posterior expectation
%  DCM.Cp    - posterior covariances
%  DCM.beta  - prior expectation of reduced parameters (default: 0)
%              NB: beta  = 'full' uses full prior expectations
%  DCM.gamma - prior variance    of reduced parameters (default: 0)
%              NB: gamma = 'full' uses full prior variances
%--------------------------------------------------------------------------
DCM.beta  = 'full';              % retain prior expectations
DCM.gamma = 0;                   % but remove any prior uncertainty
[RCM,BMR] = spm_dcm_bmr_all(DCM,'P');

% This analysis suggests that we could simplify the model by fixing a few
% of the parameters to their prior values. Any parameter that has a zero
% posterior in the lower left bar chart could have been fixed to its prior
% expectation to minimise model complexity, without compromising model
% accuracy. Note that the (log) model evidence or marginal likelihood
% corresponds to accuracy minus complexity.


%% An example of Bayesian model reduction to test hypotheses
%==========================================================================
% This section illustrates the use of Bayesian model reduction to compare
% the evidence for models with and without a key parameter defined in terms
% of the priors over model parameters. In this example, we compare models
% with and without a coupling from anthropomorphic activity to (slow)
% meteorological states.

% full priors and posteriors
%----------------------------------------------------------------------
pE     = DCM.M.pE;
pC     = DCM.M.pC;
qE     = DCM.Ep;
qC     = DCM.Cp;

% reduced priors (that suppress the anthropomorphic effect)
%----------------------------------------------------------------------
rE        = pE;
rC        = pC;
rE.P.L(1) = rE.P.L(1) + log(1/10);

% evaluate the (log) evidence of the reduced model
%----------------------------------------------------------------------
F   = spm_log_evidence(qE,qC,pE,pC,rE,rC);

disp('F (variational evidence lower bound) in natural units')
disp(F)

% In this instance, the log evidence is negative. This means the reduced
% model (without an anthropomorphic drive to climate change) has
% substantially less evidence than a model with an anthropomorphic effect.
% Here, the reduced model was specified by reducing the prior expectation
% of coupling between human activity and climate by an order of magnitude.
% Note that because we are dealing with log scale parameters, this
% reduction corresponds to a subtraction (here, of the log of 1/10).
% Generally speaking, a reduction in log evidence of three or more is taken
% as strong evidence against the reduced hypothesis. This is because the
% corresponding Bayes factor (i.e., odds ratio) is a less than one in 20:
% i.e., one in exp(3). In short, we can conclude there is (very)
% strong evidence for an anthropomorphic effect, under this model.
    
    
%% Sensitivity analysis

%==========================================================================
% This section illustrates the use of dynamic causal modelling to ask which
% parameters are key in shaping outcomes. Here, we ask which are the most
% important parameters of the dynamics (the flow parameters in Ep.P) in
% determining cumulative malnutrition (U = 11) over the next 10 years.

% sensitivity analysis in terms of partial derivatives
%--------------------------------------------------------------------------
spm_figure('GetWin','sensitivity analysis'); clf;

period = {'01-01-2020','01-01-2030'};
M.T    = period;
dYdP   = spm_diff(@(P,M,U)spm_CLIMATE_gen(P,M,U),Ep,M,11,1);

% cumulative effects over time
%--------------------------------------------------------------------------
DY    = sum(dYdP);
i     = spm_fieldindices(Ep,'P');

% plot resultsfor the parameters in Ep.P
%--------------------------------------------------------------------------
subplot(2,1,1)
bar(DY(i(1:7)))
title('Sensitivity analysis','FontSize',16), box off
xlabel('parameters');
ylabel('cumulative child malnutrition','FontSize',16), box off

% The first parameter is the sensitivity of the slow meteorological (i.e.
% climate) state to fluctuations in anthropological activity. If this
% sensitivity is positive, then increasing the influence of human activity
% on climate will increase malnutrition.


%% An example of predictive modelling
%==========================================================================
% This section illustrates the use of dynamic causal modelling to rollout
% into the future under different scenarios. In this example, we ask what
% would happen if the link between anthropomorphic activity and climate was
% broken. More specifically, we can compare posterior predictive densities
% over future outcomes with and without a coupling between anthropological
% activity and the slow aetiological's date by setting the first parameter
% to 0

% time period to consider
%--------------------------------------------------------------------------
period = {'01-01-1990','01-01-2035'};

% scenario: reduction in the reproduction rate on 01-01-2022
%--------------------------------------------------------------------------
Q          = spm_vecfun(Ep.P.L,@exp);
Q(1)       = 0;

NPI.period = period;
NPI.param  = 'L';
NPI.Q      = Q;
NPI.dates  = {'01-01-2020',period{2}};

% plot forecasts of hidden states and crop production
%==========================================================================
spm_figure('GetWin','scenario modelling'); clf;
%--------------------------------------------------------------------------
M.T     = period;
u       = 11;
[Y,X,T] = spm_CLIMATE_gen(Ep,M,u);
spm_CLIMATE_plot(Y,X,u,T,A);

% the effect of intervention (NPI) on latent states
%--------------------------------------------------------------------------
[Y,X,T] = spm_CLIMATE_gen(Ep,M,u,NPI);
spm_CLIMATE_plot(Y,X,u,T,A);

% and subsequent outcomes
%--------------------------------------------------------------------------
spm_figure('GetWin','with confidence intervals'); clf;
%--------------------------------------------------------------------------
CIPLOT = true;

subplot(2,1,1)
spm_CLIMATE_ci(Ep,Cp,A,u,M); hold on
spm_CLIMATE_ci(Ep,Cp,A,u,M,NPI);

return



