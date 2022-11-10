function [Y,X,dates] = spm_CLIMATE_gen(P,M,U,NPI)
% Generate predictions and hidden states of a CLIMATE model
% FORMAT [Y,X,T] = spm_CLIMATE_gen(P,M,U,NPI)
% P    - model parameters
% M    - model structure (M.T - dates or data structure with dates)
% U    - indices of output variables to generate
% NPI  - intervention
%     NPI(i).period = {'dd-mm-yyyy','dd-mm-yyyy'}; % dates to evaluate
%     NPI(i).param  = {'xyz',...};                 % parameter name
%     NPI(i).Q      = (value1,...);                % parameter value
%     NPI(i).dates  = {'dd-mm-yyyy','dd-mm-yyyy'}; % dates of interevention
%
% Y: (output variables)
%     'Average temperature';
%     'Extreme drought';
%     'Exceptional drought';
%     'Rainfall';
%     'Crop production';
%     'Irrigation';
%     'Fertiliser use';
%     'Milk production';
%     'Food prices';
%     'Economic activity';
%     'Income in exposed sector';
%     'Childhood malnutrition';
%     'Crop yield';
%
% X: (latent states)
%      Meteorological (fast)
%      Meteorological (fast)
%      Meteorological (slow)
%      Anthropological activity
%      Primary sector activity
%      Yield
%      Crop production
%      Irrigation
%      Crop resources (fertilisation)
%      Food production
%      Food price
%      Malnutrition
%
% This function returns outcomes Y and their latent states or causes X,
% given the parameters of a generative model P. Generative models of this
% (state space) sort have two parts. The first part concerns fluctuations
% in latent states specified in terms of equations of motion (technically,
% ordinary differential equations). The second part concerns the mapping
% from the latent states to the observable outcomes. This means the
% parameters of the generative model can be divided into the parameters of
% the equations of motion (e.g., rate or time constants) and the parameters
% of the likelihood mapping. In this instance, the likelihood mapping is
% from latent states to log transformed outcomes and is a simple linear
% mapping – such that the coefficients can be thought of as constants and
% regression coefficients.
%
% The parameters of the equations of motion are slightly more complicated
% and define the system at hand, in terms of which latent states can
% influence others. Because we want to evaluate the posterior predictive
% density over future states, we have to specify everything in terms of
% parameters (in the absence of any outside or endogenous inputs). This
% means everything has to be specified in terms of (time invariant)
% parameters, including the initial states at a specified time (d0). This
% also means that one models different scenarios (e.g., interventions) in
% terms of changes in parameters over particular time points, that can be
% specified in an optional argument (NPI).
%
% This model contains 12 latent states and (by coincidence) 12 outputs. The
% 12 latent states are coupled to each other dynamically through their
% equations of motion and then generate outcomes as mixtures of one or more
% latent states. The code below has been annotated to describe the coupling
% among states (that specifies the dynamical part of the model) and the
% coupling from states to output (that specifies the observation part of
% the model).
%
% For a detailed description of the role of each parameter please see
% spm_CLIMATE_priors.
%
% This script contains some auxiliary code (at the end) that allows one to
% examine the effects of changing various parameters by cutting and pasting
% the appropriate section. For a more mathematical rendition of the model,
% the equations of motion – and Jacobian – can be displayed in latex
% formatby putting a breakpoint in the current file (so that the sub
% function can be referenced) and then cutting and pasting the appropriate
% section into the command window.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% preliminaries
%==========================================================================

% setup and defaults
%--------------------------------------------------------------------------
if (nargin < 3) || isempty(U), U = 1:3; end         % two outcomes
if (nargin < 4), NPI = [];              end         % interventions

% get end date; either specified explicitly or via a data structure M.T
%--------------------------------------------------------------------------
if isstruct(M.T)
    D   = M.T;
    d2  = max(spm_vec(D.date));
    U   = [D.U];
else
    
    % otherwise, generate a continuous timeseries from M.T(1) to M.T(2)
    %----------------------------------------------------------------------
    try
        d1 = datenum(M.T{1},'dd-mm-yyyy');
        d2 = datenum(M.T{2},'dd-mm-yyyy');
    catch
        d1 = M.T(1);
        d2 = M.T(2);
    end
end

% dates (assuming the data are sampled at the beginning of each month)
%--------------------------------------------------------------------------
years      = (1955:2155)';
months     = (1:12)';
dates      = [kron(years,ones(size(months))), kron(ones(size(years)),months)];
dates(:,3) = 1;
dates      = datenum(dates);
dates      = dates(dates < (d2 + 365)); 

% unpack and exponentiate parameters
%--------------------------------------------------------------------------
Q    = P;                                % model parameters
Q.P  = spm_vecfun(Q.P,@exp);             % exponentiate flow parameters
R    = Q;                                % baseline parameters

% solve the equations of motion or flow
%==========================================================================
% The first step is to evaluate the trajectory of hidden states. Once this
% has been done, we can then generate the predictions at each point in time
% or the required time points at which data were sampled.
%--------------------------------------------------------------------------
d     = 2;                               % number of integration steps
x     = P.x(:);                          % initial states
f     = @(x,Q)spm_CLIMATE_f(x,Q);        % equations of motion (see below)
dt    = gradient(dates)/(365.25/12);     % time intervals (one month)
for t = 1:numel(dates)
        
    % time-dependent parameters
    %======================================================================
    
    % interventions (NPI)
    %----------------------------------------------------------------------
    for j = 1:numel(NPI)
        
        % start and end dates of an intervention
        %------------------------------------------------------------------
        dstart = datenum(NPI(j).dates{1},'dd-mm-yyyy');
        dfinal = datenum(NPI(j).dates{2},'dd-mm-yyyy');
        if (dates(t) > dstart) && (dates(t) <= dfinal)
            Q.P.(NPI(j).param) = NPI(j).Q;
        else
            Q.P.(NPI(j).param) = R.P.(NPI(j).param);
        end
    end
    
    % Update latent states using local linearisation: updating in a series
    % of steps, if necessary, to accommodate changes in the Jacobian over
    % state space.
    %======================================================================
    for i = 1:d
        dfdx = spm_cat(spm_diff(f,x,Q.P,1));
        x    = x + spm_dx(dfdx,f(x,Q.P),dt(t)/d);
    end
    
    % save states at this point in time
    %----------------------------------------------------------------------
    X(t,:) = x(:)';

end


% generate outcomes
%==========================================================================
% Having integrated time-dependent latent states, it now remains to
% generate observations. In this example, the observations are linear
% mixtures of the latent states (where the predictions are of log
% transformed timeseries). For ease of reference, the outcomes and their
% causes (i.e., latent states) are listed below:
%--------------------------------------------------------------------------
% OUTCOMES:
%--------------------------------------------------------------------------
% Average temperature
% Extreme drought
% Exceptional drought
% Rainfall
% Crop production
% Irrigation
% Fertilisation
% Milk production
% Food prices
% Income in exposed sector
% Childhood malnutrition
% Crop yield
%
% and hidden STATES:
%--------------------------------------------------------------------------
% Meteorological (fast)
% Meteorological (fast)
% Meteorological (slow)
% Anthropological activity
% Primary sector activity
% Yield
% Crop production
% Irrigation
% Fertilisation
% Food production
% Food price
% Malnutrition

% Average temperature: a mixture of fastmeteorological states
%--------------------------------------------------------------------------
Y(:,1) = log(296) + P.Y.t(1) + X(:,1)*P.Y.t(2) + X(:,2)*P.Y.t(3) + X(:,3)*P.Y.t(4);

% Extreme drought: proportional to negative yield (i.e.,natural irrigation)
%--------------------------------------------------------------------------
Y(:,2) = P.Y.d(1) - X(:,6) * P.Y.d(2);

% Exceptional drought: proportional to negative yield
%--------------------------------------------------------------------------
Y(:,3) = P.Y.d(3) - X(:,6) * P.Y.d(4);

% Rainfall: a mixture of meteorological states
%--------------------------------------------------------------------------
Y(:,4) = log(1e4) + P.Y.r(1) + X(:,1)*P.Y.r(2) + X(:,2)*P.Y.r(3) + X(:,3)*P.Y.r(4);

% Crop production: proportional to latent crop production
%--------------------------------------------------------------------------
Y(:,5) = P.Y.c(1) + X(:,7) * P.Y.c(2);

% Irrigation: proportional to latent irrigation
%--------------------------------------------------------------------------
Y(:,6) = P.Y.i(1) + X(:,8) * P.Y.i(2);

% Fertiliser use: proportional to latent fertilisation
%--------------------------------------------------------------------------
Y(:,7) = P.Y.u(1) + X(:,9) * P.Y.u(2);

% Milk production: a mixture of crop and food production
%--------------------------------------------------------------------------
Y(:,8) = P.Y.m(1) + X(:,7)* P.Y.m(2) + X(:,10) * P.Y.m(3);

% Food prices: proportional to latent food price
%--------------------------------------------------------------------------
Y(:,9) = P.Y.f(1) + X(:,11) * P.Y.f(2);

% Income in exposed sector: a mixture of crop and food production
%--------------------------------------------------------------------------
Y(:,10) = P.Y.e(1) + X(:,7)* P.Y.e(2) + X(:,10) * P.Y.e(3);

% Childhood malnutrition: proportional to latent malnutrition
%--------------------------------------------------------------------------
Y(:,11) = P.Y.p(1) + X(:,12) * P.Y.p(2);

% Crop yield: proportional to the latent yield
%--------------------------------------------------------------------------
Y(:,12) = P.Y.y(1) + X(:,6) * P.Y.y(2);


% generate time averages over 6 or 12 months
%--------------------------------------------------------------------------
period = [6 6 6 6 12 6 6 12 6 12 12 6];

n     = size(Y,1);
K6    = toeplitz(sparse(1,1,1,1,n),sparse(1,1:6,1,1,n));
K6    = diag(1./sum(K6,2))*K6;
n     = size(Y,1);
K12   = toeplitz(sparse(1,1,1,1,n),sparse(1,1:12,1,1,n));
K12   = diag(1./sum(K12,2))*K12;
for i = 1:numel(period)
    if period(i) == 6
        Y(:,i) = K6*Y(:,i);
    else
        Y(:,i) = K12*Y(:,i);
    end
end

% retain specified output variables (specified by the vector U)
%==========================================================================
Y = Y(:,U);


% vectorise if data are asynchronous
%==========================================================================
% Up until this point, the data are organised as a matrix of timeseries for
% selected outcome variables. However, if one is generating predictions of
% empirical data, only certain data points (at the appropriate sampling
% times) are required. The inversion scheme assumes that these data are in
% the form of a single vector.
%--------------------------------------------------------------------------
if exist('D','var')
    
    % if there is a data structure fill in the appropriate predictions
    %----------------------------------------------------------------------
    for t = 1:numel(D)
        j      = ismember(dates,D(t).date);
        D(t).Y = Y(j,t);
    end
    
    % and vectorise
    %----------------------------------------------------------------------
    Y  = spm_vec(D.Y);
    
else
    
    % otherwise generate a timeseries for every month between d1 and d2
    %----------------------------------------------------------------------
    j     = dates >= d1 & dates <= d2;
    Y     = Y(j,:);
    X     = X(j,:);
    dates = dates(j);
    
end

return


function dxdt = spm_CLIMATE_f(x,P)
% equations of motion for a climate model
% FORMAT dxdt = spm_CLIMATE_f(x,P)
%
% x   - latent states
% P   - parameter structure, including initial states:
%
% x(1)  - 'Meteorological (fast)'
% x(2)  - 'Meteorological (fast)'
% x(3)  - 'Meteorological (slow)'
% x(4)  - 'Anthropological activity'
% x(5)  - 'Primary sector activity'
% x(6)  - 'Yield'
% x(7)  - 'Crop production'
% x(8)  - 'Irrigation'
% x(9)  - 'Resources (e.g., fertiliser)'
% x(10) - 'Food production'
% x(11) - 'Food price'
% x(12) - 'Malnutrition'
%
% This subroutine contains evaluates the flow of latent states as a
% function of the current state and the parameters of the equations of
% motion. The functional form of these equations defines the formal (i.e.
% form of) priors that lend the model a particular causal architecture. The
% specific causal among the states are specified below in terms of which
% states influence each of the 12 latent states.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% dynamics
%==========================================================================
% dxdt = f(x):
%--------------------------------------------------------------------------
dxdt     = zeros(size(x),'like',x);

% Meteorological (Annual cycle): this is effectively a harmonic oscillator
% whose amplitude is driven by a slow variable that increases with
% anthropomorphic activity – and decays to a fixed point otherwise
%--------------------------------------------------------------------------
dxdt(1)  = x(1)*x(3) + x(2)*2*pi/12;
dxdt(2)  = x(2)*x(3) - x(1)*2*pi/12;
dxdt(3)  = (x(4) - 1)*P.L(1) - (x(3) + P.L(3))*P.L(2);

% Anthropological activity: this can be regarded as the activity generated
% by a population that is a renewal process. Here, the reproduction and
% decay rates can be influenced by food prices and malnutrition,
% respectively
%--------------------------------------------------------------------------
dxdt(4)  = x(4)*(P.R(2) - P.R(4)*x(11)) - x(4)*(P.R(1) + P.R(3)*x(12));
dxdt(5)  = x(5)*P.R(2) - x(5)*(P.R(1));  % spare state

% Yield: this can be regarded as a natural irrigation state that is driven
% by a nonlinear function of a fast meteorological state (which generates
% rainfall). The nonlinear (error) function means that yield falls with
% meteorological variability
%--------------------------------------------------------------------------
dxdt(6)  = erf(P.P(1) - x(2)*P.P(2)) - x(6)*P.P(3);

% Crop production: crop production increases with (anonlinear function of)
% yield, irrigation and fertilisation
%--------------------------------------------------------------------------
dxdt(7)  = (1 + erf(x(6) - P.C(1))*P.C(2) + x(8)*P.C(4))*x(9)*P.C(3) - x(7)*P.C(5);   

% Irrigation: irrigation increases when the yield is in a particular range
%--------------------------------------------------------------------------
dxdt(8)  = exp(-((P.I(1)  - x(6))^2)*P.I(2)) - x(8)*P.I(3);

% Fertilisation: fertiliser (or, more generally, resource) use increases
% with demand; i.e., the difference between anthropological activity and
% latent crop production
%--------------------------------------------------------------------------
dxdt(9)  = (x(4) - x(7))*P.E(1) - x(9)*P.E(2);

% Food production: this increases with demand, provided crop production is
% sufficiently high.
%--------------------------------------------------------------------------
dxdt(10) = (x(4) - 1 - x(10))*x(7)*P.F(1) - x(10)*P.F(2);

% Food price: this reflects the use of resources (i.e., fertilisers)and
% demand for food; i.e., the difference between anthropological activity
% and latent food production
%--------------------------------------------------------------------------
dxdt(11) = x(9)*P.Z(1) - (x(4) - x(10) - 1)*P.Z(2) - x(11)*P.Z(3);

% Malnutrition: increases with food price and decreases with food
% production
%--------------------------------------------------------------------------
dxdt(12) = x(11)*P.M(1) - x(10)*P.M(2) - x(12)*P.M(3);


return

%% Auxiliary code:


%% illustrate the role of priors
%==========================================================================
% Cut-and-paste this section to evaluate the role of different priors in
% shaping the trajectory of latent states. If the data structure (a) is in
% working memory, the selected outcome data will be plotted over the
% predicted trajectories (as litalltle circles). If this structure is not
% available, simply omit the last argument from the call to
% spm_CLIMATE_plot

% create figure and get prior expectations (pE)
%--------------------------------------------------------------------------
spm_figure('GetWin','States'); clf;
pE      = spm_CLIMATE_priors;

% set at the integration times in M, solve the equations of motion
%--------------------------------------------------------------------------
M.T     = {'01-01-1965','01-01-2020'};    % start and finish dates
u       = 11;                             % index of outcome variable
[Y,X,T] = spm_CLIMATE_gen(pE,M,u);

% and plot the results
%--------------------------------------------------------------------------
try
    spm_CLIMATE_plot(Y,X,u,T,A);
catch
    spm_CLIMATE_plot(Y,X,u,T);
end



%% functional form of polynomial flow
%==========================================================================
% The code below illustrates the functional form of the dynamics and their
% accompanying Jacobian. Using symbolic maths, the flow and its derivative
% with respect to states are automatically evaluated. This enables one to
% evaluate the eigenvalues of the Jacobian and identify the possibility of
% exponential divergence of trajectories (i.e., chaotic dynamics).
%--------------------------------------------------------------------------

% get the states and parameters and specify the equations of motion
%--------------------------------------------------------------------------
Q     = spm_CLIMATE_priors;                   % use prior parameters
nx    = spm_length(Q.x);                      % number of states
Q     = Q.P;                                  % retrieve flow parameters
np    = spm_length(Q);                        % number of parameters

% equations of motion accommodating symbolic forms
%--------------------------------------------------------------------------
f     = @(x,P,Q)spm_CLIMATE_f(x,spm_unvec(P,Q));

syms  x [nx 1] 'real'
syms  P [np 1] 'real'
sympref('FloatingPointOutput',1);

disp(spm_unvec(P,Q))

% evaluate flow and Jacobian for display
%==========================================================================

% Equations of motion or flow (fx) and Jacobian (J)
%--------------------------------------------------------------------------
fx    = f(x,P,Q);
for i = 1:nx
    J(:,i) = diff(fx,x(i));
end

% Symbolic and latex form (for cut-and-paste into a maths editor)
%--------------------------------------------------------------------------
disp('fx = ')
disp(fx)
disp(latex(fx)), disp(' ')

disp('J = ')
disp(J)
disp(latex(J)), disp(' ')


%% repeat using numeric values for the parameters
%--------------------------------------------------------------------------
Q     = spm_vecfun(Q,@exp);
fx    = f(x,spm_vec(Q),Q);
for i = 1:nx
    J(:,i) = diff(fx,x(i));
end
v     = eig(J);

% display the functional form of the first eigenvalue
%--------------------------------------------------------------------------
disp('v(1) = ')
disp(v(1))
disp(latex(v(1))), disp(' ')

%--------------------------------------------------------------------------
% If this function can be greater than one, then for some parts of the
% state space, there is exponential divergence of trajectories. This
% implies chaotic dynamics, should that part of state space be in the
% systems contracting set. In this example, the eigenvalue is a function of
% malnutrition, which can be greater or less than one; suggesting that this
% system is chaotic, with divergence of trajectories until malnutrition
% becomes sufficiently high. At this point, the eigenvalue becomes negative
% and the trajectories will converge.


