function [Y,X] = spm_COVID_gen(P,M,U)
% generate predictions and hidden states of a COVID model
% FORMAT [Y,X] = spm_COVID_gen(P,M,U)
% P   - model parameters
% M   - model structure (requires M.T - length of timeseries)
% U   - number of output variables [default: 2] or indices e.g., [4 5]
%
% Y(:,1) - number of new deaths
% Y(:,2) - number of new cases
% Y(:,3) - CCU bed occupancy
% Y(:,4) - basic reproduction rate
% Y(:,5) - herd immunity
% Y(:,6) - ...
%
% X      - (M.T x 4) marginal densities over four factors
% location   : {'home','out','CCU','norgue'};
% infection  : {'susceptible','infected','infectious','immune'};
% clinical   : {'asymptomatic','symptoms','ARDS','death'};
% diagnostic : {'untested','waiting','positive','negative'}
%
% This function returns data Y and their latent states or causes X, given
% the parameters of a generative model. This model is a mean field
% approximation based upon population or density dynamics with certain
% conditional dependencies among the marginal densities over four factors.
% See SPM_covid_priors details. In brief, this routine transforms model
% parameters to (exponentiated) scale parameters and then generates a
% sequence of jointed densities over four factors, after assembling a state
% dependent probability transition matrix. The number in the timeseries is
% specified by M.T.
%
% Equipped with a time-dependent ensemble density, outcome measures are
% then generated as expected values. These include the rate of (new) deaths
% and cases per day. This routine can be extended to generateother
% outcomes, or indeed consider other factorisations of the probability
% transition matrices. The subroutine (spm_COVID_B) creating the
% probability transition matrices given the current states and model
% parameters defines the generative model. This model structure rests upon
% a mean field approximation to the transition probabilities that,
% crucially, depends upon (usually the marginal) densities in question.
% Working through the code below will show how this model is constructed.
%
% A more detailed description of the generative model can be found in the
% body of the script.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_COVID_gen.m 7812 2020-04-05 12:14:56Z karl $


% The generative model:
%==========================================================================
% In brief, this model generates timeseries data based on a mean field
% approximation to ensemble or population dynamics. The implicit
% probability distributions are over four latent factors, each with four
% levels or states. These factors are sufficient to generate expected
% outcomes; for example, the number of new cases or the number of people
% infected. The first factor is the location of an individual, who can be
% at home, at work, in a critical care unit (CCU) or in the morgue. The
% second factor is infection status; namely, susceptible to infection,
% infected, infectious or immune. This model assumes that there is a
% progression from a state of susceptibility to immunity, through a period
% of (pre-contagious) infection to an infectious (contagious) status. The
% third factor is the clinical status; namely, asymptomatic, symptomatic,
% acute respiratory distress syndrome (ARDS) or deceased. Again, there is
% an assumed progression from asymptomatic to ARDS, where people with ARDS
% can either recover to an asymptomatic state or not. Finally, the fourth
% factor represents the diagnostic or testing status of. An individual can
% be untested or waiting for the results of a test that can either be
% positive or negative. With this setup, one can be in one of four places,
% with any infectious status, expressing symptoms or not and having test
% results or not. Note that – in this construction – it is possible to be
% infected and yet be asymptomatic. However, the marginal distributions are
% not independent, in virtue of the dynamics that describe the transition
% among states within each factor. Crucially, the transitions within any
% factor depend upon the marginal distribution of other factors. For
% example, the probability of becoming infected, given that one is
% susceptible to infection, depends upon whether one is at home or at work.
% Similarly, the probability of developing symptoms depends upon whether
% one is infected or not. The probability of being tested depends upon
% whether one is symptomatic. These three examples are highlighted by the
% curvilinear arrows – denoting that transition probabilities are
% conditioned upon the marginal distributions over other factors. Finally,
% to complete the circular dependency, the probability of leaving home to
% go to work depends upon the number of infected people in the population –
% as a result of social distancing (please see main text). These
% conditional dependencies constitute the mean field approximation and
% enable the dynamics to be solved or integrated over time. At any one
% point in time, the probability of being in any combination of the four
% states determines what would be observed at the population level. For
% example, the occupancy of the deceased level of the clinical factor
% determines the current number of people who have recorded deaths.
% Similarly, the occupancy of the positive level of the testing factor
% determines the expected number of positive cases reported. From these
% expectations, the expected number of new cases per day can be generated.



% setup and defaults (assume new deaths and cases as outcome variables)
%--------------------------------------------------------------------------
if (nargin < 3) || isempty(U), U = 1:2; end         % two outcomes
if numel(U) == 1,              U = 1:U; end
try, M.T; catch, M.T = 180;             end         % over six months

% exponentiate parameters
%--------------------------------------------------------------------------
sP   = spm_vecfun(P,@exp);

% initial marginals (Dirichlet parameters)
%--------------------------------------------------------------------------
n    = sP.n;             % number of initial cases
N    = sP.N*1e6;         % at risk population size
m    = sP.m*N;           % herd immunity (proportion immune)
N    = N - m;           % number of susceptible cases

p{1} = [3 1 0 0]'; % location:   {'home','out','CCU','norgue'};
p{2} = [N n 0 m]'; % infection:  {'susceptible','infected','infectious','immune'};
p{3} = [1 0 0 0]'; % clinical:   {'asymptomatic','symptoms','ARDS','death'};
p{4} = [1 0 0 0]'; % testing:    {'untested','waiting','positive','negative'}

% normalise initial marginals
%--------------------------------------------------------------------------
Nf    = numel(p);
for f = 1:Nf
    p{f}  = p{f}/sum(p{f});
end

% ensemble density tensor and solve over the specified number of days
%--------------------------------------------------------------------------
x     = spm_cross(p);
for i = 1:M.T
    
    % update ensemble density, with probability dependent transitions
    %----------------------------------------------------------------------
    B    = spm_COVID_B(x,P);
    x    = spm_unvec(B*spm_vec(x),x);
    x    = x/sum(x(:));
    
    % probabilistic mappings: outcomes based on marginal densities (p)
    %======================================================================
    p     = spm_marginal(x);
    for j = 1:Nf
        X{j}(i,:) = p{j};
    end
    
    % cumulative number of deaths
    %----------------------------------------------------------------------
    Y(i,1) = N*p{3}(4);
    
    % cumulative number of positive tests
    %----------------------------------------------------------------------
    Y(i,2) = p{4}(3);

    % CCU bed occupancy
    %----------------------------------------------------------------------
    Y(i,3) = N*p{1}(3);
    
    % basic reproduction rate (R0) infection producing contacts per day (R)
    %----------------------------------------------------------------------
    ps     = squeeze(sum(x,[3,4]));
    ps     = ps(:,3)/sum(ps(:,3));                  % P(infectious | location)
    R      = (ps(1)*sP.Rin + ps(2)*sP.Rou*sP.Tin);  % E(infectious contacts)
    Y(i,4) = R*sP.trn*p{2}(1);                      % basic reproduction rate (R0)           

    
    % herd immunity
    %----------------------------------------------------------------------
    Y(i,5) = p{2}(4);
    
end

% evaluate rates (per day) from cumulative counts 
%--------------------------------------------------------------------------
i      = 1:2;
Y(:,i) = gradient(Y(:,i)')';

% retain specified output variables
%--------------------------------------------------------------------------
Y      = Y(:,U);

return


% set up a generative model (state dependent probability transition matrix)
%__________________________________________________________________________

function T = spm_COVID_B(x,P)
% state dependent probability transition matrices
% FORMAT T = spm_COVID_B(p,P)
% x      - probability distributions (tensor)
% P      - model parameters
% 
% B      - probability transition matrices
%
% this subroutine creates a transition probability tensors as a function of
% model parameters and the joint density over four factors, each with four
% levels. With one exception, the transition probabilities of anyone factor
% depend only upon another factor. The exception is the factor modelling
% clinical status, where the transition from acute respiratory distress
% (ARDS) to death depends upon infection status (infected or not
% infected) and location (in a critical careunit will not).
%__________________________________________________________________________

% marginal probabilities
%==========================================================================
p     = spm_marginal(x);

% identity matrices
%--------------------------------------------------------------------------
dim   = size(x);
I     = cell(ndims(x),1);
for i = 1:ndims(x)
    I{i} = speye(dim(i));
end

% exponentiate parameters
%--------------------------------------------------------------------------
P     = spm_vecfun(P,@exp);

% probabilistic transitions: location
%==========================================================================
% P.out                             % P(work | home)
% P.sde                             % social distancing threshold
% P.cap                             % bed availability threshold (per capita)
% bed availability and social distancing
%--------------------------------------------------------------------------
b    = cell(1,dim(1));
Psd  = (1 - p{2}(2))^(32*P.sde);    % P(social distancing)
Pcca = spm_sigma(p{1}(3),P.cap);    % P(bed capacity)
Pout = Psd*P.out;                   % P(work | home)

% marginal: location {1} | asymptomatic {3}(1)
%--------------------------------------------------------------------------
%      home       work       CCU       morgue
%--------------------------------------------------------------------------
b{1} = [(1 - Pout) 1          1          0;
        Pout       0          0          0;
        0          0          0          0;
        0          0          0          1];

% marginal: location {1}  | symptoms {3}(2)
%--------------------------------------------------------------------------
b{2} = [1          1          1          0;
        0          0          0          0;
        0          0          0          0;
        0          0          0          1];

% marginal: location {1}  | ARDS {3}(3)
%--------------------------------------------------------------------------
b{3} = [(1 - Pcca) (1 - Pcca) 0          0;
        0          0          0          0;
        Pcca       Pcca       1          0;
        0          0          0          1];

% marginal: location {1}  | deceased {3}(4)
%--------------------------------------------------------------------------
b{4} = [0          0          0          0;
        0          0          0          0;
        0          0          0          0;
        1          1          1          1];

% kroneckor form (taking care to get the order of factors right)
%--------------------------------------------------------------------------
b     = spm_cat(spm_diag(b));
b     = spm_kron({b,I{2},I{4}});
B{1}  = spm_permute_kron(b,dim([1,3,2,4]),[1,3,2,4]);
        
% probabilistic transitions: infection
%==========================================================================
% P.Rin                             % effective number of contacts: home
% P.Rou                             % effective number of contacts: work
% P.trn                             % P(transmission | contact)
% P.Tin                             % infected (pre-contagious) period
% P.Tcn                             % infectious (contagious) period
% transmission probabilities
%--------------------------------------------------------------------------
b    = cell(1,dim(2));
R    = (1 - P.trn*p{2}(3));         % P(no transmission per contact)
Pinh = R^P.Rin;                     % P(no transmission) | home
Pinw = R^P.Rou;                     % P(no transmission) | work
Kinf = exp(-1/P.Tin);
Kcon = exp(-1/P.Tcn);
    
% marginal: infection {2} | home {1}(1)
%--------------------------------------------------------------------------
%    susceptible  infected  infectious  immune
%--------------------------------------------------------------------------
b{1} = [Pinh       0          0          0;
        (1 - Pinh) Kinf       0          0;
        0          (1 - Kinf) Kcon       0;
        0          0          (1 - Kcon) 1];
    
% marginal: infection {2} | work {1}(2)
%--------------------------------------------------------------------------
b{2} = [Pinw       0          0          0;
        (1 - Pinw) Kinf       0          0;
        0          (1 - Kinf) Kcon       0;
        0          0          (1 - Kcon) 1];

% marginal: infection {2} | CCU {1}(3)
%--------------------------------------------------------------------------
b{3} = [1          0          0          0;
        0          Kinf       0          0;
        0          (1 - Kinf) Kcon       0;
        0          0          (1 - Kcon) 1];

% marginal: infection {2} | morgue {1}(4)
%--------------------------------------------------------------------------
b{4} = [0          0          0          0;
        0          0          0          0;
        0          0          0          0;
        1          1          1          1];

% kroneckor form
%--------------------------------------------------------------------------
b     = spm_cat(spm_diag(b));
b     = spm_kron({b,I{3},I{4}});
B{2}  = spm_permute_kron(b,dim([2,1,3,4]),[2,1,3,4]);


% probabilistic transitions: clinical
%==========================================================================
% P.dev                             % P(developing symptoms | infected)
% P.sev                             % P(severe symptoms | symptomatic)
% P.Tsy                             % symptomatic period
% P.Trd                             % acute RDS   period
% P.fat                             % P(fatality | CCU)
% P.sur                             % P(survival | home)
% probabilities of developing symptoms
%--------------------------------------------------------------------------
b    = cell(1,dim(3));
Psev = P.sev;                       % P(developing symptoms | infected)
Ksym = exp(-1/P.Tsy);               % acute symptomatic rate
Ksev = exp(-1/P.Trd);               % acute RDS rate
Pdev = P.dev;                       % P(symptoms  | infected)
Pfat = 1 - P.sur;                   % baseline fatality rate

% marginal: clinical {3} | susceptible {2}(1)
%--------------------------------------------------------------------------
%  asymptomatic symptomatic acute RDS deceased
%--------------------------------------------------------------------------
b{1} = [1          1          1         0;
        0          0          0         0;
        0          0          0         0;
        0          0          0         1];

% marginal: clinical {3} | infected {2}(2)
%--------------------------------------------------------------------------
b{2} = [(1 - Pdev) (1 - Ksym)*(1 - Psev) (1 - Ksev)*(1 - Pfat) 0;
        Pdev       Ksym                   0                    0;
        0          (1 - Ksym)*Psev        Ksev                 0;
        0          0                     (1 - Ksev)*Pfat       1];
    
% marginal: clinical {3} | infectious {2}(3)
%--------------------------------------------------------------------------
b{3} = [(1 - Pdev) (1 - Ksym)*(1 - Psev) (1 - Ksev)*(1 - Pfat) 0;
        Pdev       Ksym                   0                    0;
        0          (1 - Ksym)*Psev        Ksev                 0;
        0          0                     (1 - Ksev)*Pfat       1];
    
% marginal: clinical {3} | immune {2}(4)
%--------------------------------------------------------------------------
b{4} = [1          1          1         0;
        0          0          0         0;
        0          0          0         0;
        0          0          0         1];

% kroneckor form
%--------------------------------------------------------------------------
b     = spm_cat(spm_diag(b));
b     = spm_kron({b,I{1}});
b     = spm_permute_kron(b,dim([3,2,1]),[3,2,1]);

% location dependent fatalities (P.fat in CCU)
%--------------------------------------------------------------------------
for i = 3:numel(p{1}):length(b)
    k      = logical(b(:,i) == (1 - Ksev)*Pfat);
    b(k,i) = (1 - Ksev)*P.fat;
    k      = logical(b(:,i) == (1 - Ksev)*(1 - Pfat));
    b(k,i) = (1 - Ksev)*(1 - P.fat);
end


% kroneckor form
%--------------------------------------------------------------------------
B{3}  = spm_kron({b,I{4}});


% probabilistic transitions: testing
%==========================================================================
% P.tft                       % threshold:   testing capacity
% P.sen;                      % sensitivity: testing capacity
% P.del                       % delay:       testing capacity
% P.tes                       % relative probability of testing if an infected
% test availability and prevalence of symptoms
%--------------------------------------------------------------------------
b    = cell(1,dim(4));
Ptes = P.sen*spm_sigma(p{4}(2),P.tft);
Pdia = P.tes*Ptes;
Kdel = exp(-1/P.del);         % exp(-1/waiting period)

% marginal: testing {4} | susceptible {2}(1)
%--------------------------------------------------------------------------
%    not tested  waiting       +ve -ve
%--------------------------------------------------------------------------
b{1} = [(1 - Pdia) 0            0   0;
        Pdia       Kdel         0   0;
        0          0            1   0;
        0          (1 - Kdel)   0   1];

% marginal: testing {4} | infected {2}(2)
%--------------------------------------------------------------------------
b{2} = [(1 - Ptes) 0            0   0;
        Ptes       Kdel         0   0;
        0          (1 - Kdel)   1   0;
        0          0            0   1];
    
% marginal: testing {4} | infectious {2}(3)
%--------------------------------------------------------------------------
b{3} = [(1 - Ptes) 0            0   0;
        Ptes       Kdel         0   0;
        0          (1 - Kdel)   1   0;
        0          0            0   1];
    
% marginal: testing {4} | immune {2}(4)
%--------------------------------------------------------------------------
b{4} = [(1 - Pdia) 0            0   0;
        Pdia       Kdel         0   0;
        0          0            1   0;
        0          (1 - Kdel)   0   1];

% kroneckor form
%--------------------------------------------------------------------------
b     = spm_cat(spm_diag(b));
b     = spm_kron({b,I{1},I{3}});
B{4}  = spm_permute_kron(b,dim([4,2,1,3]),[3,2,4,1]);
    
% probability transition matrix
%==========================================================================
T     = 1;
for i = 1:numel(B)
    T =  T*B{i};
end


return

% Auxiliary functions
%__________________________________________________________________________

function p = spm_sigma(x,u,s)
% reverse sigmoid function
% FORMAT p = spm_sigma(p,u)
% x    - probability
% u    - threshold
% u    - sensitivity (default four)
%
% p    - probability (0 < p < 1)
%
% this function is standard sigmoid function but scales the input argument
% by the bias and flips the (biased) input. This provides a monotonically
% decreasing sigmoid function of the input that hits 50 at the threshold
% (u). The scaling ensures the probability at x = 0 is about one, for a
% suitably large sensitivity parameter s.
%--------------------------------------------------------------------------

% default sensitivity
%--------------------------------------------------------------------------
if nargin < 3, s = 4; end

% sigmoid function
%--------------------------------------------------------------------------
p = spm_phi(s*(u - x)/u);

return



