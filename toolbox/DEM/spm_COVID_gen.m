function [Y,X] = spm_COVID_gen(P,M,U)
% generate predictions and hidden states of a COVID model
% FORMAT [Y,X] = spm_COVID_gen(P,M,U)
% P   - model parameters
% M   - model structure (requires M.T - length of timeseries)
% U   - number of output variables [default: 2]
%
% Y(:,1) - number of new deaths
% Y(:,2) - number of new cases
% Y(:,3) - recovered cases
% Y(:,4) - CCU bed occupancy
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
% $Id: spm_COVID_gen.m 7809 2020-03-31 11:55:09Z karl $


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
if (nargin < 3) || isempty(U), U = 2; end         % two outcomes
try, M.T; catch, M.T = 176;           end         % over six months

% initial marginals (Dirichlet parameters)
%--------------------------------------------------------------------------
n    = exp(P.n);        % number of initial cases
N    = exp(P.N)*1e6;    % at risk population size
m    = exp(P.m)*N;      % herd immunity (proportion immune)
Pfat = exp(P.fat);      % P(fatality | CCU)
Psur = exp(P.sur);      % P(survival | home)
N    = N - m;           % number of susceptible cases

p{1} = [8 4 0 0]'; % location:   {'home','out','CCU','norgue'};
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
    Y(i,1) = p{3}(4);
    
    % cumulative number of positive tests
    %----------------------------------------------------------------------
    Y(i,2) = p{4}(3);
    
    % recovery rate times the number of people in CCU with ARDS
    %----------------------------------------------------------------------
    ps     = squeeze(sum(x,[2,4]));
    Pccu   = ps(3,3)/(sum(ps(:,3)) + eps);
    Y(i,3) = Pccu*(p{3}(4)*(1 - Pfat)/Pfat) + ...
             (1 - Pccu)*(p{3}(4)*Psur/(1 - Psur));

    % CCU bed occupancy
    %----------------------------------------------------------------------
    Y(i,4) = p{1}(3);
    
    % ...
    %----------------------------------------------------------------------
    
    
end

% evaluate rates (per day) from cumulative counts 
%--------------------------------------------------------------------------
i      = 1:3;
Y(:,i) = gradient(Y(:,i)')';

% retain specified number of output variables
%--------------------------------------------------------------------------
Y        = N*Y(:,1:U);

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
    I{i} = eye(dim(i));
end

% exponentiate parameters
%--------------------------------------------------------------------------
P = spm_unvec(exp(spm_vec(P)),P);

% probabilistic transitions: location
%==========================================================================
% P.out = 0.3;                      % P(going work | home)
% P.Tcu = 10;                       % period in CCU
% P.sde = 1/128;                    % social distancing threshold
% P.u_b = 16/100000;                % bed availability threshold (per capita)
b       = cell(1,dim(1));

% bed availability and social distancing
%--------------------------------------------------------------------------
Psd  = (1 - p{2}(2))^(32*P.sde);    % P(social distancing)
Pcca = spm_sigma(p{1}(3),P.u_b);    % P(bed available)
Pout = Psd*P.out;                   % P(going work | home)

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
% P.Rin = 2;                        % effective number of contacts: home
% P.Rou = 16;                       % effective number of contacts: work
% P.trn = .8;                       % P(transmission)
% P.Tin = 5;                        % infected (pre-contagious) period
% P.Tcn = 5;                        % contagious period
b       = cell(1,dim(2));

% marginal: infection p{2}
%--------------------------------------------------------------------------
R    = (1 - P.trn*p{2}(3));         % P(no transmission per contact)
Pinf = 1 - R^P.Rin;                 % 1 - P(no transmission) | home
Pcon = 1 - R^P.Rou;                 % 1 - P(no transmission) | work
Kinf = exp(-1/P.Tin);
Kcon = exp(-1/P.Tcn);
    
% marginal: infection {2} | home {1}(1)
%--------------------------------------------------------------------------
%    susceptible  infected  infectious  immune
%--------------------------------------------------------------------------
b{1} = [(1 - Pinf) 0          0          0;
        Pinf       Kinf       0          0;
        0          (1 - Kinf) Kcon       0;
        0          0          (1 - Kcon) 1];
    
% marginal: infection {2} | work {1}(2)
%--------------------------------------------------------------------------
b{2} = [(1 - Pcon) 0          0          0;
        Pcon       Kinf       0          0;
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
% P.dev = 1/2;                      % P(developing symptoms | infected)
% P.sev = 1/5;                      % P(severe symptoms | symptomatic)
% P.Tsy = 8;                        % symptomatic period
% P.Trd = 5;                        % acute RDS   period
% P.fat = 0.2;                      % P(fatality | CCU)
% P.rat = 0.8;                      % P(fatality | home)
b       = cell(1,dim(3));

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
for i = 1:size(b,2)
    if rem(i,numel(p{1}) + 1) == 3
        k      = logical(b(:,i) == (1 - Ksev)*Pfat);
        b(k,i) = (1 - Ksev)*P.fat;
        k      = logical(b(:,i) == (1 - Ksev)*(1 - Pfat));
        b(k,i) = (1 - Ksev)*(1 - P.fat);
    end
end

% kroneckor form
%--------------------------------------------------------------------------
B{3}  = spm_kron({b,I{4}});


% probabilistic transitions: testing
%==========================================================================
% P.u_t = 1/128;           % threshold:   testing capacity
% P.sen = 1/32;           % sensitivity: testing capacity
% P.Tts = 1/2;             % delay:       testing capacity
% P.tes = 1/8;             % relative probability of testing if an infected
b       = cell(1,dim(4));

% test availability and prevalence of symptoms
%--------------------------------------------------------------------------
Ptes = P.sen*spm_sigma(p{4}(2),P.u_t);
Pdia = P.tes*Ptes;
Ktri = exp(-1/P.Tts);    % exp(-1/waiting period)

% marginal: testing {4} | susceptible {2}(1)
%--------------------------------------------------------------------------
%    not tested  waiting       +ve -ve
%--------------------------------------------------------------------------
b{1} = [(1 - Pdia) 0            0   0;
        Pdia       Ktri         0   0;
        0          0            1   0;
        0          (1 - Ktri)   0   1];

% marginal: testing {4} | infected {2}(2)
%--------------------------------------------------------------------------
b{2} = [(1 - Ptes) 0            0   0;
        Ptes       Ktri         0   0;
        0          (1 - Ktri)   1   0;
        0          0            0   1];
    
% marginal: testing {4} | infectious {2}(3)
%--------------------------------------------------------------------------
b{3} = [(1 - Ptes) 0            0   0;
        Ptes       Ktri         0   0;
        0          (1 - Ktri)   1   0;
        0          0            0   1];
    
% marginal: testing {4} | immune {2}(4)
%--------------------------------------------------------------------------
b{4} = [(1 - Pdia) 0            0   0;
        Pdia       Ktri         0   0;
        0          0            1   0;
        0          (1 - Ktri)   0   1];

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



