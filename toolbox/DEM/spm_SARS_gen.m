function [Y,X,Z] = spm_SARS_gen(P,M,U,NPI)
% Generate predictions and hidden states of a COVID model
% FORMAT [Y,X,Z] = spm_COVID_gen(P,M,U)
% P    - model parameters
% M    - model structure (M.T - length of timeseries or data structure)
% U    - number of output variables [default: 2] or indices e.g., [4 5]
% Z{t} - joint density over hidden states at the time t
%
% Y(:,1)  - Daily deaths (28 days)
% Y(:,2)  - Daily confirmed cases
% Y(:,3)  - Mechanical ventilation
% Y(:,4)  - Reproduction ratio (R)
% Y(:,5)  - Seropositive immunityy (%)
% Y(:,6)  - PCR testing rate
% Y(:,7)  - Contagion risk (%)
% Y(:,8)  - Prevalence {%}
% Y(:,9)  - Daily contacts
% Y(:,10) - Daily incidence
% Y(:,11) - Number infected 
% Y(:,12) - Number symptomatic
% Y(:,13) - Mobility (%)
% Y(:,14) - Workplace (%)
% Y(:,15) - Certified deaths
% Y(:,16) - Hospital admissions
% Y(:,17) - Hospital deaths
% Y(:,18) - Non-hospital deaths
% Y(:,19) - Deaths (>60 years)
% Y(:,20) - Deaths (<60 years)
% Y(:,21) - IFR (infection) (%)
% Y(:,22) - IFR (symptoms) (%)
% Y(:,23) - Daily vaccinations
%
% X       - (M.T x 4) marginal densities over four factors
% location   : {'home','out','CCU','morgue','isolation'};
% infection  : {'susceptible','infected','infectious','immune','resistant'};
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
% and cases per day. This routine can be extended to generate other
% outcomes, or indeed consider other factorisations of the probability
% transition matrices. The subroutine (spm_COVID_T) creating the
% probability transition matrices given the current states and model
% parameters defines the generative model. This model structure rests upon
% a mean field approximation to the transition probabilities that,
% crucially, depends upon (usually the marginal) densities in question.
% Working through the code below will show how this model is constructed.
%
% A more detailed description of the generative model can be found in the
% body of the script.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: spm_SARS_gen.m 8029 2020-12-05 13:37:31Z karl $


% The generative model:
%==========================================================================
% In brief, this model generates timeseries data based on a mean field
% approximation to ensemble or population dynamics. The implicit
% probability distributions are over four latent factors, each with several
% levels or states. These factors are sufficient to generate expected
% outcomes; for example, the number of new cases or the number of people
% infected. The first factor is the location of an individual, who can be
% at home, at work, in a critical care unit (CCU), self isolated or in the
% morgue. The second factor is infection status; namely, susceptible to
% infection, infected, infectious or immune. In addition, we include a
% resistant state that does not participate in the transmission of the
% virus. This model assumes that there is a progression from a state of
% susceptibility to immunity, through a period of (pre-contagious)
% infection to an infectious (contagious) status. The third factor is the
% clinical status; namely, asymptomatic, symptomatic, acute respiratory
% distress syndrome (ARDS) or deceased. Again, there is an assumed
% progression from asymptomatic to ARDS, where people with ARDS can either
% recover to an asymptomatic state or not. Finally, the fourth factor
% represents the diagnostic or testing status of. An individual can be
% untested or waiting for the results of a test that can either be positive
% or negative. With this setup, one can be in one of five places, with any
% infectious status, expressing symptoms or not and having test results or
% not. Note that - in this construction - it is possible to be infected and
% yet be asymptomatic. However, the marginal distributions are not
% independent, in virtue of the dynamics that describe the transition among
% states within each factor. Crucially, the transitions within any factor
% depend upon the marginal distribution of other factors. For example, the
% probability of becoming infected, given that one is susceptible to
% infection, depends upon whether one is at home or at work. Similarly, the
% probability of developing symptoms depends upon whether one is infected
% or not. The probability of being tested depends upon whether one is
% symptomatic. Finally, to complete the circular dependency, the
% probability of leaving home to go to work depends upon the number of
% infected people in the population - as a result of social distancing
% (please see main text). These conditional dependencies constitute the
% mean field approximation and enable the dynamics to be solved or
% integrated over time. At any one point in time, the probability of being
% in any combination of the four factors determines what would be observed
% at the population level. For example, the occupancy of the deceased level
% of the clinical factor determines the current number of people who
% contribute to daily deaths. Similarly, the occupancy of the positive
% level of the testing factor determines the daily positive cases reported.

% References
% neutralising antibodies : https://www.nature.com/articles/s41586-020-2012-7
%--------------------------------------------------------------------------
% seropositive: https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30196-1/fulltext
%--------------------------------------------------------------------------
% For 16 patients with serum samples available 14 days or longer after
% symptom onset, rates of seropositivity were 94% for anti-NP IgG (n=15),
% 88% for anti-NP IgM (n=14), 100% for anti-RBD IgG (n=16), and 94% for
% anti-RBD IgM (n=15). Anti-SARS-CoV-2-NP or anti-SARS-CoV-2-RBD IgG levels
% correlated with virus neutralisation titre (R2>0.9). No genome mutations
% were detected on serial samples.
%--------------------------------------------------------------------------
% immunity : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2271881/
%--------------------------------------------------------------------------
% In this group, antibody concentrations started to increase 1 week after
% inoculation and reached a maximum about 1 week later. Thereafter antibody
% titres slowly declined. Although concentrations were still slightly
% raised 1 year later, this did not always prevent reinfection when
% volunteers were then challenged with the homologous virus.
%--------------------------------------------------------------------------
% long-lasting immunity : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2851497/
%--------------------------------------------------------------------------
% Among 176 patients who had had severe acute respiratory syndrome (SARS),
% SARS-specific antibodies were maintained for an average of 2 years, and
% significant reduction of immunoglobulin G-positive percentage and titers
% occurred in the third year. Thus, SARS patients might be susceptible to
% reinfection >3 years after initial exposure.

% setup and defaults (assume new deaths and cases as outcome variables)
%--------------------------------------------------------------------------
if (nargin < 3) || isempty(U), U = 1:2; end         % two outcomes
if (nargin < 4), NPI = [];              end         % interventions
try, M.T; catch, M.T = 180;             end         % over six months

% deal with data structures (asynchronous timeseries)
%--------------------------------------------------------------------------
if isstruct(M.T)
    D   = M.T;
    d   = spm_vec(D.date);
    if isfield(M,'date')
        d0 = datenum(M.date,'dd-mm-yyyy');
    else
        d0 = min(d);
    end
    M.T = max(d) - d0 + 1;
    d   = d0:max(d);
end

% exponentiate parameters
%--------------------------------------------------------------------------
Q    = spm_vecfun(P,@exp);

% initial marginals (Dirichlet parameters)
%--------------------------------------------------------------------------
n    = Q.n;                % number of initial cases
N    = Q.N*1e6;            % population size
m    = N - N*Q.o;          % number of unexposed cases
r    = Q.r*N;              % number of resistant cases
s    = N - n - r;          % number of susceptible cases
h    = (N - m)*3/4;        % number at home
w    = (N - m)*1/4;        % number at work
p{1} = [h w 0 m 0 0]';     % location 
p{2} = [s n 0 0 r]';       % infection 
p{3} = [1 0 0 0]';         % clinical 
p{4} = [1 0 0 0]';         % testing

% normalise initial marginals
%--------------------------------------------------------------------------
Nf    = numel(p);
for f = 1:Nf
    p{f}  = p{f}/sum(p{f});
end

% initial ensemble density
%--------------------------------------------------------------------------
x     = spm_cross(p);
R     = P;
R.n   = -16;
R.t   = 0;
for i = 1:8
    T     = spm_COVID_T(x,R);
    x     = spm_unvec(T*spm_vec(x),x);
    x     = x/sum(x(:));
end

% ensemble density tensor and solve over the specified number of days
%--------------------------------------------------------------------------
Y     = zeros(M.T,32);
for i = 1:M.T
    
    % time-dependent parameters
    %======================================================================
    
    % nonpharmacological interventions
    %----------------------------------------------------------------------
    for j = 1:numel(NPI)
        
        dstart = datenum(NPI(j).dates{1},'dd-mm-yyyy') - datenum(NPI(j).period{1},'dd-mm-yyyy');
        dfinal = datenum(NPI(j).dates{2},'dd-mm-yyyy') - datenum(NPI(j).period{1},'dd-mm-yyyy');
        if (i > dstart) && (i <= dfinal)
            if ischar(NPI(j).param)
                P.(NPI(j).param) = log(NPI(j).Q);
            else
                for k = 1:numel(NPI(j).param)
                    P.(NPI(j).param{k}) = log(NPI(j).Q(k));
                end
            end
        else
            if ischar(NPI(j).param)
                P.(NPI(j).param) = log(Q.(NPI(j).param));
            else
                for k = 1:numel(NPI(j).param)
                    P.(NPI(j).param{k}) = log(Q.(NPI(j).param{k}));
                end
            end
        end
    end
        
    % start of trace and track
    %----------------------------------------------------------------------
    if isfield(M,'TTT')
        if i > M.TTT
            if isfield(M,'FTT')
                P.ttt = log(M.FTT);
            else
                P.ttt = log(Q.ttt) + 2;
            end
        else
            P.ttt = log(Q.ttt);
        end
    end
    
    % circuit breaker
    %----------------------------------------------------------------------
    if isfield(M,'CBT')
        if (i > M.CBT) && (i <= (M.CBT + M.CBD))
            P.sde = log(Q.sde) - log(4);
        else
            P.sde = log(Q.sde);
        end
    end
    
    % update ensemble density (x)
    %----------------------------------------------------------------------
    P.t   = log(i);
    T     = spm_COVID_T(x,P);
    x     = spm_unvec(T*spm_vec(x),x);
    x     = x/sum(x(:));
    
    % incidence (per day) (r)
    %----------------------------------------------------------------------
    R     = T - diag(diag(T));
    r     = spm_unvec(R*spm_vec(x),x);
    u     = spm_marginal(r);
        
    % marginal densities (p)
    %----------------------------------------------------------------------
    p     = spm_marginal(x);
    for j = 1:Nf
        X{j}(i,:) = p{j};
    end
    
    % outcomes
    %======================================================================
    S   = (1 + cos(2*pi*i/365))/2; % seasonal fluctuations
    
    % number of daily deaths
    %----------------------------------------------------------------------
    if isfield(Q,'dc')
        Y(i,1) = N * (S*Q.dc(1) + (1 - S)*Q.dc(2)) * p{3}(4);
    else
        Y(i,1) = N * p{3}(4);
    end

    % number of daily (positive) tests
    %----------------------------------------------------------------------
    Y(i,2) = N * p{4}(3);

    % CCU bed occupancy (mechanical ventilation)
    %----------------------------------------------------------------------
    if isfield(Q,'mv')
        Y(i,3) = N * (S*Q.mv(1) + (1 - S)*Q.mv(2)) * p{1}(3);
    else
        Y(i,3) = N * p{1}(3);
    end

    % effective reproduction ratio (R) (based on infection prevalence)
    %----------------------------------------------------------------------
    Y(i,4) = p{2}(2) + p{2}(3);
    
    % seropositive immunity (%)
    %----------------------------------------------------------------------
    Y(i,5) = p{2}(4)*100;
    
    % total number of daily tests (positive or negative)
    %----------------------------------------------------------------------
    Y(i,6) = N * (p{4}(3) + p{4}(4));
    
    % probability of contracting virus (in a class of 15)
    %----------------------------------------------------------------------
    Y(i,7) = (1 - (1 - Q.trn*p{2}(3))^15)*100;
    
    % prevalence of (contagious) infection (%)
    %----------------------------------------------------------------------
    Y(i,8) = p{2}(3)*100;
    
    % number of people at home, asymptomatic, untested but infected
    %----------------------------------------------------------------------
    Y(i,9) = N * x(1,2,1,1);
    
    % incidence of new cases
    %----------------------------------------------------------------------
    Y(i,10) = 100 * u{2}(3);
    
    % number of infected people
    %----------------------------------------------------------------------
    Y(i,11) = N * (p{2}(2) + p{2}(3));
    
    % number of symptomatic people
    %----------------------------------------------------------------------
    Y(i,12) = N * p{3}(2);

    % mobility (% normal)
    %----------------------------------------------------------------------
    q  = p{1}(2)/(1 - p{1}(4));
    if isfield(Q,'mo')
        Y(i,13) = 100 * Q.mo(1) * q^Q.mo(2);
    else
        Y(i,13) = 100 * q;
    end

    % work (% normal)
    %----------------------------------------------------------------------
    if isfield(Q,'wo')
        Y(i,14) = 100 * Q.wo(1) * q^Q.wo(2);
    else
        Y(i,14) = 100 * q;
    end 

    % certified deaths per day
    %----------------------------------------------------------------------
    Y(i,15) = N * p{3}(4);

    % hospital admissions (symptomatic/ARDS people in hospital/CCU)
    %----------------------------------------------------------------------
    q  = squeeze(spm_sum(r,[2,4]));
    q  = sum(sum(q([3,6],[2,3])));
    if isfield(P,'ho')
        Y(i,16) = N * (S*Q.ho(1) + (1 - S)*Q.ho(2)) * q;
    else
        Y(i,16) = N * q;
    end

    % excess deaths in hospital/CCU
    %----------------------------------------------------------------------
    q       = squeeze(spm_sum(x,[2,4]));
    Y(i,17) = N * sum(q([3,6],4));
    
    % excess deaths not in hospital
    %----------------------------------------------------------------------
    Y(i,18) = N * sum(q([1,2,4,5],4));
    
    % excess deaths > 60 and < 60 (as a function of place of death)
    %----------------------------------------------------------------------
    if isfield(Q,'ag')
        Qag     = spm_softmax(P.ag);
        Y(i,19) = N * Qag(1,:) * q([3,5,6],4);
        Y(i,20) = N * Qag(2,:) * q([3,5,6],4);
    end
    
    % incidence of new infections
    %----------------------------------------------------------------------
    Y(i,21) = u{2}(2);
    
    % incidence of new symptomatic cases
    %----------------------------------------------------------------------
    Y(i,22) = u{3}(2);
    
    % incidence of vaccinations
    %----------------------------------------------------------------------
    q       = squeeze(spm_sum(x,[2,4]));
    Y(i,23) = N * Q.vac*q(6,1);
    
    % joint density if requested
    %----------------------------------------------------------------------
    if nargout > 2
        Z{i} = x;
    end

end

% effective reproduction ratio: exp(K*Q.Tcn): K = dln(N)/dt
%--------------------------------------------------------------------------
Y(:,4)  = exp((Q.Tcn)*gradient(log(Y(:,4))));

% infection fatality ratio (infection and symptomatic cases)
%--------------------------------------------------------------------------
Y(:,21) = (100/N) * cumsum(Y(:,1))./cumsum(Y(:,21));
Y(:,22) = (100/N) * cumsum(Y(:,1))./cumsum(Y(:,22));

% retain specified output variables
%--------------------------------------------------------------------------
Y       = Y(:,U);

% vectorise if data are asynchronous
%--------------------------------------------------------------------------
if exist('D','var')
    for i = 1:numel(D)
        j      = ismember(d,D(i).date);
        D(i).Y = Y(j,i); 
    end
    Y  = spm_vec(D.Y);
end

return
