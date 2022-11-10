function [Y,X,Z] = spm_COVID_gen(P,M,U)
% Generate predictions and hidden states of a COVID model
% FORMAT [Y,X,Z] = spm_COVID_gen(P,M,U)
% P    - model parameters
% M    - model structure (requires M.T - length of timeseries)
% U    - number of output variables [default: 2] or indices e.g., [4 5]
% Z{t} - joint density over hidden states at the time t
%
% Y(:,1)  - number of new deaths
% Y(:,2)  - number of new cases
% Y(:,3)  - CCU bed occupancy
% Y(:,4)  - effective reproduction rate (R)
% Y(:,5)  - population immunity (%)
% Y(:,6)  - total number of tests
% Y(:,7)  - contagion risk (%)
% Y(:,8)  - prevalence of infection (%)
% Y(:,9)  - number of infected at home, untested and asymptomatic
% Y(:,10) - new cases per day
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
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


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
try, M.T; catch, M.T = 180;             end         % over six months


% exponentiate parameters
%--------------------------------------------------------------------------
Q    = spm_vecfun(P,@exp);

% initial marginals (Dirichlet parameters)
%--------------------------------------------------------------------------
n    = Q.n;                % number of initial cases
N    = Q.N*1e6;            % population size
m    = (1 - Q.m)*N;        % number of sequestered cases
r    = Q.r*N;              % number of resistant cases
s    = N - n - r;          % number of susceptible cases
h    = (N - m)*3/4;        % number at home
w    = (N - m)*1/4;        % number at work
p{1} = [h w 0 m 0]';       % location 
p{2} = [s n 0 0 r]';       % infection 
p{3} = [1 0 0 0]';         % clinical 
p{4} = [1 0 0 0]';         % testing

% normalise initial marginals
%--------------------------------------------------------------------------
Nf    = numel(p);
for f = 1:Nf
    p{f}  = p{f}/sum(p{f});
end

% ensemble density tensor and solve over the specified number of days
%--------------------------------------------------------------------------
u     = sparse(U,1,1,16,1);
q     = p;
x     = spm_cross(p);
for i = 1:M.T
    
    % time-dependent parameters
    %======================================================================
    
    % baseline testing
    %----------------------------------------------------------------------
    if isfield(M,'R')
        try
            P.bas = log(Q.bas * M.R(i + 2));
        catch
            P.bas = log(Q.bas * max(M.R));
        end
    else
        if isfield(Q,'tts')
            P.bas = log(Q.bas + Q.sus*spm_phi((i - 32*Q.ont)/Q.tts));
        else
            P.bas = log(Q.bas + Q.sus*spm_phi((i - 32*Q.ont)/16));
        end
    end

    % start of trace and track
    %----------------------------------------------------------------------
    if isfield(M,'TTT')
        P.ttt = log(Q.ttt) + log(spm_phi((i - M.TTT)/8));
    end
    
    % early lockdown
    %----------------------------------------------------------------------
    if isfield(M,'TT')
        P.sde = log(Q.sde) - spm_phi((M.TT - i)/8)*4;
    end
    
    % update ensemble density, with probability dependent transitions
    %----------------------------------------------------------------------
    if isfield(M,'r')
        B = spm_COVID_B(x,P,M.r{i});
    else
        B = spm_COVID_B(x,P);
    end
    x     = spm_unvec(B*spm_vec(x),x);
    x     = x/sum(x(:));
    
    % probabilistic mappings: outcomes based on marginal densities (p)
    %======================================================================
    p     = spm_marginal(x);
    for j = 1:Nf
        X{j}(i,:) = p{j};
    end
    
    % number of daily deaths
    %----------------------------------------------------------------------
    Y(i,1) = N*p{3}(4);
    
    % number of daily (positive) tests
    %----------------------------------------------------------------------
    Y(i,2) = N*p{4}(3);

    % CCU bed occupancy
    %----------------------------------------------------------------------
    Y(i,3) = N*p{1}(3);
    
    % effective reproduction rate (R) (based on infection prevalence)
    %----------------------------------------------------------------------
    Y(i,4) = p{2}(2) + p{2}(3);
    
    % population immunity (%)
    %----------------------------------------------------------------------
    Y(i,5) = p{2}(4)*100;
    
    % total number of daily tests (positive or negative)
    %----------------------------------------------------------------------
    Y(i,6) = N*(p{4}(3) + p{4}(4));
    
    % probability of contracting virus (in a class of 15)
    %----------------------------------------------------------------------
    Y(i,7) = (1 - (1 - Q.trn*p{2}(3))^15)*100;
    
    % prevalence of (contagious) infection (%)
    %----------------------------------------------------------------------
    Y(i,8) = p{2}(3)*100;
    
    % number of people at home, asymptomatic, untested but infected
    %----------------------------------------------------------------------
    Y(i,9) = N*x(1,2,1,1);
    
    % number of new cases
    %----------------------------------------------------------------------
    if u(10)
        q       = x;
        B       = B - diag(diag(B));
        q       = spm_unvec(B*spm_vec(q),q);
        q       = spm_marginal(q);
        Y(i,10) = N*q{2}(3);
    end
    
    
    % joint density if requested
    %----------------------------------------------------------------------
    if nargout > 2
        Z{i} = x;
    end

end

% effective reproduction ratio: exp(K*Q.Tcn): K = dln(N)/dt
%--------------------------------------------------------------------------
Y(:,4) = exp(Q.Tcn*gradient(log(Y(:,4))));

% retain specified output variables
%--------------------------------------------------------------------------
Y      = Y(:,U);

return
