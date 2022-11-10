function [Y,X] = spm_COVID_S(P,M,U)
% Generate predictions and hidden states of a COVID model
% FORMAT [Y,X] = spm_COVID_S(P,M,U)
% P   - model parameters
% M   - model structure (requires M.T - length of timeseries)
% U   - number of output variables [default: 2] or indices e.g., [4 5]
%
% Y(:,1) - number of new deaths
% Y(:,2) - number of new cases
% Y(:,3) - CCU bed occupancy
% Y(:,4) - effective reproduction rate (R)
% Y(:,5) - population immunity (%)
% Y(:,6) - total number of tests
% Y(:,7) - contagion risk (%)
% Y(:,8) - prevalence of infection (%)
% Y(:,9) - number of infected at home, untested and asymptomatic
%
% X      - (M.T x 4) marginal densities over four factors
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

% setup and defaults (assume new deaths and cases as outcome variables)
%--------------------------------------------------------------------------
try, M.T; catch, M.T = 180; end                           % over six months

% demgraphy (Wikipedia)
%--------------------------------------------------------------------------
Age = [
    0   3318449 6.26
    5   2972632 5.61
    10  3080929 5.81
    15  3340265 6.30
    20  3595321 6.78
    25  3650881 6.89
    30  3509221 6.62
    35  3549116 6.69
    40  3885934 7.33
    45  3879815 7.32
    50  3400095 6.41
    55  2996992 5.65
    60  3172277 5.98
    65  2508154 4.73
    70  2044129 3.86
    75  1669345 3.15
    80  1258773 2.37
    85  776311  1.46
    90  403817  0.76];

% bin age ranges and normalise
%--------------------------------------------------------------------------
age = Age(:,3);
pop = [sum(age(1:4)) sum(age(5:8)) sum(age(9:12)) sum(age(13:16)) sum(age(17:end))];
pop = pop/sum(pop);
ns  = numel(pop);



% initialise ensemble densities
%==========================================================================
Q     = spm_vecfun(P,@exp);         % scale parameters
R.n   = [0 0 1 0 0];                % initial cases
R.N   = pop;                        % proportion of census population
for i = 1:ns
    
    % initial marginals (Dirichlet parameters)
    %----------------------------------------------------------------------
    n    = Q.n*R.n(i);              % number of initial cases
    N    = Q.N*R.N(i)*1e6;          % population size
    m    = (1 - Q.m)*N;             % number of sequestered cases
    r    = Q.r(i)*N;                % number of resistant cases
    s    = N - n - r;               % number of susceptible cases
    h    = (N - m)*3/4;             % number at home
    w    = (N - m)*1/4;             % number at work
    p{1} = [h w 0 m 0]';            % location
    p{2} = [s n 0 0 r]';            % infection
    p{3} = [1 0 0 0]';              % clinical
    p{4} = [1 0 0 0]';              % testing

    
    % normalise initial marginals
    %----------------------------------------------------------------------
    Nf    = numel(p);
    for f = 1:Nf
        p{f} = p{f}/sum(p{f});
    end
    
    % ensemble density tensor
    %----------------------------------------------------------------------
    x{i}  = spm_cross(p);
    
end

% identity matrices (precomputed for convenience)
%--------------------------------------------------------------------------
dim   = size(x{1});
I     = cell(ndims(x{1}),1);
for i = 1:ndims(x{1})
    I{i} = speye(dim(i));
end

% solve for density dynamics
%==========================================================================
% https://ourworldindata.org/mortality-risk-covid IFR
%--------------------------------------------------------------------------


% group specific parameters
%--------------------------------------------------------------------------
R.sev = linspace(-8,0,ns);          % log severe symptom rate
R.res = linspace(0,-4,ns);          % log proportion with innate immunity
R.out = log([1 1 1 1/2 1/64]);      % log P(work | home)

% generate timeseries
%--------------------------------------------------------------------------
m     = x{1};                       % initialise marginal over group
for t = 1:M.T
    
    for i = 1:ns
        
        % group-dependent parameters
        %==================================================================
        N     = Q.N*R.N(i)*1e6;               % population size
        S     = P;                            % population log parameters
        S.out = S.out + R.out(i);             
        S.res = S.res + R.res(i);
        S.sev = S.sev + R.sev(i);
        
        % time-dependent parameters
        %==================================================================
        
        % testing
        %------------------------------------------------------------------
        S.bas = log(Q.bas + Q.sus*spm_phi((t - 32*Q.ont)/16));
                      
        % marginal prevalence and occupancy
        %------------------------------------------------------------------
        q    = spm_sum(m,[3 4]);
        q    = q(1:3,:);
        q    = q/sum(q(:));
        Prev = sum(q(:,2));                  % prevalence of infection
        Pcco = sum(q(3,:));                  % CCU occupancy
                
        % coupling between groups (contact)
        %==================================================================
        Pinh  = 1;                           % P(no transmission) | home
        Pinw  = 1;                           % P(no transmission) | work
        for j = 1:ns
            q    = spm_sum(x{j},[3 4]);
            ph   = q(1,:)/sum(q(1,:));       % infection prevalence | home
            pw   = q(2,:)/sum(q(2,:));       % infection prevalence | work
            Pinh = Pinh*(1 - Q.trn*ph(3))^(Q.Rin(i,j));
            Pinw = Pinw*(1 - Q.trn*pw(3))^(Q.Rou(i,j));
        end

        % update ensemble density, with probability dependent transitions
        %------------------------------------------------------------------
        B     = spm_COVID_SB(S,I,dim,Prev,Pcco,Pinh,Pinw);
        x{i}  = spm_unvec(B*spm_vec(x{i}),x{i});
        x{i}  = x{i}/sum(x{i}(:));
        
        % probabilistic mappings: outcomes based on marginal densities (p)
        %==================================================================
        p     = spm_marginal(x{i});
        for f = 1:Nf
            X{f,i}(t,:) = p{f};
        end
        
        % number of daily deaths for group s
        %------------------------------------------------------------------
        D(t,i)  = N*p{3}(4);
        
        % number of daily (positive) tests
        %------------------------------------------------------------------
        C(t,i)  = N*p{4}(3);
        
        % effective reproduction rate (R) (based on infection prevalence)
        %----------------------------------------------------------------------
        RR(t,i) = p{2}(2) + p{2}(3);
        
    end
    
    % marginals over groups
    %----------------------------------------------------------------------
    m     = 0;
    for i = 1:ns
        m = m + R.N(i)*x{i};
    end
    m     = m/sum(m(:));
    
end

% effective reproduction ratio: exp(K*Q.Tcn): K = dln(N)/dt
%----------------------------------------------------------------------
for i = 1:ns
    RR(:,i) = exp(Q.Tcn*gradient(log(RR(:,i) + exp(-8))));
end

% output variables
%--------------------------------------------------------------------------
Y     = [D C RR];
Y     = Y(:,U);



return
