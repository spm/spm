function [Y,X] = spm_COVID_US(P,M,U)
% Generate predictions and hidden states of a multi-region COVID model
% FORMAT [Y,X] = spm_COVID_US(P,M,U)
% P   - model parameters
% M   - model structure (requires M.T - length of timeseries)
% U   - number of output variables [default: 2] or indices e.g., [4 5]
%
% Y(:,1) - number of new deaths
% Y(:,2) - number of new cases
% Y(:,3) - CCU bed occupancy
% Y(:,4) - working days
% Y(:,5) - herd immunity
% Y(:,6) - ...
%
% X{i}    - (M.T x 4) marginal densities over four factors for region i
% location   : {'home','out','CCU','morgue'};
% infection  : {'susceptible','infected','infectious','immune'};
% clinical   : {'asymptomatic','symptoms','ARDS','deceased'};
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
% body of spm_COVID_gen.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% prior connectivity
%--------------------------------------------------------------------------
% According to the Bureau of Transportation Statistics
% (http://www.transtats.bts.gov/), a total of 631,939,829 passengers
% boarded domestic flights in the United States in the year 2010.  This
% averages to 1.73 million passengers flying per day. The population of the
% United States is about 327.2 million: A ratio of 190.
%--------------------------------------------------------------------------

% setup and defaults (assume new deaths and cases as outcome variables)
%--------------------------------------------------------------------------
if (nargin < 3) || isempty(U), U = 1:2; end  % two outcomes
try, M.T; catch, M.T = 365;             end  % over one year

% exponentiate parameters
%--------------------------------------------------------------------------
S     = spm_vecfun(P,@exp);                  % scale parameters

% initial distributions
%--------------------------------------------------------------------------
N     = [M.data.pop];                        % population size
for j = 1:numel(M.data)
    
    n    = exp(M.Q(j).Ep.n);                 % number of initial cases
    m    = exp(M.Q(j).Ep.m);                 % proportion of exposed cases
    r    = exp(M.Q(j).Ep.r);                 % proportion of resistant cases
    n    = S.n(j)*n/N(j);                    % proportion of initial cases
    s    = (1 - n - r);                      % proportion susceptible
    h    = m*3/4;                            % proportion at home
    w    = m*1/4;                            % proportion at work
    p{1} = [h w 0 (1 - m) 0]';               % location
    p{2} = [s n 0 0 r]';                     % infection
    p{3} = [1 0 0 0]';                       % clinical
    p{4} = [1 0 0 0]';                       % testing
    
    % normalise initial marginals
    %----------------------------------------------------------------------
    x{j} = spm_cross(p);                     % initial joint distribution
    
end

% generate outcomes for the specified number of days
%==========================================================================
R     = p;                                   % marginal over regions
for i = 1:M.T
    
    % for this region
    %----------------------------------------------------------------------
    for j = 1:numel(M.data)
        
        % coupling between states
        %------------------------------------------------------------------
        xj    = x{j}(2,:,1,1);
        for k = 1:numel(M.data)
            
            % transport between marginals (at work, asymptomatic and
            % untested)
            %--------------------------------------------------------------
            xk  = x{k}(2,:,1,1);
            xj  = xj*(1 - S.erc(j,k)) + N(k)/N(j)*xk*S.erc(j,k);
            
        end
        xj(xj < 0)    = 0;
        x{j}(2,:,1,1) = xj;
        
        % renormalise joint and marginal densities
        %------------------------------------------------------------------
        x{j} = x{j}/sum(x{j}(:));
        p    = spm_marginal(x{j});

        % state dependent transitions
        %==================================================================
        Q     = M.Q(j).Ep;
        
        % time specific parameters
        %------------------------------------------------------------------
        Q.bas = log(exp(Q.bas) + exp(Q.sus)*spm_phi((i - 32*exp(Q.ont))/16));
        
        % region specific parameters
        %------------------------------------------------------------------
        Q.fed = P.fed;
        
        % update probability distribution
        %------------------------------------------------------------------
        B     = spm_COVID_B(x{j},Q,R);
        x{j}  = spm_unvec(B*spm_vec(x{j}),x{j});
        
        
        % probabilistic mappings: outcomes based on marginal densities (p)
        %==================================================================
        for f = 1:numel(p)
            X{f,j}(i,:) = p{f};
        end
        
        % cumulative number of deaths
        %------------------------------------------------------------------
        Y(i,1,j) = N(j)*p{3}(4);
        
        % cumulative number of positive tests
        %------------------------------------------------------------------
        Y(i,2,j) = N(j)*p{4}(3);
        
        % CCU bed occupancy
        %------------------------------------------------------------------
        Y(i,3,j) = N(j)*p{1}(3);
        
        % working days
        %------------------------------------------------------------------
        Y(i,4,j) = N(j)*p{1}(2);
        
        % herd immunity (proportion)
        %------------------------------------------------------------------
        Y(i,5,j) = p{2}(4);
        
    end
    
    % marginals over regions
    %----------------------------------------------------------------------
    s     = 0;
    for j = 1:numel(M.data)
        s = s + x{j};
    end
    R     = spm_marginal(s/sum(s(:)));
    
end


% retain specified output variables
%--------------------------------------------------------------------------
Y      = Y(:,U,:);

return

