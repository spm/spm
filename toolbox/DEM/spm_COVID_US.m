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
% $Id: spm_COVID_US.m 7838 2020-04-23 17:40:45Z karl $

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
    
    if j == 1
        n = S.n/N(j);                        % proportion of initial cases
    else
        n = 0;
    end
    r    = exp(M.Q(j).Ep.r);                 % proportion of resistant cases
    s    = (1 - n - r);                      % proportion susceptible
    p{1} = [3 1 0 0]'/4;                     % location
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
m     = 1:4;                                 % infection status of exchange
for i = 1:M.T
    
    % for this region
    %----------------------------------------------------------------------
    for j = 1:numel(M.data)
        
        % coupling between states
        %------------------------------------------------------------------
        for k = j:numel(M.data)
            
            % transport between marginals
            %--------------------------------------------------------------
            dNk           = x{k}(2,m,1,1)*N(k)*S.erc(j,k);
            dNj           = x{j}(2,m,1,1);
            dNj           = dNj * sum(dNk(:))/sum(dNj(:));
            dN            = (dNj - dNk);
            
            % from j to k
            %--------------------------------------------------------------
            x{k}(2,m,1,1) = x{k}(2,m,1,1) + dN/N(k);
            x{j}(2,m,1,1) = x{j}(2,m,1,1) - dN/N(j);
            
        end
        
        % joint and marginal densities
        %------------------------------------------------------------------
        x{j}(x{j} < 0) = 0;
        x{j} = x{j}/sum(x{j}(:));
        p    = spm_marginal(x{j});

        % state dependent transitions
        %==================================================================
        
        % region specific parameters
        %------------------------------------------------------------------
        Q     = M.Q(j).Ep;
        Q.fed = P.fed;
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

% evaluate rates (per day) from cumulative counts
%--------------------------------------------------------------------------
for i = 1:2
    for j = 1:size(Y,3)
        Y(:,i,j) = gradient(Y(:,i,j));
    end
end

% retain specified output variables
%--------------------------------------------------------------------------
Y      = Y(:,U,:);

return

