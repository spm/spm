function xp = spm_dirichlet_exceedance(alpha,Nsamp)
% Compute exceedance probabilities for a Dirichlet distribution
% FORMAT xp = spm_dirichlet_exceedance(alpha,Nsamp)
% 
% Input:
% alpha     - Dirichlet parameters
% Nsamp     - number of samples used to compute xp [default = 1e6]
% 
% Output:
% xp        - exceedance probability
%__________________________________________________________________________
%
% This function computes exceedance probabilities, i.e. for any given model
% k1, the probability that it is more likely than any other model k2.  
% More formally, for k1=1..Nk and for all k2~=k1, it returns p(x_k1>x_k2) 
% given that p(x)=dirichlet(alpha).
% 
% Refs:
% Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ
% Bayesian Model Selection for Group Studies. NeuroImage (in press)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny & Klaas Enno Stephan
% $Id: spm_dirichlet_exceedance.m 3118 2009-05-12 17:37:32Z guillaume $

if nargin < 2
    Nsamp = 1e6;
end

Nk = length(alpha);

% Perform sampling in blocks
%--------------------------------------------------------------------------
blk = ceil(Nsamp*Nk*8 / 2^28);
blk = floor(Nsamp/blk * ones(1,blk));
blk(end) = Nsamp - sum(blk(1:end-1));

xp = zeros(1,Nk);
for i=1:length(blk)
    
    % Sample from univariate gamma densities then normalise
    % (see Dirichlet entry in Wikipedia or Ferguson (1973) Ann. Stat. 1,
    % 209-230)
    %----------------------------------------------------------------------
    r = zeros(blk(i),Nk);
    for k = 1:Nk
        r(:,k) = spm_gamrnd(alpha(k),1,blk(i),1);
    end
    sr = sum(r,2);
    for k = 1:Nk
        r(:,k) = r(:,k)./sr;
    end
    
    % Exceedance probabilities:
    % For any given model k1, compute the probability that it is more
    % likely than any other model k2~=k1
    %----------------------------------------------------------------------
    [y, j] = max(r,[],2);
    xp = xp + histc(j, 1:Nk)';
    
end
xp = xp / Nsamp;
