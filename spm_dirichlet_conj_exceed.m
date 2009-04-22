function xp = spm_dirichlet_conj_exceed(alpha,Nsamp)
% Compute exceedance probabilities for a Dirichlet distribution
% FORMAT xp = spm_dirichlet_conj_exceed(alpha,Nsamp)
% 
% Input:
% alpha     - dirichlet parameters
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
% $Id: spm_dirichlet_conj_exceed.m 3076 2009-04-22 12:11:35Z guillaume $

if nargin < 2
    Nsamp = 1e6;
end

Nk = length(alpha);

% Sample from univariate gamma densities then normalise
% (see Dirichlet entry in Wikipedia or Ferguson (1973) Ann. Stat. 1,
% 209-230)
%==========================================================================
r  = zeros(Nsamp,Nk); % if "out of memory" error, reduce Nsamp
for k = 1:Nk
    r(:,k) = spm_gamrnd(alpha(k),1,Nsamp,1);
end
sr = sum(r,2);
for k = 1:Nk
    r(:,k) = r(:,k)./sr;
end


% Exceedance probabilities:
% For any given model k1, compute the probability that it is more likely
% than any other model k2
%==========================================================================
[y, i] = max(r,[],2);
xp     = histc(i, 1:Nk)' / Nsamp;
