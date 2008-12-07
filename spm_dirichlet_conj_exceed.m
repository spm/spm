function [xp] = spm_dirichlet_conj_exceed (alpha,Nsamp)
% Computes exceedance probabilities, i.e. for any given model k1, the 
% probability that it is more likely than any other model k2.  
% More formally, for k1=1..Nk and for all k2~=k1, it returns p(x_k1>x_k2) 
% given that p(x)=dirichlet(alpha).
% 
% FORMAT [xp] = spm_dirichlet_conj_exceed (alpha,Nsamp)
% 
% INPUT:
% alpha   - dirichlet parameters
% Nsamp   - number of samples used to compute xp (default = 1e6)
% 
% OUTPUT:
% xp      - exceedance probability
% 
% REFERENCE:
% Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ
% Bayesian Model Selection for Group Studies. NeuroImage (under review)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny & Klaas Enno Stephan
% $Id: spm_dirichlet_conj_exceed.m 2533 2008-12-07 20:43:53Z klaas $


if nargin < 2 | isempty(Nsamp)
    Nsamp = 1e6;
end

% Sample from univariate gamma densities then normalise
% (see Dirichlet entry in Wikipedia or Ferguson (1973) Ann. Stat. 1,
% 209-230)
%==========================================================================
Nk = length(alpha);
for k = 1:Nk,
    r(:,k) = gamrnd(alpha(k),1,Nsamp,1);
end
sr = sum(r,2);
for k = 1:Nk,
    r(:,k) = r(:,k)./sr;
end


% Exceedance probabilities:
% For any given model k1, compute the probability that it is more likely
% than any other model k2
%==========================================================================
xp = [];
for k1 = 1:Nk,
    eval_str  = ['length(find('];
    logic_str = [];
    ind       = setdiff([1:Nk],[k1]);
    for k2 = ind,
        add_str = ['r(:,' int2str(k1) ')>r(:,' int2str(k2),')'];
        if (k2 < Nk) & ~((k1 == Nk) & (k2 == Nk-1))
            and_str = [' & '];
        else
            and_str = [''];
        end
        logic_str   = [logic_str,add_str,and_str];
    end
    eval_str = [eval_str,logic_str,'))/Nsamp'];
    xp(k1)   = eval(eval_str);
end
xp = xp';


return
