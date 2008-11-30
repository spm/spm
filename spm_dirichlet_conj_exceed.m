function [xp] = spm_dirichlet_conj_exceed (alpha,Nsamp)
% Return p(x1>x2 & x1>x3 & ... & x1 > xk) for p(x)=dirichlet(alpha)
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

% Will Penny
% $Id: spm_dirichlet_conj_exceed.m 2507 2008-11-30 14:45:22Z klaas $


if nargin < 2 | isempty(Nsamp)
    Nsamp = 1e6;
end

% Sample from univariate gamma densities then normalise
% (see Dirichlet entry in Wikipedia or Ferguson (1973) Ann. Stat. 1,
% 209-230)
Nk=length(alpha);
for k=1:Nk,
    r(:,k)=gamrnd(alpha(k),1,Nsamp,1);
end
sr=sum(r,2);
for k=1:Nk,
    r(:,k)=r(:,k)./sr;
end

eval_str=['length(find('];
logic_str=[];
for k=2:Nk,
    add_str=['r(:,1)>r(:,',int2str(k),')'];
    if k<Nk
        and_str=[' & '];
    else
        and_str=[''];
    end
    logic_str=[logic_str,add_str,and_str];
end
eval_str=[eval_str,logic_str,'))/Nsamp'];
xp=eval(eval_str);
            
return
