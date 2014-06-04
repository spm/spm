function [logBF,F] = spm_bms_anova (y,group,prior)
% Log Bayes factor against null for one-way ANOVA
% FORMAT [logBF,F] = spm_bms_anova (y,group,prior)
%
% y         [n x 1] data vector
% group     [n x 1] vector with elements 1,2,3 etc. indicating group
%           membership
% prior     'jzs' (default) or 'unit'
%
% logBF     LogBayesFactor in favour of alternative
%           logBF < -3 : Accept null (no effect)
%           logBF > +3 : Accept alternative (an effect)
% F         F-statistic
%
% Bayesian ANOVA from [1]
% [1] Wetzels et al 2012, A default Bayesian Hypothesis test
% for ANOVA designs, American Statistical Association, 66(2), 104-111.
%
% For a single group this function calls spm_bms_ttest.m
%_______________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_bms_anova.m 6038 2014-06-04 15:22:42Z will $

try
    prior=prior;
catch
    prior='jzs';
end

k=length(unique(group));
if k==1
    logBF=spm_bms_ttest(y,prior);
    return
end

% Make design matrix
% X         [n x k] design matrix for one-way ANOVA (k-levels)
n=length(y);
X=zeros(n,k);
for i=1:n,
    X(i,group(i))=1;
end

beta_est=pinv(X)*y;
yhat=X*beta_est;
e=y-yhat;
R2=1-(std(e)^2/std(y)^2); % proportion of variance explained

F=((n-k)/(k-1))*(R2/(1-R2));

switch prior
    case 'jzs',
        logBF = jzs_logbf_anova(R2,n,k);
    case 'unit',
        % Unit information prior
        g=n;
        
        % From equation 1 in [1]
        BF=(1+g)^(0.5*(n-k-1));
        BF=BF*(1+g*(1-R2))^(-0.5*(n-1));
        logBF=log(BF);
    otherwise
        disp('Unknown prior in spm_bms_anova');
end

%-------------------------------------------------------
function [logbf] = jzs_logbf_anova(R2,N,k)
% JZS Log Bayes factor for one-way ANOVA
%
% See equation (3) in [1]

logt1=0.5*log(N/2)-gammaln(0.5);

gmax=1000;
t2=quad(@(g)jzs_integrand_anova(g,R2,N,k),0,gmax);

logbf=logt1+log(t2);

%--------------------------------------------------------
function [j] = jzs_integrand_anova(g,R2,N,k)

t1=(1+g).^(0.5*(N-k-1));
t2=(1+g*(1-R2)).^(-0.5*(N-1));
t3=(g.^(-1.5)).*exp(-N./(2*g));
j=t1.*t2.*t3;