function [logbf,t] = spm_bms_ttest (y,prior)
% Log Bayes Factor against null for one sample t-test
% FORMAT [logbf,t] = spm_bms_ttest (y,prior)
% 
% y         [N x 1] data vector
% prior     'jzs' (default), 'unit'
%
% logbf     Log Bayes Factor 
% t         t-statistic
%
% Default prior is 'jzs'. The different priors are described in
% [1] Rouder et al. (2009) Psych Bull Rev, 16(2),225-237.
%
% These tests consider the quantity d = mean / standard deviation. 
% The unit information prior, 'unit', uses a zero mean unit variance 
% Gaussian prior over d (the prior variance of d, sig_d^2, is unity).
%
% The (Jeffreys-Zellner-Snow) JZS prior, 'jzs', uses a Cauchy over d 
% and an inverse chi^2 over sig_d^2.
%_______________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_bms_ttest.m 6038 2014-06-04 15:22:42Z will $

try 
    pr=prior;
catch
    pr='jzs';
end

N=length(y);
m=mean(y);
s=std(y);
sem=s/sqrt(N);
t=m/sem;
v=N-1;

switch pr,
    case 'unit',
        rs=1;
    case 'jzs',
        logbf = jzs_logbf(t,v,N);
        return
    otherwise
        disp('Unknown prior in bayes_ttest');
        return
end

% For Unit Information prior, 'unit':
%
% The minimum value of logbf obtains when t=0
% This minimum value is logbf=-0.5*log(1+N*r)
% For r=1 logbf=-3 can be gotten with N=405 samples.

% See note 5 at end of [1]
t2=t^2;
nr=1+N*rs^2;
log_num=-0.5*N*log(1+t2/v);
log_denom1=-0.5*log(nr);
log_denom2=-0.5*N*log(1+t2/(nr*v));

logbf=log_num-log_denom1-log_denom2;
logbf=-logbf;


%-------------------------------------------------------
function [logbf] = jzs_logbf(t,v,N)
% JZS Log Bayes factor 
%
% See equation (1) in 
% Rouder et al. (2009) Psych Bull Rev, 16(2),225-237.

t2=t^2;
num=(1+t2/v)^(-0.5*(v+1));

gmax=1000;
denom=quad(@(g)jzs_integrand(g,t2,N,v),0,gmax);

logbf=log(num/denom);
logbf=-logbf;

%--------------------------------------------------------
function [j] = jzs_integrand(g,t2,N,v)

ng=1+N*g;
t1=ng.^(-0.5);
t2=(1+t2./(ng*v)).^(0.5*(v+1));
t3=(2*pi)^(-0.5)*(g.^(-1.5)).*exp(-1./(2*g));
j=t1.*t2.*t3;
