function [stats,Y,X] = spm_mci_mvnpost (post,method,verbose,max_lag)
% Are MCMC samples consistent with Gaussian posterior ?
% FORMAT [stats,Y,X] = spm_mci_mvnpost (post,method,verbose,max_lag)
%
% post      posterior data structure from spm_mci_post
% method    'ESS' or 'thinning'
% verbose   create plots
% max_lag   maximum potential lag for MAR model
% 
% stats     (multivariate) normal test statistics
%           See spm_mci_mvntest.m
% Y         uncorrelated posterior samples
% X         original posterior samples
% 
% Run Gaussianity test on Markov chain samples
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_mci_mvnpost.m 6697 2016-01-27 14:57:28Z spm $

try, pmax=max_lag; catch, pmax=10; end
try, meth=method; catch, meth='ESS'; end
try, ver=verbose; catch, ver=1; end


j=post.ind;
X=post.P(:,j)';
Nj=length(j);
Np=size(post.P,1);
        
if any(std(X)==0)
    disp('Warning from spm_mci_mvnpost: some parameters have zero variance');
    stats=[];
    return
end

switch meth
    case 'ESS',
        %C=cov(post.P(:,post.ind)');
        for p=1:Np,
            ess(p)=spm_mci_ess(post.P(p,j));
            %v(p)=C(p,p)*ess(p)/Nj;
        end
        %post.Cp=diag(v);
        %mess=mean(ess);
        mess=ceil(min(ess));
        stats = spm_mci_mvntest(X,mess);
        
        if ver
            disp('ESS method');
            disp('Effective Sample Sizes:');
            disp(ess);
        end
        
    case 'thinning',
        for p=1:Np,
            [tmp,m(p)]=spm_mci_ess(post.P(p,j));
        end
        lag=max(m)+1;
        Y=X(1:lag:end,:);
        stats = spm_mci_mvntest(Y);
        if ver
            disp('Thinning method');
            disp(sprintf('Using only every %d-th sample',lag));
        end
        
        
    otherwise
        disp('Unknown method in spm_mci_mvnpost.m');
        return
end