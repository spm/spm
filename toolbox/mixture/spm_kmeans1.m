function [mix1] = spm_kmeans1 (y,k)
% K-means clustering for 1-dimensional data
% FORMAT [mix1] = spm_kmeans1 (y,k)
% 
% y          [1 x N] data vector
% k          Number of components
%
% mix1       Returned model
%
% -------------------------------------------------------
% The fields in mix1 are:
% k                The number of components
% m                Vector of means, m=[m_1,m_2,...,m_k]
% v                Vector of variances, v=[v_1,v_2,..,v_k]
% pi               Vector of mixing proportions, pi=[pi_1,pi_2,..,pi_k]
% nloops           Number of iterations used
% assign           Which class data points are assigned to
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_kmeans1.m 1143 2008-02-07 19:33:33Z spm $

y=y(:)';
N=length(y);

% Spread seeds evenly according to CDF
[x,i]=sort(y);
seeds=[1,2*ones(1,k-1)]*N/(2*k);
seeds=ceil(cumsum(seeds));

last_i=ones(1,N);
m=x(seeds);
for loops=1:100,  
 for j=1:k,
   d(j,:)=(y-m(j)).^2;
 end
 [tmp,i]=min(d);
 if sum(i-last_i)==0
   % If assignment is unchanged
   break;
 else
   % Recompute centres
   for j=1:k,
     m(j)=mean(y(i==j));
   end
   last_i=i;
 end
end  

% Compute variances and mixing proportions
for j=1:k,
   v(j)=mean((y(i==j)-m(j)).^2);
   pi(j)=length(y(i==j))/N;
end


mix1.v=v;
mix1.m=m;
mix1.pi=pi;
mix1.k=k;
mix1.nloops=loops;
mix1.assign=i;