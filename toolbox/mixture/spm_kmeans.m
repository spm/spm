function [priors,means,covs,post] = spm_kmeans (y,k,method,return_cov)
% K-means clustering
% FORMAT [priors,means,covs,post] = spm_kmeans (y,k,method,return_cov)
%
% y             [N x d] data matrix containing N samples of d-dim data
% k             number of clusters
% method        'uniform', 'points' or 'random' (default) seeding
% return_cov    Set to 1 to return class covariances. Zero otherwise.
%               (default is 1).
%
% priors        [1 x k] vector of class prior probabilities
% means         [k x d] matrix of class means
% covs          [d x d x k] matrix containing class covariances. This
%               matrix is empty if return_covs=0
% post          [N x k] matrix of class labels
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_kmeans.m 1143 2008-02-07 19:33:33Z spm $

if nargin < 3 | isempty(method)
   method='random';
end

if nargin < 4 | isempty(return_cov)
    return_cov=1;
end

[N,d]=size(y);

switch lower(method),
    case 'uniform',
        % Spread seeds out evenly in d-space
        for i=1:d,
            delta(i)=(max(y(:,i))-min(y(:,i)))/(k+1);
        end
        for m=1:k,
            means(m,:)=m*delta;
        end,
    case 'points',
        % Pick k data points at random. Use as class means
        rp=randperm(N);
        for m=1:k,
            means(m,:)=y(rp(m),:);
        end
    case 'fixed-points',
        % Pick first k data points. Use as class means
        rp=[1:N];
        for m=1:k,
            means(m,:)=y(rp(m),:);
        end
    case 'random',
        mu=mean(y);
        C=cov(y);
        means=spm_samp_gauss(mu,C,k);
    otherwise
        disp('Unknown initialisation method');   
end

data=y;
[nmeans, dim] = size(means);
% Matrix to make unit vectors easy to construct
id = eye(nmeans);

[ndata, dimx] = size(data);
niters=100;
tol=0.0001;
% Main loop of algorithm
for n = 1:niters

  % Save old means to check for termination
  old_means = means;
  
  % Calculate posteriors based on existing means
  d2 = (ones(nmeans, 1) * sum((data.^2)', 1))' + ...
        ones(ndata, 1) * sum((means.^2)',1) - ...
        2.*(data*(means'));
  % Assign each point to nearest centre
  [minvals, index] = min(d2', [], 1);
  post = id(index,:);

  num_points = sum(post, 1);
  % Adjust the means based on new posteriors
  for j = 1:nmeans
    if (num_points(j) > 0)
      means(j,:) = sum(data(find(post(:,j)),:), 1)/num_points(j);
    end
  end

  % Error value is total squared distance from cluster means
  e = sum(minvals);
  if n > 1
    % Test for termination
    if max(max(abs(means - old_means))) < tol & ...
        abs(old_e - e) < tol
      break;
    end
  end
  old_e = e;
end

for m=1:k,
   x=y(find(post(:,m)),:);
   Nm=size(x,1);
   if return_cov
       if Nm > 1
           covs(:,:,m)=cov(x-ones(Nm,1)*means(m,:));
       else
           covs(:,:,m)=0.0001*diag(ones(1,d));
       end
   else
       covs=[];
   end
   priors(m)=sum(post(:,m))/N;
end
