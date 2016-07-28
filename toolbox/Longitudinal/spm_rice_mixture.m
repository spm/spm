function [mg,nu,sig] = spm_rice_mixture(h,x,K)
% Fit a mixture of Ricians to a histogram
% FORMAT [mg,nu,sig] = rice_mixture(h,x,K)
% h   - histogram counts
% x   - bin positions (plot(x,h) to see the histogram)
% K   - number of Ricians
% mg  - integral under each Rician
% nu  - "mean" parameter of each Rician
% sig - "standard deviation" parameter of each Rician
%
% An EM algorithm is used, which involves alternating between computing
% belonging probabilities, and then the parameters of the Ricians.
% The Koay inversion technique is used to compute the Rician parameters
% from the sample means and standard deviations. This is described at
% http://en.wikipedia.org/wiki/Rician_distribution
%_______________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_rice_mixture.m 6844 2016-07-28 20:02:34Z john $

mg  = ones(K,1)/K;
nu  = (0:(K-1))'*max(x)/(K+1);
sig = ones(K,1)*max(x)/K;

m0 = zeros(K,1);
m1 = zeros(K,1);
m2 = zeros(K,1);
ll = -Inf;
for iter=1:10000,
    p  = zeros(numel(x),K);
    for k=1:K,
        % Product Rule
        % p(class=k, x | mg, nu, sig) = p(class=k|mg) p(x | nu, sig, class=k)
        p(:,k) = mg(k)*ricepdf(x(:),nu(k),sig(k)^2);
    end

    % Sum Rule
    % p(x | mg, nu, sig) = \sum_k p(class=k, x | mg, nu, sig)
    sp  = sum(p,2)+eps;
    oll = ll;
    ll  = sum(log(sp).*h(:)); % Log-likelihood
    if ll-oll<1e-8*sum(h), break; end

    %fprintf('%g\n',ll);
    %md = mean(diff(x));
    %plot(x(:),p,'--',x(:),h/sum(h)/md,'b.',x(:),sp,'r'); drawnow

    % Bayes Rule
    % p(class=k | x, mg, nu, sig) = p(class=k, x | mg, nu, sig) / p(x | mg, nu, sig)
    p = bsxfun(@rdivide,p,sp);


    % Compute moments from the histograms, weighted by the responsibilities (p).
    for k=1:K,
        m0(k) = sum(p(:,k).*h(:));              % Number of voxels in class k
        m1(k) = sum(p(:,k).*h(:).*x(:));        % Sum of the intensities in class k
        m2(k) = sum(p(:,k).*h(:).*x(:).*x(:));  % Sum of squares of intensities in class k
    end

    mg = m0/sum(m0); % Mixing proportions
    for k=1:K,
        mu1 = m1(k)./m0(k);                                % Mean 
        mu2 = (m2(k)-m1(k)*m1(k)/m0(k)+1e-6)/(m0(k)+1e-6); % Variance

        % Compute nu & sig from mean and variance
        [nu(k),sig(k)] = moments2param(mu1,mu2);
    end
    %disp([nu'; sig'])
end
%_______________________________________________________________________

%_______________________________________________________________________

function [nu,sig] = moments2param(mu1,mu2)
% Rician parameter estimation (nu & sig) from mean (mu1) and variance
% (mu2) via the Koay inversion technique.
% This follows the scheme at
% https://en.wikipedia.org/wiki/Rice_distribution#Parameter_estimation_.28the_Koay_inversion_technique.29
% This Wikipedia description is based on:
% Koay, C.G. and Basser, P. J., Analytically exact correction scheme
% for signal extraction from noisy magnitude MR signals,
% Journal of Magnetic Resonance, Volume 179, Issue = 2, p. 317–322, (2006)

r     = mu1/sqrt(mu2);
theta = sqrt(pi/(4-pi));
if r>theta,
    for i=1:256,
        xi    = 2+theta^2-pi/8*exp(-theta^2/2)*((2+theta^2)*besseli(0,theta^2/4)+theta^2*besseli(1,theta^2/4))^2;
        g     = sqrt(xi*(1+r^2)-2);
        if abs(theta-g)<1e-6, break; end
        theta = g;
    end
    sig = sqrt(mu2)/sqrt(xi);
    nu  = sqrt(mu1^2+(xi-2)*sig^2);
else
    nu  = 0;
    sig = (2^(1/2)*(mu1^2 + mu2)^(1/2))/2;
end
%_______________________________________________________________________

%_______________________________________________________________________

function p = ricepdf(x,nu,sig2)
% Rician PDF
% p = ricepdf(x,nu,sig2)
% https://en.wikipedia.org/wiki/Rice_distribution#Characterization
p      = x./sig2.*exp(-(x.^2+nu.^2)/(2*sig2));
msk    = find(p>0); % Done this way to prevent division of 0 by Inf
p(msk) = p(msk).*besseli(0,x(msk)*nu/sig2);

