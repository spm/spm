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
% $Id: spm_rice_mixture.m 4922 2012-09-13 12:45:46Z john $

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
        p(:,k) = mg(k)*ricepdf(x(:),nu(k),sig(k)^2);
    end
    sp  = sum(p,2)+eps;
    oll = ll;
    ll  = sum(log(sp).*h(:));
    if ll-oll<1e-8*sum(h), break; end

    %fprintf('%g\n',ll);
    %md = mean(diff(x));
    %plot(x(:),p,'--',x(:),h/sum(h)/md,'b.',x(:),sp,'r'); drawnow

    p = bsxfun(@rdivide,p,sp);


    for k=1:K,
        m0(k) = sum(p(:,k).*h(:));
        m1(k) = sum(p(:,k).*h(:).*x(:));
        m2(k) = sum(p(:,k).*h(:).*x(:).*x(:));
    end
    mg = m0/sum(m0);
    for k=1:K,
        mu1 = m1(k)./m0(k);
        mu2 = (m2(k)-m1(k)*m1(k)/m0(k)+1e-6)/(m0(k)+1e-6);
        [nu(k),sig(k)] = moments2param(mu1,mu2);
    end
    %disp([nu'; sig'])
end
%_______________________________________________________________________

%_______________________________________________________________________

function [nu,sig] = moments2param(mu1,mu2)
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
p      = x./sig2.*exp(-(x.^2+nu.^2)/(2*sig2));
msk    = find(p>0); % Done this way to prevent division of 0 by Inf
p(msk) = p(msk).*besseli(0,x(msk)*nu/sig2);

