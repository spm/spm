function [sig0,alph] = spm_smohist(t0,lam)
% Smooth a histogram
% FORMAT [sig,alpha] = spm_smohist(t,lam)
% t     - a column vector, or matrix of column vectors containing
%         histogram data.
% lam   - regularisation parameter, or vector of regularisation
%         parameters for each column of t.
% sig   - smoothed probability density function(s)
%         (columns sum to 1).
% alpha - logarithm of sig.
%__________________________________________________________________________
%
% Maximises: -sum(log(sig).*t) + 0.5*a'*G*a
% where: sig = exp(a)/sum(exp(a))
%   and: G = lam*L'*L - where L is the Laplacian operator.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_smohist.m 3864 2010-05-05 17:21:20Z john $

sig0 = zeros(size(t0));
n  = size(t0,1);
if nargin<2,
    lam = zeros(size(t0,2),1);
    x   = (1:n)';
    for k=1:size(t0,2),
        t  = t0(:,k)+eps;
        mu = sum(t.*x)./sum(t);
        vr = sum(t.*(x-mu).^2)/sum(t);
        lam(k) = vr;
    end
end

alph = log(t0+1);

% Regularisation
G0 = spdiags(repmat([-1 2 -1],n,1),[-1 0 1],n,n);
G0(1,1)     = 1;
G0(end,end) = 1;
G0          = G0'*G0;

% Prevent over/underflow
constr      = log(realmax)-1;

for k=1:size(t0,2),
    t   = t0(:,k) + 2*sum(t0(:,k))/realmax + eps;
    G   = G0*lam(k);
    am  = alph(:,k);
    sig = exp(am);
    sig = sig/sum(sig);
    L   = spdiags(ones(n,1)*(sum(t) + lam(k))*1e-8,0,n,n);
    ll  = Inf;
    for it=1:60,
        gr  = sum(t)*sig - t + G*am;
        W   = spdiags(sig*sum(t) + abs(gr)*1e-8,0,n,n);
        H   = W + G + L;
        da  = H\gr;

        % Try to prevent over/underflow
        sc  = 1;
        an  = am - da;
        while any(abs(an)>constr)
            sc = sc*0.5;
            an = am - sc*da;
        end
        am  = an - (max(an)+min(an))/2;
        am  = max(am,-constr);
        am  = min(am, constr);

        % Softmax function
        sig = exp(am);
        sig = sig/sum(sig);

        if ~rem(it,4)
            oll = ll;
            ll  = -sum(log(sig+1e-6).*t) + 0.5*am'*G*am;
            if oll-ll<n*1e-8, break; end
        end
       %if ~rem(it,1), plot(1:n,t/sum(t),'g.',1:n,sig,'k-'); drawnow; end
       %fprintf('%d\t%g\t%g\n', it, gr'*gr, ll);
   end
   alph(:,k) = am;
   sig0(:,k) = sig;
end

