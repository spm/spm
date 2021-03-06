function [sig0,alph] = spm_smohist(t0,lam,alph)
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

% John Ashburner
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


sig0 = zeros(size(t0));
n  = size(t0,1);
if nargin<2
    lam = zeros(size(t0,2),1);
    x   = (1:n)';
    for k=1:size(t0,2)
        t  = t0(:,k)+eps;
        mu = sum(t.*x)./sum(t);
        vr = sum(t.*(x-mu).^2)/sum(t);
        lam(k) = vr;
    end
end

if nargin<3
    alph = log(t0+1);
end

% Regularisation
G0          = spdiags(repmat([-1 2 -1],n,1),[-1 0 1],n,n);
G0(1,1)     = 1;
G0(end,end) = 1;
G0          = G0'*G0;

for k=1:size(t0,2)
    %fprintf('\n%d) ');
    t   = t0(:,k);
    G   = G0*lam(k);
    am  = alph(:,k);
    mx  = max(am);
    sig = exp(am-mx);
    sig = sig/sum(sig);
    L   = spdiags(ones(n,1)*(sum(t) + lam(k))*1e-6,0,n,n);
    ll  = Inf;
    st  = sum(t);
    for it=1:60
        gr  = sum(t)*sig - t + G*am;
        W   = spdiags(sig*sum(t) + abs(gr)*1e-8,0,n,n);
        H   = W + G + L;
        da  = H\gr;
        am  = am - da;

        % Softmax function
        mx  = max(am);
        sig = exp(am-mx);
        lse = log(sum(sig)) + mx;
        sig = sig/sum(sig);
        if ~rem(it,4)
            oll = ll;
           %ll  = -sum(log(sig+1e-9).*t) + 0.5*am'*G*am;
            ll  = -(am'*t - sum(t)*lse)  + 0.5*am'*G*am;
            if oll-ll<st*1e-5, break; end
        end
   end
   alph(:,k) = am;
   sig0(:,k) = sig;
end
