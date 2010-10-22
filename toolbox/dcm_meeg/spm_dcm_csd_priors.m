function [pE] = spm_dcm_csd_priors(M,U,Y,k)
% Optimisation of priors
% FORMAT [pE] = spm_dcm_csd_priors(M,U,Y,k)
%__________________________________________________________________________
%
% M.IS - function name f(P,M,U) - generative model
%        This function specifies the nonlinear model: 
%        y = Y.y = IS(P,M,U) + X0*P0 + e
%        were e ~ N(0,C).  For dynamic systems this would be an integration
%        scheme (e.g. spm_int). spm_int expects the following:
%
%     M.f  - f(x,u,P,M)
%     M.g  - g(x,u,P,M)
%       x  - state variables
%       u  - inputs or causes
%       P  - free parameters
%       M  - fixed functional forms and parameters in M
%
% M.FS - function name f(y,M)   - feature selection
%        This [optional] function performs feature selection assuming the
%        generalized model y = FS(y,M) = FS(IS(P,M,U),M) + X0*P0 + e
%
% M.P  - starting estimates for model parameters [optional]
%
% M.pE - prior expectation  - E{P}   of model parameters
% M.pC - prior covariance   - Cov{P} of model parameters
%
% M.hE - prior expectation  - E{h}   of log-precision parameters
% M.hC - prior covariance   - Cov{h} of log-precision parameters
%
% U.u  - inputs
% U.dt - sampling interval
%
% Y.y  - outputs
% Y.dt - sampling interval for outputs
% Y.X0 - Confounds or null space      (over size(y,1) bins or all vec(y))
% Y.Q  - q error precision components (over size(y,1) bins or all vec(y))
%
% k    - indices of parameter vector to search over
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_csd_priors.m 4095 2010-10-22 19:37:51Z karl $
 
% reduce system to one source
%==========================================================================
 
% reduce data to s single channel
%--------------------------------------------------------------------------
y     = 0;
for i = 1:length(Y.y)
    for j = 1:size(Y.y{i},2)
        y = y + Y.y{i}(:,j,j);
    end
end
y     = real(y);
 
% reduce model to a single source
%--------------------------------------------------------------------------
S.dipfit    = M.dipfit;
S.dipfit.Ns = 1;
S.dipfit.Nc = 1;
model       = S.dipfit.model;
[pE,pC]     = spm_dcm_neural_priors({1,1,1},[],1,model);
[pE,pC]     = spm_L_priors(S.dipfit,pE,pC);
[pE,pC]     = spm_ssr_priors(pE,pC);
[x,f]       = spm_dcm_x_neural(pE,model);

% single source model
%--------------------------------------------------------------------------
S.IS = 'spm_csd_mtf';
S.g  = 'spm_gx_erp';
S.f  = f;
S.x  = x;
S.n  = length(spm_vec(x));
S.pE = pE;
S.m  = 1;
S.Hz = M.Hz;

% unpack covariance
%--------------------------------------------------------------------------
if isstruct(pC);
    pC = spm_diag(spm_vec(pC));
end
 
% find parameters that most affect expected frequency
%==========================================================================
try
    k;
catch
    dfdp  = spm_diff('spm_csd_chf',S.pE,S,U,1);
    dfdp  = pC*spm_vec(dfdp);
    [q k] = sort(-abs(dfdp));
end
 
% and compute csd over two of them
%--------------------------------------------------------------------------
s     = linspace(-3,3,8);
P     = [];
X     = [];
for i = 1:length(s)
    for j = 1:length(s);
        P(:,end + 1) = spm_vec(S.pE);
        P(k(1),end)  = P(k(1),end) + sqrt(pC(k(1),k(1)))*s(i);
        P(k(2),end)  = P(k(2),end) + sqrt(pC(k(2),k(2)))*s(j);
        [G w]        = spm_csd_mtf(spm_unvec(P(:,end),S.pE),S,U);
        X(:,end + 1) = spm_vec(G);
    end
end
 
% find the best match with average response
%==========================================================================
for i = 1:size(P,2)
    L      = [w.^0 1./w X(:,i)/norm(X(:,i))];
    b(:,i) = (L\y);
    g(:,i) = L*b(:,i);
    ssq(i) = sum((y - g(:,i)).^2);
end
[q i] = min(ssq);
P     = P(:,i);
G     = g(:,i);
B     = b(:,i);
 
subplot(2,1,1)
plot(w,y,w,G,':')
title('optimised priors')
xlabel('frequency')
ylabel('power')
axis square
drawnow
 
% adjust original (neuronal) priors
%--------------------------------------------------------------------------
pE     = M.pE;
pS     = spm_unvec(P,S.pE);
feilds = fieldnames(pS);
for  i = 1:length(feilds)
    try
        pF = getfield(pS,feilds{i});
        pM = getfield(pE,feilds{i});
        if isvector(pM)
            if any(pM - pF)
                pF = pM - pM + pF;
                pE = setfield(pE,feilds{i},pF);
                fprintf('paramters %s optimised\n',feilds{i})
            end
        end
    end
end

 
return
 
 
% show expected frequency as a function of two parameters
%==========================================================================
dfdp  = spm_diff('spm_csd_chf',S.pE,S,U,1);
dfdp  = pC*spm_vec(dfdp);
[q k] = sort(-abs(dfdp));
s     = linspace(-4,4,64);
for i = 1:length(s)
    for j = 1:length(s);
        P        = spm_vec(S.pE);
        P(k(1))  = P(k(1)) + sqrt(pC(k(1),k(1)))*s(i);
        P(k(2))  = P(k(2)) + sqrt(pC(k(2),k(2)))*s(j);
        P        = spm_unvec(P,S.pE);
        Ew(i,j)  = spm_csd_chf(P,S,U);
    end
end
 
subplot(2,2,1)
imagesc(s,s,Ew)
title('expected frequency')
xlabel(sprintf('parameter %i',k(2)))
ylabel(sprintf('parameter %i',k(1)))
axis square
 
subplot(2,2,2)
plot(s,Ew)
title('expected frequency')
xlabel(sprintf('parameter %i',k(2)))
ylabel(sprintf('parameter %i',k(1)))
axis square
