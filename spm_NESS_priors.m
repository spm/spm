function [pE,pC] = spm_NESS_priors(n,K,V,W)
% priors for a NESS genertive model
% FORMAT [pE,pC] = spm_NESS_priors(n,K,V,W)
%--------------------------------------------------------------------------
% n   = 3;                    % number of states
% K   = 1 + 1;                % order (+1) of polynomial expansion (suprisal)
% V   = 512                   % prior precision over parameters
% W   = 32;                   % precision of random fluctuations
%
% pE  - prior expectation
% pC  - prior covariances
%
% This routine returns the prior expectations and covariances of the
% parameters (i.e., polynomial coefficients) of a nonequilibrium
% steady-state model. This model is based upon the Helmholtz-Hodge
% decomposition of the solution to density dynamics of any system that
% possesses a pullback attractor. The requisite flow operators Q and
% surprisal S are parameterised to 2nd and fourth order, respectively. The
% surprisal is parameterised with even terms to 4th order by creating a
% surprisal kernel K and then taking the outer product H = K’*K, where the
% kernel K is parameterised to first-order. The solenoidal
% (non-dissipative) parts of the flow operator are parameterised to 2nd
% order in the states.
%
% It is assumed that the dissipative part of the flow (i.e., the amplitude
% of random fluctuations G) is spherical and constant. This corresponds to
% the inverse precision of the random fluctuations that can also be
% specified as M.W in the accompanying model M.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging

% assume second-order approximation (for solenoidal flow)
%--------------------------------------------------------------------------
L     = 3;

% get (polynomial) expansion: L > K = 3
%--------------------------------------------------------------------------
o     = (1:L) - 1;
for i = 2:n
    o = repmat(o,1,L);
    o = [o; kron((1:L) - 1,ones(1,L^(i - 1)))];
end
k  = sum(o) < L;
o  = o(:,k);
    
% get parameters
%--------------------------------------------------------------------------
nb     = size(o,2);
nQ     = n*n/2 + n/2;

pE.Qp  = zeros(nb*nQ,1);  % polynomial coefficients for solenoidal operator
pE.Rp  = zeros(1,n);      % polynomial coefficients for surprisal mean
pE.Sp  = zeros(nb,n,n);   % polynomial coefficients for surprisal kernel
pE.W   = W;               % precision of random fluctuations


pC.Qp  = pE.Qp;
pC.Rp  = pE.Rp;
pC.Sp  = pE.Sp;
pC.W   = 0;

% constraints on solenoidal operator (G is modelled by M.W)
%--------------------------------------------------------------------------
[ks,kq]  = spm_NESS_indices(o,K);
k        = find(kq);
pC.Qp(k) = 1/V;

% constraints on mean
%--------------------------------------------------------------------------
k        = 1:n;
pC.Rp(k) = 0; %%%

% constraints on potential (surprisal)
%--------------------------------------------------------------------------
k     = find(ks);
for i = 1:n
    for j = 1:n
        if i == j
           pE.Sp(1,i,j) = 1; %%% exp(-j);
        end
        if j == i   %%%
           pC.Sp(k,i,j) = 1/V;
        end
    end
end

% priors
%--------------------------------------------------------------------------
pC    = diag(spm_vec(pC));

return

function [ks,kq,kg] = spm_NESS_indices(o,K)
% constraints on polynomial coefficients of dynamical systems
% FORMAT [ks,kq,kg] = spm_NESS_indices(o,K);
% o - matrix of orders for polynomial expansion
% K - upper bound on order for surprisal parameters
%
% ks  - indices for surprisal   parameters
% kq  - indices for solenoidal  parameters
% kg  - indices for dissipative parameters
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% constraints on potential parameters due to dissipative flow
%==========================================================================
[n,nb] = size(o);                           % number of basis functions

% terms for suprisal kernel
%--------------------------------------------------------------------------
ks    = sum(o) < K;                         % polynomial order constraints

% solenoidal part of Q (i.e., R)
%--------------------------------------------------------------------------
kq    = [];
for i = 1:n
    for j = i:n
        if j == i
            kq = [kq false(1,nb)];
        else
            kq = [kq  true(1,nb)];
        end
    end
end

% dissipative part of Q (i.e., G, which is modelled by W)
%--------------------------------------------------------------------------
kg    = [];
for i = 1:n
    for j = i:n
        if j == i
            kg = [kg  true(1,nb)];
        else
            kg = [kg false(1,nb)];
        end
    end
end

return

