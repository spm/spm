function [pE,pC] = spm_NESS_priors(n,L,V,W,J,R,C)
% priors for a NESS genertive model
% FORMAT [pE,pC] = spm_NESS_priors(n,L,V,[W,J,R,C])
%--------------------------------------------------------------------------
% n   - number of states (n)
% L   - order (+1) of polynomial expansion
% V   - prior covariance over parameters
% W   - precisions of random fluctuations
% J   - contraints on Jacobian [ones(n,n)]
% R   - coupling to mean (first order) [zeros(n,n)]
% C   - covaraince of NESS [eye(n,n)]
%
% pE  - prior expectation
% pC  - prior covariances (matrix)
% sC  - prior covariances (structure)
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

% checks and defaults
%--------------------------------------------------------------------------
if nargin < 5,  J = ones(n,n);  end
if nargin < 6,  R = zeros(n,n); end
if nargin < 7,  C = eye(n,n);   end

if isscalar(C), C = eye(n,n)*C; end
if isscalar(W), W = eye(n,n)*W; end

if isnumeric(V)
    P.Qp  = V;            % covariance of parameters (solenoidal)
    P.Sp  = V;            % covariance of parameters (surprisal)
    P.Rp  = V;            % covariance of parameters (expectation)
    V     = P;
end

% get (polynomial) expansion: assume first-order approximation
%--------------------------------------------------------------------------
o     = (1:L) - 1;
for i = 2:n
    o = repmat(o,1,L);
    o = [o; kron((1:L) - 1,ones(1,L^(i - 1)))];
end
k     = sum(o) < L;
o     = o(:,k);
nb    = size(o,2);

% get parameters
%--------------------------------------------------------------------------
pE.Qp = zeros(nb,n,n);    % polynomial coefficients for solenoidal operator
pE.Rp = zeros(nb,n);      % polynomial coefficients for surprisal mean
pE.Sp = zeros(nb,n,n);    % polynomial coefficients for surprisal kernel
pE.W  = W;                % precision of random fluctuations

pC.Qp = pE.Qp;
pC.Rp = pE.Rp;
pC.Sp = pE.Sp;
pC.W  = W*0;

% constraints on solenoidal operator (G is modelled by M.W)
%--------------------------------------------------------------------------
for i = 1:n
    for j = 1:n
        if J(i,j)
            d     = ~J(i,:);
            for k = 1:nb

                % if j influences i and upper diagonal part of Q
                %----------------------------------------------------------
                if j > i && ~any(o(d,k))
                    pC.Qp(k,i,j) = V.Qp;
                end

            end
        end
    end
end

% constraints on potential (surprisal)
%--------------------------------------------------------------------------
for i = 1:n
    for j = 1:n
        if J(i,j)

            % disallowed influences
            %--------------------------------------------------------------
            d     = ~J(i,:);
            for k = 1:nb

                % if j influences i 
                %----------------------------------------------------------
                if i == j
                    pE.Sp(1,i,j) = sqrt(1/C(i,j));      % precision of NESS
                end
                if j == i && ~any(o(d,k))   % j >= i for nonorthogonal NESS
                    pC.Sp(k,i,j) = V.Sp;
                end

            end
        end
    end
end




% constraints on mean
%--------------------------------------------------------------------------
pC.Rp(1,:) = V.Rp;

% constraints on mean (coupling to enslaved states)
%--------------------------------------------------------------------------
for i = 1:n
    for j = 1:n
        if R(i,j)
            for k = 1:nb

                % if this basis is first order in j
                %----------------------------------------------------------
                if sum(o(:,k)) == 1 && o(j,k)
                    pC.Rp(k,i) = V.Rp;
                end
            end
        end
    end
end

qE    = [];
qC    = [];
for i = 1:n
    for j = i:n
        for k = 1:nb
            qE(end + 1) = pE.Qp(k,i,j);
            qC(end + 1) = pC.Qp(k,i,j);
        end
    end
end
pE.Qp = qE;
pC.Qp = qC;

return
