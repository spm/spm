function [Ce,h,W,u] = spm_reml(Cy,X,Q,TOL,TOS);
% REML estimation of covariance components from Cov{y}
% FORMAT [Ce,h,W,u] = spm_reml(Cy,X,Q,TOL,TOS);
%
% Cy  - (m x m) data covariance matrix y*y'  {y = (m x n) data matrix}
% X   - (m x p) design matrix
% Q   - {1 x q} covariance components
% TOL - Tolerance for convergence [norm of gradient df/dh; default = 1e-6]
% TOS - Tolerance for SVD of curvature [ddf/dhdh;          default = 1e-6]
%
% Ce  - (m x m) estimated errors = h(1)*Q{1} + h(2)*Q{2} + ...
% h   - (q x 1) hyperparameters
% W   - (q x q) W*n = precision of hyperparameter estimates
% u   - {1 x p} estimable components C{i} = u(1,i)*Q{1} + u(2,i)*Q{2} +...
%___________________________________________________________________________
% %W% John Ashburner, Karl Friston %E%

% set tolerance if not specified
%---------------------------------------------------------------------------
if nargin < 4, TOL = 1e-6; end
if nargin < 5, TOS = 1e-6; end

% ensure X is not rank deficient
%---------------------------------------------------------------------------
X     = full(X);
X     = orth(X);
X     = sparse(X);

% find estimable components (encoded in the precision matrix W)
%---------------------------------------------------------------------------
m     = length(Q);
n     = length(Cy);
W     = zeros(m,m);
for i = 1:m
	RQ{i}  = Q{i} - X*(X'*Q{i});
end
for i = 1:m
for j = i:m
	dFdhh  = sum(sum(RQ{i}.*RQ{j}'));
	W(i,j) = dFdhh;
	W(j,i) = dFdhh;
end
end

% eliminate inestimable components
% NB: The threshold for normalized eigenvalues is 1e-6 in spm_svd
%---------------------------------------------------------------------------
u     = spm_svd(W,TOS)
for i = 1:size(u,2)
	C{i}  = sparse(n,n);
	for j = 1:m
		C{i} = C{i} + Q{j}*u(j,i);
	end
end
Q     = C;

% initialize hyperparameters (assuming Cov{e} = 1}
%---------------------------------------------------------------------------
m     = length(Q);
dFdh  = zeros(m,1);
W     = zeros(m,m); 
C     = [];
for i = 1:m
	C = [C Q{i}(:)];
end
I     = speye(n,n);
h     = inv(C'*C)*(C'*I(:));

% Iterative EM
%---------------------------------------------------------------------------
for k = 1:32

	% Q are variance components		
	%------------------------------------------------------------------
	Ce    = sparse(n,n);
	for i = 1:m
		Ce = Ce + h(i)*Q{i};
	end
	iCe   = inv(Ce);

	% E-step: conditional covariance cov(B|y) {Cby}
	%===================================================================
        iCeX  = iCe*X;
        Cby   = inv(X'*iCeX);

	% M-step: ReML estimate of hyperparameters 
	%===================================================================

	% Gradient dFd/h (first derivatives)
	%-------------------------------------------------------------------
	P     = iCe  - iCeX*Cby*iCeX';
	PCy   = Cy*P'- speye(n,n);
	for i = 1:m

		% dF/dh = -trace(dF/diCe*iCe*Q{i}*iCe) = 
		%---------------------------------------------------
		PQ{i}   = P*Q{i};
		dFdh(i) = sum(sum(PCy.*PQ{i}))/2;

	end

	% Expected curvature E{ddF/dhh} (second derivatives)
	%-------------------------------------------------------------------
	for i = 1:m
		for j = i:m

			% ddF/dhh = -trace{P*Q{i}*P*Q{j}}
			%---------------------------------------------------
			dFdhh  = sum(sum(PQ{i}.*PQ{j}))/2;
			W(i,j) = dFdhh;
			W(j,i) = dFdhh;
		end
	end

	% Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
	%-------------------------------------------------------------------
	dh    = inv(W)*dFdh;
	h     = h + dh;

	% Convergence (or break if there is only one hyperparameter)
	%===================================================================
	w     = dFdh'*dFdh;
	if w < TOL | m == 1, break, end
	fprintf('%-30s: %i %30s%e\n','  ReML Iteration',k,'...',full(w));
end

% estimate of cov{e}
%---------------------------------------------------------------------------
Ce    = sparse(n,n);
for i = 1:m
	Ce = Ce + h(i)*Q{i};
end

% rotate hyperparameter esimates and precision back
%---------------------------------------------------------------------------
h     = u*h;
W     = u*W*u';
