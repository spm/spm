function [Ce,h,W] = spm_reml(Cy,X,Q,TOL);
% REML estimation of covariance components from Cov{y}
% FORMAT [Ce,h,W] = spm_reml(Cy,X,Q,TOL);
%
% Cy  - (m x m) data covariance matrix y*y'  {y = (m x n) data matrix}
% X   - (m x p) design matrix
% Q   - {1 x q} covariance constraints
% TOL - Tolerance {default = 1e-6}
%
% Ce  - (m x m) estimated errors = h(1)*Q{1} + h(2)*Q{2} + ...
% h   - (q x 1) hyperparameters
% W   - (q x q) W*n = precision of hyperparameter estimates 
%___________________________________________________________________________
% %W% John Ashburner, Karl Friston %E%

% set tolerance if not specified
%---------------------------------------------------------------------------
if nargin < 4, TOL = 1e-6; end

% ensure X is not rank deficient
%---------------------------------------------------------------------------
X     = spm_svd(X);

% find independent partitions (encoded in d)
%---------------------------------------------------------------------------
m     = length(Q);
d     = sparse(m,m);
for i = 1:m
	XQX{i} = X'*Q{i}*X;
end
for i = 1:m
for j = i:m
	o      = nnz(XQX{i}*XQX{j});
	d(i,j) = o;
	d(j,i) = o;
end
end

% recursive calling for (p) seperable partitions 
%---------------------------------------------------------------------------
h     = zeros(m,1);
n     = length(Cy);
Ce    = sparse(n,n); 
u     = spm_svd(d);
p     = size(u,2);
if p > 1

	% cycle over (p) partitions
	% j = indices of {1 x b} Q  for ith partition
	% q = indices of (n x n) Cy for ith partition
	%-------------------------------------------------------------------
	for i = 1:p

		fprintf('%-30s: %i %30s\n','  ReML Partition',i,'...');

		% find subset of bases (indexed by j) and partitions (by q)
		%-----------------------------------------------------------
		j     = find(u(:,i));
		b     = length(j);
		C     = {};
		q     = sparse(1,n);
		for k = 1:b
			q = q + any(Q{j(k)});
		end
		q     = find(q);
		for k = 1:b
			C{k} = Q{j(k)}(q,q);
		end

		% ReML
		%-----------------------------------------------------------
		[Cep,hp] = spm_reml(Cy(q,q),X(q,:),C);
		Ce(q,q)  = Ce(q,q) + Cep;
		h(j)     = hp;
	end

	return
end

% If one partition proceed to ReML
%===========================================================================

% initialize hyperparamters (h)
%---------------------------------------------------------------------------
for i = 1:m
	h(i)   = any(diag(Q{i}));
end

% Iterative EM
%---------------------------------------------------------------------------
dFdh  = zeros(m,1);
W     = zeros(m,m);
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
		if d(i,j)

			% ddF/dhh = -trace{P*Q{i}*P*Q{j}}
			%---------------------------------------------------
			ddFdhh = sum(sum(PQ{i}.*PQ{j}'))/2;
			W(i,j) = ddFdhh;
			W(j,i) = ddFdhh;
		end
		end
	end

	% Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
	%-------------------------------------------------------------------
	dh    = W\dFdh;
	h     = h + dh;

	% Convergence (or break if there is only one hyperparameter)
	%===================================================================
	w     = dFdh'*dFdh;
	if w < TOL, break, end
	fprintf('%-30s: %i %30s%e\n','  ReML Iteration',k,'...',full(w));
	if m == 1,   break, end
end
