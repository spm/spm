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

% find independent partitions (encoded in W)
%---------------------------------------------------------------------------
m     = length(Q);
n     = length(Cy);
d     = zeros(m,m);
RQR   = sparse(n,n);
for i = 1:m
	RQ{i} = Q{i} - X*(X'*Q{i});
	RQR   = RQR | RQ{i};
end
for i = 1:m
for j = i:m
	q      = trace(RQ{i}*RQ{j});
	d(i,j) = q;
	d(j,i) = q;
end
end

% recursive calling for (p) separable partitions 
%---------------------------------------------------------------------------
u     = spm_blk_diag(RQR);
h     = zeros(m,1);
W     = zeros(m,m); 
Ce    = sparse(n,n); 
p     = size(u,2);
if p > 1
	for i = 1:p

		% find subset of bases (indexed by J) and partitions (by q)
		%-----------------------------------------------------------
		J     = [];
		C     = {};
		q     = find(u(:,i));
		for j = 1:m
			if nnz(Q{j}(q,q))
				C{end + 1} = Q{j}(q,q);
				J(end + 1) = j;
			end
		end

		% ReML
		%-----------------------------------------------------------
		fprintf('%-30s- %i\n','  ReML Partition',i);
		[Cep,hp,Wp] = spm_reml(Cy(q,q),X(q,:),C,1e-16);
		Ce(q,q)     = Ce(q,q) + Cep;
		h(J)        = hp;
		W(J,J)      = Wp;
		
	end

	return
end

% If one partition proceed to ReML
%===========================================================================

% initialize hyperparameters
%---------------------------------------------------------------------------
for i = 1:m
	h(i) = any(diag(Q{i}));
end
h     = d*pinv(d)*h;

% Iterative EM
%---------------------------------------------------------------------------
dFdh  = zeros(m,1);
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
	dh    = pinv(W)*dFdh;
	h     = h + dh;

	% Convergence (or break if there is only one hyperparameter)
	%===================================================================
	w     = dFdh'*dFdh;
	if w < TOL, break, end
	fprintf('%-30s: %i %30s%e\n','  ReML Iteration',k,'...',full(w));
	if m == 1,   break, end
end
