function [C,P] = spm_PEB(y,P)
% PEB estimation for hierarchical linear  models - FULL EFFICIENT version
% FORMAT [C,P] = spm_PEB(y,P)
% y       - (n x 1)     response variable
%
% PRIOR SPECIFICATION OF MODEL
% P{i}.X  - (n x m)     ith level design matrix i.e:
%                       contraints on the form of E{b{i - 1}}
% P{i}.C  - {q}(n x n)  ith level contraints on the form of Cov{e{i}} i.e:
%                       contraints on the form of Cov{b{i - 1}}
% c{i}    - {q}         level-specific contrasts
%
% POSTERIOR OR CONDITIONAL ESTIMATES
%
% C{i}.E  - (n x 1)     conditional expectation E{b{i - 1}|y}
% C{i}.C  - (n x n)     conditional covariance  Cov{b{i - 1}|y} = Cov{e{i}|y}
% C{i}.M  - (n x n)     ML estimate of Cov{b{i - 1}} = Cov{e{i}}
% C{i}.h  - (q x 1)     ith level hyperparameters for covariance:
%                       [inv]Cov{e{i}} = P{i}.h(1)*P{i}.C{1} +  ...
%
% HYPERPARAMETER ESTIMATES
%
% P{1}.h  - ReML Hyperparameter estimates
% P{1}.W  - ReML precisions = ddln(p(y|h)/dhdh
%
% If P{i}.C is not a cell the covariance at that level is assumed to be kown
% and Cov{e{i}} = P{i}.C (i.e. the hyperparameter is fixed at 1)
%
% If P{n}.C is not a cell this is taken to indicate that a full Bayesian
% estimate is required where P{n}.X is the prior expectation and P{n}.C is
% the known prior covariance.  For consistency, with PEB, this is implemented
% by setting b{n} = 1 through appropriate constraints at level {n + 1}.
%
% To implement non-hierarchical Bayes with priors on the parameters use
% a two level model setting the second level design matrix to zero.
%___________________________________________________________________________
%
% Returns the moments of the posterior p.d.f. of the parameters of a 
% hierarchical linear observation model under Gaussian assumptions
%
%                            y = X{1}*b{1} + e{1}
%                         b{1} = X{2}*b{2} + e{2}
%                                 ...
%
%                     b{n - 1} = X{n}*b{n} + e{n}
%
% e{n} ~ N{0,Ce{n}} 
%
% using Parametic Emprical Bayes (PEB)
%
% Ref: Dempster A.P., Rubin D.B. and Tsutakawa R.K. (1981) Estimation in
% covariance component models.  J. Am. Stat. Assoc. 76;341-353
%___________________________________________________________________________
% %W% Karl Friston, John Ashburner %E%


% number of levels (p)
%---------------------------------------------------------------------------
M     = 32;				% maximum number of iterations
p     = length(P);

% check covariance constraints - assume i.i.d. errors conforming to X{i}
%---------------------------------------------------------------------------
for i = 1:p
	if ~isfield(P{i},'C')
		[n m] = size(P{i}.X);
		if i == 1
			P{i}.C            = {speye(n,n)};
		else
			for j = 1:m
				k         = find(P{i}.X(:,j));
				P{i}.C{j} = sparse(k,k,1,n,n);
			end
		end
	end

end

% Construct augmented non-hierarchical model
%===========================================================================

% design matrix and indices
%---------------------------------------------------------------------------
I     = {0};
J     = {0};
K     = {0};
XX    = [];
X     = 1;
for i = 1:p

	% design matrix
	%-------------------------------------------------------------------
	X     = X*P{i}.X;
	XX    = [XX X];

	% indices for ith level parameters
	%-------------------------------------------------------------------
	[n m] = size(P{i}.X);
	I{i}  = [1:n] + I{end}(end);
	J{i}  = [1:m] + J{end}(end);

end

% augment design matrix and data
%---------------------------------------------------------------------------
n        = size(XX,2);
XX       = [XX; speye(n,n)];
y        = [y; sparse(n,1)];

% last level constraints 
%---------------------------------------------------------------------------
n        = size(P{p}.X,2);
I{p + 1} = [1:n] + I{end}(end);
q        = I{end}(end);
Cb       = sparse(q,q);
if ~iscell(P{end}.C)

	% Full Bayes: (i.e. Cov(b) = 0, <b> = 1)
	%-------------------------------------------------------------------
	y( I{end})        = sparse(1:n,1,1);
	Cb(I{end},I{end}) = sparse(1:n,1:n,1e-6);
else

	% Empirical Bayes: uniform priors (i.e. Cov(b) = Inf, <b> = 0)
	%-------------------------------------------------------------------
	Cb(I{end},I{end}) = sparse(1:n,1:n,1e+6);
end


% assemble augmented  constraints Q: [inv]Cov{e} = Cb + h(i)*Q{i} + ...
%===========================================================================
if ~isfield(P{1},'Q')

	% covariance contraints Q on Cov{e{i}} = Cov{b{i - 1}}
	%-------------------------------------------------------------------
	h     = [];
	Q     = {};
	for i = 1:p

		% collect constraints on prior covariances - [inv]Cov{e{i}}
		%-----------------------------------------------------------
		if iscell(P{i}.C)
			m     = length(P{i}.C);
			for j = 1:m
				h          = [h; any(diag(P{i}.C{j}))];
				[u v s]    = find(P{i}.C{j});
				u          = u + I{i}(1) - 1;
				v          = v + I{i}(1) - 1;
				Q{end + 1} = sparse(u,v,s,q,q);
			end

			% indices for ith level hyperparameters
			%---------------------------------------------------
			K{i}  = [1:m] + K{end}(end);


		% unless they are known - augment Cb  
		%-----------------------------------------------------------
		else
			[u v s] = find(P{i}.C + speye(length(P{i}.C))*1e-6);
			u       = u + I{i}(1) - 1;
			v       = v + I{i}(1) - 1;
			Cb      = Cb + sparse(u,v,s,q,q);

			% indices for ith level hyperparameters
			%---------------------------------------------------
			K{i}  = [];

		end

	end

	% note overlapping bases - requiring 2nd order M-Step derivatives 
	%-------------------------------------------------------------------
	m     = length(Q);
	d     = sparse(m,m);
	for i = 1:m
		XQX{i} = XX'*Q{i}*XX;
	end
	for i = 1:m
	for j = i:m
		o      = nnz(XQX{i}*XQX{j});
		d(i,j) = o;
		d(j,i) = o;
	end
	end

	P{1}.Cb = Cb;
	P{1}.Q  = Q;
	P{1}.h  = h;
	P{1}.K  = K;
	P{1}.d  = d;

end
Cb    = P{1}.Cb;
Q     = P{1}.Q;
h     = P{1}.h;
K     = P{1}.K;
d     = P{1}.d;


% Iterative EM
%---------------------------------------------------------------------------
m     = length(Q);
dFdh  = zeros(m,1);
W     = zeros(m,m);
for k = 1:M

	% inv(Cov(e)) - iCe(h)
	%-------------------------------------------------------------------
	Ce    = Cb;
	for i = 1:m
		Ce = Ce + h(i)*Q{i};
	end
	iCe   = inv(Ce);


	% E-step: conditional mean E{B|y} and covariance cov(B|y)
	%===================================================================
        iCeX  = iCe*XX;
        Cby   = inv(XX'*iCeX);
	B     = Cby*(iCeX'*y);


	% M-step: ReML estimate of hyperparameters (if m > 0)
	%===================================================================
	if m == 0, break, end

	% Gradient dFd/h (first derivatives)
	%-------------------------------------------------------------------
	Py    = iCe*(y - XX*B);
	iCeXC = iCeX*Cby;
	for i = 1:m

		% dF/dh = -trace(dF/diCe*iCe*Q{i}*iCe) = 
		%-----------------------------------------------------------
		PQ{i}   = iCe*Q{i} - iCeXC*(iCeX'*Q{i});
		dFdh(i) = -trace(PQ{i}) + Py'*Q{i}*Py;
	end

	% Expected curvature E{ddF/dhh} (second derivatives)
	%-------------------------------------------------------------------
	for i = 1:m
		for j = i:m
		if d(i,j)

			% ddF/dhh = -trace{P*Q{i}*P*Q{j}}
			%---------------------------------------------------
			ddFdhh = sum(sum(PQ{i}.*PQ{j}'));
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
	if w < 1e-6, break, end
	fprintf('%-30s: %i %30s%e\n','  PEB Iteration',k,'...',full(w));
end

% place hyperparameters in P{1} and output structure for {n + 1}
%---------------------------------------------------------------------------
P{1}.h             = h;
P{1}.W             = W;
C{p + 1}.E         = B(J{p});
C{p + 1}.M         = Cb(I{end},I{end});

% recursive computation of conditional means E{b|y} 
%---------------------------------------------------------------------------
for i = p:-1:2
        C{i}.E     = B(J{i - 1}) + P{i}.X*C{i + 1}.E;
end

% conditional covariances Cov{b|y} and ReML esimtates of Ce{i) = Cb{i - 1}
%---------------------------------------------------------------------------
for i = 1:p
        C{i + 1}.C = Cby(J{i},J{i});
        C{i}.M     = Ce(I{i},I{i});
        C{i}.h     = h(K{i});
end

% warning
%---------------------------------------------------------------------------
if k == M, warning('maximum number of iterations exceeded'), end
