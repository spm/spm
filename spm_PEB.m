function [C,P] = spm_PEB(y,P)
% PEB estimation for hierarchical linear observation models
% FORMAT [C,P] = spm_PEB(y,P)
% y       - (n x 1)     response variable
%
% PRIOR SPECIFICATION OF MODEL
% P{i}.X  - (n x m)     ith level design matrix i.e:
%                       prior contraints on the form of E{b{i - 1}}
% P{i}.C  - {q}(n x n)  ith level contraints on the form of Cov{e{i}} i.e:
%                       prior contraints on the form of Cov{b{i - 1}}
%
% POSTERIOR OR CONDITIONAL ESTIMATES
%
% C{i}.E  - (n x 1)     conditional expectation E{b{i - 1}|y}
% C{i}.C  - (n x n)     conditional covariance  Cov{b{i - 1}|y} = Cov{e{i}|y}
% C{i}.M  - (n x n)     ML estimate of Cov{b{i - 1}} = Cov{e{i}}
% C{i}.h  - (q x 1)     ith level hyperparameters for covariance:
%                       Cov{e{i}} = P{i}.h(1)*P{i}.C{1} +  ...
%
% if P{n}.C is not a cell the last structure is taken to specify the
% expectation and known covariance of the penultimate level
%___________________________________________________________________________
%
% Returns the moments of the posterior p.d.f. of the parameters of a 
% hierarchical linear observation model under Gaussian assumptions
%
%                           y    = X{1}*b{1} + e{1}
%                           b{1} = X{2}*b{2} + e{2}
%                                     ...
%
%                           b{n} = X{n}*b{n} + e{n}
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

% If full Bayes - only estimate hyperparameters upto the penultimate level
%---------------------------------------------------------------------------
if ~iscell(P{end}.C)
	p = p - 1;
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


% last level constraints 
%---------------------------------------------------------------------------
n        = size(P{p}.X,2);
I{p + 1} = [1:n] + I{end}(end);
if ~iscell(P{end}.C)

	% Full Bayes: remove expectation P{n}.X from y
	%-------------------------------------------------------------------
	y                 = y - X*P{end}.X;
	C{p + 1}.E        = P{end}.X;

	% and set Cb to P{n}.C    (i.e. Cov(b) = P{n}.C, <b> = P{n}.X)
	%-------------------------------------------------------------------
	Cp                = P{end}.C + sparse(1:n,1:n,1e-8,n,n);
	C{p + 1}.M        = Cp;
else

	% Empirical Bayes: uniform priors (i.e. Cov(b) = Inf, <b> = 0)
	%-------------------------------------------------------------------
	Cp                = sparse(1:n,1:n,1e8);
	C{p + 1}.E        = sparse(1:n,1,0);
	C{p + 1}.M        = Cp;
end

% put prior covariance in Cb
%---------------------------------------------------------------------------
q        = I{end}(end);
Cb       = sparse(q,q);
Cb(I{end},I{end}) = Cp;


% augment design matrix and data
%---------------------------------------------------------------------------
n        = size(XX,2);
XX       = [XX; speye(n,n)];
y        = [y; sparse(n,1)];


% assemble augmented constraints Q 
%===========================================================================
if ~isfield(P{1},'iQ')

	% covariance contraints Q on Cov{e{i}} = Cov{b{i - 1}}
	%-------------------------------------------------------------------
	h     = [];
	Q     = {};
	for i = 1:p

		% collect constraints on prior covariances - Cov{e{i}}
		%-----------------------------------------------------------
		m     = length(P{i}.C);
		for j = 1:m
			h                 = [h; any(diag(P{i}.C{j}))];
			[u v s]           = find(P{i}.C{j});
			u                 = u + I{i}(1) - 1;
			v                 = v + I{i}(1) - 1;
			Q{end + 1}        = sparse(u,v,s,q,q);
		end

		% indices for ith level hyperparameters
		%-----------------------------------------------------------
		K{i}  = [1:m] + K{end}(end);
	end


	% overlapping bases
	%-------------------------------------------------------------------
	m     = length(Q);
	for i = 1:m
	for j = i:m
		d(i,j)  = nnz(Q{i}*Q{j});
	end
	end

	P{1}.Q  = Q;
	P{1}.h  = h;
	P{1}.K  = K;
	P{1}.d  = d;

end
Q     = P{1}.Q;
h     = P{1}.h;
K     = P{1}.K;
d     = P{1}.d;


% Iterative EM estimation
%---------------------------------------------------------------------------
m     = length(Q);
dFdh  = zeros(m,1);
W     = zeros(m,m);
for k = 1:M

	% assemble estimate of inv(Cov(e)) - iCe
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


	% M-step: ML estimate of hyperparameters
	%===================================================================

	% Ce:  using dF/diCe = r*r' + X*Cby*X' and diCe/dh(i) = iCe*Q{i}*iCe
	%-------------------------------------------------------------------
	r     = y - XX*B;
	Cr    = iCe*r;
	for i = 1:m

		% dF/dh = trace(dF/diCe*iCe*Q{i}*iCe) = 
		%----------------------------------------------------------
		QC{i}   = Q{i}*iCe;
		dFdh(i) = trace(QC{i}) - Cr'*Q{i}*Cr - ...
			  sum(sum(Cby.*(iCeX'*Q{i}*iCeX)));
	end
	for i = 1:m
		for j = i:m
		if d(i,j)

			% -ddF/dhh = trace{Q{i}*iCe*Q{j}*iCe}
			%--------------------------------------------------
			ddFdhh = sum(sum(QC{i}.*QC{j}));
			W(i,j) = ddFdhh;
			W(j,i) = ddFdhh;

		end
		end
	end

	% Upadate dh = -ddF/dhh*dF/dh
	%-------------------------------------------------------------------
	dh    = W\dFdh;
	h     = h - dh;

	% Convergence
	%===================================================================
	w     = dh'*dh;
	fprintf('%-30s: %i %30s%e\n','  PEB Iteration',k,'...',full(w));
	if w < 1e-6, break, end

end

% place hyperparameters in P{1}
%---------------------------------------------------------------------------
P{1}.h     = h;

% conditional means E{b|y}
%---------------------------------------------------------------------------
B          = Cby*iCeX'*y;
C{p + 1}.E = B(J{p}) + C{p + 1}.E;
for i = p:-1:2
        C{i}.E     = B(J{i - 1}) + P{i}.X*C{i + 1}.E;
end

% conditional covariances Cov{b|y} and ML esimtates of Ce{i) = Cb{i - 1}
%---------------------------------------------------------------------------
for i = 1:p
        C{i + 1}.C = Cby(J{i},J{i});
        C{i}.M     = Ce(I{i},I{i});
        C{i}.h     = h(K{i});
end


% warning
%---------------------------------------------------------------------------
if k == M, warning('maximum number of iterations exceeded'), end

