function varargout = spm_SpUtil(varargin)
% Space matrix utilities
% FORMAT varargout = spm_SpUtil(action,varargin)
%
%_______________________________________________________________________
%
% spm_SpUtil is a multi-function function containing various utilities
% for Design matrix and contrast construction and manipulation. In
% general, it accepts design matrices as plain matrices or as space
% structures setup by spm_sp.
%
% Many of the space utilities are computed using an SVD of the design
% matrix. The advantage of using space structures is that the svd of
% the design matrix is stored in the space structure, thereby saving
% unnecessary repeated computation of the SVD. This presents a
% considerable efficiency gain for large design matrices.
%
% Note that when space structures are passed as arguments is is
% assummed that their basic fields are filled in. See spm_sp for
% details of (design) space structures and their manipulation.
%
% ======================================================================
%
% FORMAT i = spm_SpUtil('isCon',x,c)
% Tests whether weight vectors specify contrasts
% x   - Design matrix X, or space structure of X
% c   - contrast matrix (I.e. matrix of contrast weights, contrasts in columns)
%       Must have column dimension matching that of X
%       [defaults to eye(size(X,2)) to test uniqueness of parameter estimates]
% i   - logical row vector indiciating estimability of contrasts in columns of c
%
% A linear combination of the parameter estimates is a contrast if and
% only if the weight vector is in the space spanned by the rows of X.
%
% The algorithm works by regressing the contrast weight vectors using
% design matrix X' (X transposed). Any contrast weight vectors will be
% fitted exactly by this procedure, leaving zero residual. Parameter
% tol is the tolerance applied when searching for zero residuals.
%
% Christensen R (1996)
%	"Plane Answers to Complex Questions"
%	 2nd Ed. Springer-Verlag, New York
%
%                           ----------------
%
% FORMAT i = spm_SpUtil('allCon',x,c)
% Tests whether all weight vectors specify contrasts:
% Same as all(spm_SpUtil('isCon',x,c)).
%
%                           ----------------
%
% FORMAT r = spm_SpUtil('ConR',x,c)
% Assess orthogonality of contrasts (wirit the data)
% x   - Design matrix X, or space structure of X
% c   - contrast matrix (I.e. matrix of contrast weights, contrasts in columns)
%       Must have column dimension matching that of X
%       [defaults to eye(size(X,2)) to test independence of parameter estimates]
% r   - Contrast correlation matrix, of dimension the number of contrasts.
%
% For the general linear model Y = X*B + E, a contrast weight vector c
% defines a contrast c*B. This is estimated by c*b, where b are the
% least squares estimates of B given by b=pinv(X)*Y. Thus, c*b = w*Y,
% where weight vector w is given by w=c*pinv(X); Since the data are
% assummed independent, two contrasts are indpendent if the
% corresponding weight vectors are orthogonal.
%
% r is the matrix of normalised inner products between the weight
% vectors corresponding to the contrasts. For iid E, r is the
% correlation matrix of the contrasts.
%
% The logical matrix ~r will be true for orthogonal pairs of contrasts.
% 
%                           ----------------
%
% FORMAT r = spm_SpUtil('ConO',x,c)
% Assess orthogonality of contrasts (wirit the data)
% x   - Design matrix X, or space structure of X
% c   - contrast matrix (I.e. matrix of contrast weights, contrasts in columns)
%       Must have column dimension matching that of X
%       [defaults to eye(size(X,2)) to test uniqueness of parameter estimates]
% r   - Contrast orthogonality matrix, of dimension the number of contrasts.
%
% This is the same as ~spm_SpUtil('ConR',X,c), but uses a quicker
% algorithm by looking at the orthogonality of the subspaces of the
% design space which are implied by the contrasts:
%       r = abs(c*X'*X*c')<tol
% 
%                           ----------------
%
% FORMAT c = spm_SpUtil('FCon',x,i0)
% Return F-contrast for specified design matrix partition
% x   - Design matrix X, or space structure of X
% i0  - column indices of null hypothesis design matrix
%
% This functionality returns a rank n mxp matrix of contrasts suitable
% for an extra-sum-of-squares F-test comparing the design X, with a
% reduced design. The design matrix for the reduced design is X0 =
% X(:,i0), a reduction of n degrees of freedom.
%
% The algorithm, due to J-B, and derived from Christensen, computes the
% contrasts as an orthonormal basis set for the rows of the
% hypothesised redundant columns of the design matrix, after
% orthogonalisation with respect to X0. For non-unique designs, there
% are a variety of ways to produce equivalent F-contrasts. This method
% produces contrasts with non-zero weights only for the hypothesised
% redundant columns.
% 
%                           ----------------
%
% FORMAT [X1,X0] = spm_SpUtil('cTestSp',X,c)
% Orthogonalised partitioning of design space implied by F-contrast
% x   - Design matrix X, or space structure of X
% c   - contrast matrix (I.e. matrix of contrast weights, contrasts in columns)
%       Must have column dimension matching that of X
% X1  - contrast space - design matrix corresponding according to contrast
%       (orthogonalised wirit X0)
% X0  - matrix reduced according to null hypothesis
%       (of same size as X but rank deficient)
%
% ( Note that unless X0 is reduced to a set of linearely independant   )
% ( vectors, c will only be contained in the null space of X0.  If X0  )
% ( is "reduced", then the "parent" space of c must be reduced as well )
% ( for c to be the actual null space of X0.                           )
%
% This functionality returns a design matrix subpartition whose columns
% span the hypothesised null design space of a given contrast. Note
% that X1 is orthogonal(ised) to X0, reflecting the situation when an
% F-contrast is tested using the extra sum-of-squares principle (when
% the extra distance in the hypothesised null space is measured
% orthogonal to the space of X0).
%
% Note that the null space design matrix will probably not be a simple
% sub-partition of the full design matrix, although the space spanned
% will be the same.
%
%                           ----------------
%
% FORMAT X1 = spm_SpUtil('iTestSp',X,i0)
% x   - Design matrix X, or space structure of X
% i0  - Columns of X that make up X0 - the reduced model (Ho:B1=0)
% X1  - Hypothesised null design space, i.e. that part of X orthogonal to X0
% This offers the same functionality as the 'cTestSp' option, but for
% simple reduced models formed from the columns of X.
%
%                           ----------------
%
% FORMAT [trRV,trRVRV] = spm_SpUtil('trRV',x[,V])
% trace(RV) & trace(RVRV) - used in df calculation
% x      - Design matrix X, or space structure of X
% V      - V matrix [defult eye] (trRV == trRVRV if V==eye, since R idempotent)
% trRV   - trace(R*V),     computed efficiently
% trRVRV - trace(R*V*R*V), computed efficiently
% This uses the Karl's cunning understanding of the trace:
%              (tr(A*B) = sum(sum(A'*B)).
% If the space of X is set, then algorithm uses x.u to avoid extra computation.
%
%                           ----------------
%
% FORMAT [trMV, trMVMV]] = spm_SpUtil('trMV',x[,V])
% trace(MV) & trace(MVMV) if two ouput arguments.
% x      - Design matrix X, or space structure of X
% V      - V matrix [defult eye] (trMV == trMVMV if V==eye, since M idempotent)
% trMV   - trace(M*V),     computed efficiently
% trMVMV - trace(M*V*M*V), computed efficiently
% Again, this uses the Karl's cunning understanding of the trace:
%              (tr(A*B) = sum(sum(A'.*B)).
% If the space of X is set, then algorithm uses x.u to avoid extra computation.
%
%                           ----------------
%
% FORMAT O = spm_SpUtil('BetaRc',x,c)
% Extra sum of squares matrix O for beta's from contrast
% x   - Design matrix X, or space structure of X
% c   - contrast matrix (I.e. matrix of contrast weights, contrasts in columns)
%       Must have column dimension matching that of X
% O   - Matrix such that b'*O*b = extra sum of squares for F-test of contrast c
%
%                           ----------------
%
% FORMAT Mp = spm_SpUtil('Mpc',x,c)
% warning('spm_SpUtil(''Mpc'',... appears to give the wrong answer!') %-**
% Extra sum of squares matrix O for data from contrast
% x   - Design matrix X, or space structure of X
% c   - contrast matrix (I.e. matrix of contrast weights, contrasts in columns)
%       Must have column dimension matching that of X
% Mp  - Matrix such that Y'*Mp*Y = extra sum of squares for F-test of contrast c
%
%                           ----------------
%
% FORMAT MpX1 = spm_SpUtil('MpX1',x,i1)
% warning('spm_SpUtil(''MpX1'',... appears to give the wrong answer!') %-**
% Extra sum of squares matrix for data from X1 partition
% x   - Design matrix X, or space structure of X
% i1  - Columns of X corresponding to X1 partition X = [X1,X0] & with
%       parameters B = [B1;B0]. Ho:B1=0
% MpX1 - Matrix such that Y'*MpX1*Y = extra sum of square
%        (I.e. ResSS(B0) - ResSS(B))
%
%                           ----------------
%
% FORMAT MpX0 = spm_SpUtil('MpX0',x,i0)
% warning('spm_SpUtil(''MpX1'',... appears to give the wrong answer!') %-**
% Extra sum of squares matrix for data from X0 partition
% x   - Design matrix X, or space structure of X
% i0  - Columns of X corresponding to X0 partition X = [X1,X0] & with
%       parameters B = [B1;B0]. Ho:B1=0
% MpX0 - Matrix such that Y'*MpX0*Y = extra sum of square
%        (I.e. ResSS(B0) - ResSS(B))
%
%                           ----------------
%
% FORMAT b = spm_SpUtil('cxpequi',x,c)
% x   - Design matrix X, or space structure of X
% c   - contrast matrix (I.e. matrix of contrast weights, contrasts in columns)
%       Must have column dimension matching that of X
% b   - True is c is a spanning set for space of X
%       (I.e. if contrast and space test the same thing)
%
%                           ----------------
%
% FORMAT [df1,df2] = spm_SpUtil('edf',x,i0,V)
% (effective) df1 and df2 the residual df for the projector onto the
% null space of x' (residual forming projector) and the numerator of
% the F-test where i0 are the columns for the null hypothesis model.
% x   - Design matrix X, or space structure of X
% i0  - Columns of X corresponding to X0 partition X = [X1,X0] & with
%       parameters B = [B1;B0]. Ho:B1=0
% V   - V matrix
%
%
%_______________________________________________________________________
% %W% Andrew Holmes, Jean-Baptiste Poline %E%

%-Format arguments
%-----------------------------------------------------------------------
if nargin==0, error('do what? no arguments given...')
	else, action = varargin{1}; end



switch lower(action), case {'iscon','allcon','conr','cono'}
%=======================================================================
% i = spm_SpUtil('isCon',x,c)
if nargin==1, varargout={[]}; end;

if spm_sp('isspc',varargin{2})		%-Working with a space structure
	if nargin==2, c=eye(size(varargin{2}.X,2)); else, c=varargin{3}; end
	switch lower(action)
	case 'iscon'
		varargout = {    spm_sp('isinspp',varargin{2},c) };
	case 'allcon'
		varargout = {all(spm_sp('isinspp',varargin{2},c))};
	case 'conr'
		if size(c,1) ~= size(varargin{2}.X,2) 
			error('Contrast not of the right size'), end
		%-Compute inner products of data weight vectors
		% (c'b = c'*pinv(X)*Y = w'*Y
		% (=> w*w' = c'*pinv(X)*pinv(X)'*c = c'*pinv(X'*X)*c
		r   = c'*spm_sp('pinvxpx',varargin{2})*c;
		%-normalize by "cov(r)" to get correlations
		r   = r./(sqrt(diag(r))*sqrt(diag(r))');
		r(abs(r)<varargin{2}.tol)=0;		%-set near-zeros to zero
		varargout = {r};			%-return r
	case 'cono'
		%-This is the same as ~spm_SpUtil('ConR',x,c), and so returns
		% the contrast orthogonality (though not their corelations).
		varargout = {abs(c'*varargin{2}.X'*varargin{2}.X*c)...
							<varargin{2}.tol};
	end

else					%-Working with a raw design matrix
	if nargin==2, c=eye(size(varargin{2},2)); else, c=varargin{3}; end
	if size(c,1) ~= size(varargin{2},2)
		error('Contrast not of the right size'), end
	tol=max(size(varargin{2}))*norm(varargin{2})*eps;
	switch lower(action)
	case 'iscon'
		varargout={    all(c-varargin{2}'*pinv(varargin{2}')*c<tol,1) };
	case 'allcon'
		varargout={all(all(c-varargin{2}'*pinv(varargin{2}')*c<tol,1))};
	case 'conr'
		%-Compute inner products of data weight vectors for contrast
		% (c'b = c'*pinv(X)*Y = w'*Y
		% (=> w*w' = c'*pinv(X)*pinv(X)'*c = c'*pinv(X'*X)*c
		r = c'*pinv(varargin{2}'*varargin{2})*c;
		%-normalize by "cov(r)" to get correlations
		r = r./(sqrt(diag(r))*sqrt(diag(r))');
		r(abs(r)<tol)=0;			%-set near-zeros to zero
		varargout = {r};			%-return r
	case 'cono'
		%-This is the same as ~spm_SpUtil('ConR',X,c), and so returns
		% the contrast orthogonality (though not their corelations).
		varargout = {abs(c'*varargin{2}'*varargin{2}*c)<tol};
	end
end % (if spm_sp('isspc'...)



case 'fcon'        %-Compute F-contrast given partition of design matrix
%=======================================================================
% c = spm_SpUtil('FCon',x,i0)

%-Argument checks
%-----------------------------------------------------------------------
if nargin<3, error('insufficient arguments'), else, i0 = varargin{3}; end
bSp = spm_sp('isspc',varargin{2});
if bSp, p=size(varargin{2}.X,2); else, p=size(varargin{2},2); end
if all(ismember(i0,[0,1])) & length(i0(:))==p, i0=find(i0); end
if ~isempty(i0) & any(floor(i0)~=i0) | any(i0<1) | any(i0>p)
	error('logical mask or vector of column indices required')
end

%-Computation
%-----------------------------------------------------------------------
%-Work out compliment of i0 - the columns hypothesised redundant
i1 = setdiff([1:p],i0);
if isempty(i1), varargout = {[]}; warning('empty contrast space'), return, end

if bSp
	if isempty(i0), varargout={varargin{2}.X'*varargin{2}.X}; return, end
	xX0 = spm_sp('set',varargin{2}.X(:,i0));	%-Space of Ho model
	varargout = {varargin{2}.X'*spm_sp('res',xX0,varargin{2}.X(:,i1))};
else
	if isempty(i0), varargout={varargin{2}'*varargin{2}}; return, end
	xX0 = spm_sp('set',varargin{2}(:,i0)); 		%-Space of Ho model
	varargout = {varargin{2}'  *spm_sp('res',xX0,varargin{2}(:,i1))};
end
%-Note : the result can be orthogonolized if necessary


case 'ctestsp'       %-Orthogonalised partitioning implied by F-contrast
%=======================================================================
% [X1,X0] = spm_SpUtil('cTestSp',X,c)

if nargin<3, error('Insufficient arguments'), end
if ~spm_sp('isspc',varargin{3})
	c = spm_sp('set',varargin{3});
else
	c = varargin{3};
end

if spm_sp('isspc',varargin{2})
	if size(c.X,1) ~= size(varargin{2}.X,2)
		error('matrix & contrast dimensions don''t match')
	end
	rk = varargin{2}.rk; if rk==0, error('Rank is null'), end
	
	%X_0 = varargin{2}.X*orth(spm_sp('res',c,varargin{2}.X'));
	X0 = varargin{2}.X*spm_sp('res',c,varargin{2}.v(:,1:rk));
	X1 = varargin{2}.X*c.X - X0*pinv(X0)*varargin{2}.X*c.X;
	varargout = {X1,X0};
else
	if size(c.X,1) ~= size(varargin{2},2)
		error('matrix & contrast dimensions don''t match')
	end
	
	X0 = varargin{2}*spm_sp('res',c,orth(varargin{2}'));
	X1 = varargin{2}*spm_sp('oP',c,orth(varargin{2}'));
	X1 = varargin{2}*c.X - X0*pinv(X0)*varargin{2}*c.X;
	varargout = {X1,X0};
end


case 'itestsp'       %-Returns space tested whilst keeping size of X(i0)
%=======================================================================
% X1 = spm_SpUtil('iTestSp',X,i0)

%-Argument checks
%-----------------------------------------------------------------------
if nargin<3, error('insufficient arguments'), else, i0 = varargin{3}; end
bSp = spm_sp('isspc',varargin{2});
if bSp, p=size(varargin{2}.X,2); else, p=size(varargin{2},2); end
if all(ismember(i0,[0,1])) & length(i0(:))==p, i0=find(i0); end
if ~isempty(i0) & any(floor(i0)~=i0) | any(i0<1) | any(i0>p)
	error('logical mask or vector of column indices required')
end
i1   = setdiff([1:p],i0);	%-Indices of columns hypothesised redundant

%-Computation
%-----------------------------------------------------------------------
if spm_sp('isspc',varargin{2})
	if isempty(i1)
		varargout = {spm_sp('oP',varargin{2})}; 
	else 
		varargout = {varargin{2}.X(:,i1) - ...
			varargin{2}.X(:,i0)*pinv(varargin{2}.X(:,i0))*...
			varargin{2}.X(:,i1)};
	end
else 
	if isempty(i1)
		varargout = {varargin{2}*pinv(varargin{2})}; 
	else 
		varargout = {varargin{2}(:,i1) - ...
			varargin{2}(:,i0)*pinv(varargin{2}(:,i0))*...
			varargin{2}(:,i1)};
	end
end


case 'trrv'                      %-Traces for (effective) df calculation
%=======================================================================
% [trRV,trRVRV]= spm_SpUtil('trRV',x[,V])

if nargin == 1, error('insufficient arguments'), end

if nargin == 2	 			%-no V specified: trRV == trRVRV
	if spm_sp('isspc',varargin{2}) 
		varargout = {size(varargin{2}.X,1) - varargin{2}.rk, ...
			size(varargin{2}.X,1) -varargin{2}.rk};
	else
		rk = rank(varargin{2});
		varargout = {size(varargin{2},1) - rk, ...
			size(varargin{2},1) -rk}; 
	end
   
else 					%-V provided: assume correct!

    if spm_sp('isspc',varargin{2})
	rk = varargin{2}.rk; 
        if rk==0, error('Rank is null'), end

	if nargout==1
		%-only trRV needed
		trMV = sum(sum( varargin{2}.u(:,1:rk)' .* ...
			(varargin{2}.u(:,1:rk)'*varargin{3}) ));
		varargout = { trace(varargin{3}) - trMV};
	else
		%-trRVRV is needed as well
		MV = varargin{2}.u(:,1:rk)*(varargin{2}.u(:,1:rk)'*varargin{3});
		%-NB: It's possible to avoid computing MV; with MV = u*(u'*V) 
		trMV = trace(MV);
		trRVRV = sum(sum(varargin{3}.*varargin{3})) - ...
			2*sum(sum(varargin{3}.*MV)) + ...
			sum(sum(MV'.*MV));
		varargout = {(trace(varargin{3}) - trMV), trRVRV};
	end
    else 
	if nargout==1
		%-only trRV is needed
		trMV = sum(sum(varargin{2}'.*(pinv(varargin{2})*varargin{3})));
		varargout = {trace(varargin{3}) -trMV};
	else
		%-trRVRV is needed as well
		MV     = varargin{2}*(pinv(varargin{2})*varargin{3});
		trMV   = trace(MV);					
		trRVRV = sum(sum(varargin{3}.*varargin{3})) - ...
			2*sum(sum(varargin{3}.*MV)) + ...
			sum(sum(MV'.*MV));
		varargout = {(trace(varargin{3}) -trMV), trRVRV};
	end
    end
end


case 'trmv'                     %-Traces for (effective) Fdf calculation
%=======================================================================
% [trMV, trMVMV]] = spm_SpUtil('trMV',x[,V])

if nargin == 1, error('insufficient arguments'), end

if nargin == 2	 			%-no V specified: trMV == trMVMV
	if spm_sp('isspc',varargin{2})
		varargout = {varargin{2}.rk, varargin{2}.rk};
	else
		rk = rank(varargin{2}); varargout = {rk,rk};
	end
   
else 					%- V provided, and assumed correct !

    if spm_sp('isspc',varargin{2})
	rk = varargin{2}.rk; 
        if rk==0, error('Rank is null'), end

	if nargout==1
		%-only trMV needed
		trMV = sum(sum(varargin{2}.u(:,1:rk)' .* ...
			(varargin{2}.u(:,1:rk)'*varargin{3}) ));
		varargout = {trMV};
	else 
		%-trMVMV is needed as well
		MV = varargin{2}.u(:,1:rk)*(varargin{2}.u(:,1:rk)'*varargin{3});
		%-(See note on saving mem. for MV in 'trRV')
		trMVMV = sum(sum(MV'.*MV));
		varargout = {trace(MV), trMVMV};
	end

    else 
	if nargout==1
		%-only trMV is needed then use : 
		trMV = sum(sum(varargin{2}'.*(pinv(varargin{2})*varargin{3})))
		varargout = {trMV};
	else
		MV     =  varargin{2}*(pinv(varargin{2})*varargin{3});
		trMVMV = sum(sum(MV'.*MV));
		varargout = {trace(MV), trMVMV};
	end
    end
end



case 'betarc'     %-Extra sum of squares matrix for beta's from contrast
%=======================================================================
% O = spm_SpUtil('BetaRc',x,c)

if nargin<3, error('insufficient arguments'), end
c = varargin{3};
if spm_sp('isspc',varargin{2})
	if ~spm_sp('isinspp',varargin{2},varargin{3}), ...
		error('contrast not estimable'), end
	varargout = { c*pinv(c'*spm_sp('pinvXpX',varargin{2})*c)*c' };
	%-(Note there's a better way: first sets the space of c...
else
	varargout = { c*pinv(c'*pinv(varargin{2}'*varargin{2})*c)*c' };
end


case 'mpc'          %-Extra sum of squares matrix for data from contrast
%=======================================================================
% Mp = spm_SpUtil('Mpc',x,c)
warning('spm_SpUtil(''Mpc'',... appears to give the wrong answer!') %-**

if nargin<3, error('insufficient arguments'), end
c = varargin{3};
if spm_sp('isspc',varargin{2})
	if ~spm_sp('isinspp',varargin{2},varargin{3}), 
		error('contrast not estimable'), end
	varargout = {varargin{2}.X * ...
		c*pinv(c'*spm_sp('xpx', varargin{2})*c)*c' * ...
		varargin{2}.X'};
	%-(Note there's a better way: first setting the space of c...
else
	varargout = {varargin{2} * ...
		c*pinv(c'*varargin{2}'*varargin{2}*c)*c' * ...
		varargin{2}'};
end


case 'mpx1'     %-Extra sum of squares matrix for data from X1 partition
%=======================================================================
% MpX1 = spm_SpUtil('MpX1',x,i1)
warning('spm_SpUtil(''MpX1'',... appears to give the wrong answer!') %-**

if nargin<3, error('insufficient arguments'), else, i1 = varargin{3}; end
bSp = spm_sp('isspc',varargin{2});
if bSp, p=size(varargin{2}.X,2); else, p=size(varargin{2},2); end
if all(ismember(i1,[0,1])) & length(i1(:))==p, i1=find(i1); end
if ~isempty(i1) & any(floor(i1)~=i1) | any(i1<1) | any(i1>p)
	error('logical mask or vector of column indices required')
end

if spm_sp('isspc',varargin{2})
	varargout = {varargin{2}.X(:,i1) * ...
		pinv(varargin{2}.X(:,i1)'*varargin{2}.X(:,i1)) * ...
		varargin{2}.X(:,i1)'};
else 
	varargout = {varargin{2}(:,i1) * ...
		pinv(varargin{2}(:,i1)'*varargin{2}(:,i1)) * ...
		varargin{2}(:,i1)'};
end


case 'mpx0'     %-Extra sum of squares matrix for data from X0 partition
%=======================================================================
% MpX0 = spm_SpUtil('MpX0',x,i0)
warning('spm_SpUtil(''MpX0'',... appears to give the wrong answer!') %-**

%-Argument checks
%-----------------------------------------------------------------------
if nargin<3, error('insufficient arguments'), else, i0 = varargin{3}; end
bSp = spm_sp('isspc',varargin{2});
if bSp, p=size(varargin{2}.X,2); else, p=size(varargin{2},2); end
if all(ismember(i0,[0,1])) & length(i0(:))==p, i0=find(i0); end
if ~isempty(i0) & any(floor(i0)~=i0) | any(i0<1) | any(i0>p)
	error('logical mask or vector of column indices required')
end
i1 = setdiff([1:p],i0);		%-Indices of columns hypothesised redundant

%-Computation
%-----------------------------------------------------------------------
if spm_sp('isspc',varargin{2})
	if isempty(i1)
		varargout = {spm_sp('oP',varargin{2})};
	else 
		P = varargin{2}.X(:,i1) - ...
			varargin{2}.X(:,i0)*pinv(varargin{2}.X(:,i0)) * ...
			varargin{2}.X(:,i1);
		varargout = {P*pinv(P)};
	end
else 
	if isempty(i1)
		varargout = {varargin{2}*pinv(varargin{2})};
	else 
		P = varargin{2}(:,i1) - ...
			varargin{2}(:,i0)*pinv(varargin{2}(:,i0)) * ...
			varargin{2}(:,i1);
		varargout = {P*pinv(P)};
	end
end


case 'cxpequi'                       %-Equality of subspace and contrast
%=======================================================================
% b = spm_SpUtil('cxpequi',x,c)

if nargin<3, error('insufficient arguments'), else, c = varargin{3}; end

if spm_sp('isspc',varargin{2})
	varargout = {all(spm_sp('isinspp',varargin{2},c)) & ...
		all(spm_sp('isinsp',spm_sp('set',c), varargin{2}.X'))}; 
else 
	rkc = rank(c);
	if size(c,1) ~= size(varargin{2},2)
		error('dimensions don''t match'), end
	varargout = { 	(rank([varargin{2}' c]) == rkc) & ...
			(rkc == rank(varargin{2}')) };
end


case 'edf'                              %-Effective F degrees of freedom
%=======================================================================
% [df1,df2] = spm_SpUtil('edf',x,i0,V)
%-Argument checks
%-----------------------------------------------------------------------
if nargin<4, error('insufficient arguments'), else, i0 = varargin{3}; end
bSp = spm_sp('isspc',varargin{2});
if bSp, p=size(varargin{2}.X,2); else, p=size(varargin{2},2); end
if all(ismember(i0,[0,1])) & length(i0(:))==p, i0=find(i0); end
if ~isempty(i0) & any(floor(i0)~=i0) | any(i0<1) | any(i0>p)
	error('logical mask or vector of column indices required')
end

[trRV,trRVRV]    = spm_SpUtil('trRV', varargin{2}, varargin{4});
[trMpV,trMpVMpV] = spm_SpUtil('trMV',spm_SpUtil('iTestSp',varargin{2}, i0),...
			varargin{4});
varargout = {trMpV^2/trMpVMpV, trRV^2/trRVRV};


otherwise
%=======================================================================
error('Unknown action string in spm_SpUtil')



%=======================================================================
end
