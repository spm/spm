function varargout = spm_DesUtil(varargin)
% Design matrix utilities
% FORMAT varargout = spm_DesUtil(action,varargin)
%
%_______________________________________________________________________
%
% spm_DesUtil is a multi-function function containing various utilities
% for Design matrix construction and manipulation
%
% ======================================================================
%
% FORMAT i = spm_DesUtil('pds',v,m,n)
%
% Patterned data setting function - inspired by MINITAB's "SET" command
% v - base pattern vector
% m - (scalar natural number) #replications of elements of v [default 1]
% n - (scalar natural number) #repeats of pattern [default 1]
% i - resultant pattern vector, with v's elements replicated m times,
%     the resulting vector repeated n times.
%
%                           ----------------
%
% spm_DesUtil('pds'... is a simple utility for patterned data setting,
% inspired by MINITAB's "SET" command. It is particularly useful for
% creating structured indicator vectors.
%
% The vector v has it's elements replicated m times, and the resulting
% vector is itself repeated n times, giving a resultant vector i of
% length n*m*length(v)
%
% Examples:
%     spm_DesUtil('pds',1:3)       % = [1,2,3]
%     spm_DesUtil('pds',1:3,2)     % = [1,1,2,2,3,3]
%     spm_DesUtil('pds',1:3,2,3)   % = [1,1,2,2,3,3,1,1,2,2,3,3,1,1,2,2,3,3]
% NB: MINITAB's "SET" command has syntax n(v)m:
%
% ======================================================================
%
% FORMAT i = spm_DesUtil('isCon',X,c,tol)
% Tests whether weight vectors specify contrasts
% X   - design matrix
% c   - contrast matrix (I.e. matrix of contrast weights, contrasts in rows)
%       Must have row dimension matching that of X
%       [defaults to eye(size(X,2)) to test uniqueness of parameter estimates]
% tol - Tolerance for computation [default max(size(X))*norm(X)*eps, as in rank]
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
% ======================================================================
%
% FORMAT i = spm_DesUtil('allCon',X,c)
% Tests whether all weight vectors specify contrasts:
% Same as all(spm_DesUtil('isCon',X,c)), but works directly.
%
% Since a contrast must have weight vector contained in the space
% spanned by the rows of the design matrix X, [X;c] will have the same
% rank as X if all the weight vectors of c specify contrasts.
%
% ======================================================================
%
% FORMAT r = spm_DesUtil('ConR',X,c,tol)
% Assess orthogonality of contrasts (wirit the data)
% X   - design matrix
% c   - contrast matrix (I.e. matrix of contrast weights, contrasts in rows)
%       Must have row dimension matching that of X
%       [defaults to eye(size(X,2)) to test uniqueness of parameter estimates]
% tol - Tolerance for computation [default max(size(X))*norm(X)*eps, as in rank]
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
% ======================================================================
%
% FORMAT r = spm_DesUtil('ConO',X,c,tol)
% Assess orthogonality of contrasts (wirit the data)
% X   - design matrix
% c   - contrast matrix (I.e. matrix of contrast weights, contrasts in rows)
%       Must have row dimension matching that of X
%       [defaults to eye(size(X,2)) to test uniqueness of parameter estimates]
% tol - Tolerance for computation [default max(size(X))*norm(X)*eps, as in rank]
% r   - Contrast orthogonality matrix, of dimension the number of contrasts.
%
% This is the same as ~spm_DesUtil('ConR',X,c), but uses a quicker
% algorithm by looking at the orthogonality of the subspaces of the
% design space which are implied by the contrasts:
%       r = abs(c*X'*X*c')<tol
% 
% ======================================================================
%
% FORMAT c = spm_DesUtil('FCon',X,i0)
% Return F-contrast for specified design matrix partition
% X   - design matrix
% i0  - column indices of null hypothesis design matrix
%
% This functionality returns a nxp matrix of contrasts suitable for an
% extra-sum-of-squares F-test comparing the design X, with a reduced
% design. The design matrix for the reduced design is X0 = X(:,i0), a
% reduction of n degrees of freedom.
%
% The algorithm, due to J-B, and derived from Christensen, computes the
% contrasts as an orthonormal basis set for the rows of the
% hypothesised redundant columns of the design matrix, after
% orthogonalisation with respect to X0. For non-unique designs, there
% are a variety of ways to produce equivalent F-contrasts. This method
% produces contrasts with non-zero weights only for the hypothesised
% redundant columns.
% 
% ======================================================================
%
% FORMAT [X0,X1] = spm_DesUtil('cNullSpace',X,c)
% Contrast null space for a design matrix and contrast
% X   - design matrix
% c   - contrast
% X0  - null design space - matrix reduced according to null hypothesis
% X1  - contrast space - design matrix corresponding according to contrast
%       (orthogonalised wirit X0)
%
% This functionality returns a design matrix subpartition whose columns
% span the null space of a given contrast.
%
% The algorithm (andrew's), computes the compliment (null space)
% contrast as an orthonormal basis set for the rows of the design
% matrix orthogonalised with respect to the contrast c, and then
% multiplying the columns of the design matrix by this null space
% contrast. The resulting mini design matrix has columns spanning the
% null design space implied by the contrast.
%
% Note that the null space design matrix will probably not be a simple
% sub-partition of the full design matrix, althouth the space spanned
% will be the same.
%_______________________________________________________________________
% %W% Andrew Holmes, Jean-Baptiste Poline %E%

%-Format arguments
%-----------------------------------------------------------------------
if nargin==0, error('do what? no arguments given...')
	else, action = varargin{1}; end



switch lower(action), case 'pds'
%=======================================================================
%-Condition arguments
%-----------------------------------------------------------------------
% i = spm_DesUtil('pds',v,m,n)
if nargin<4, n=1; else, n=varargin{4}; end
if nargin<3, m=1; else, m=varargin{3}; end
if nargin<2, varargout={[]}, return, else, v=varargin{2}; end
if any([size(n),size(m)])>1, error('n & m must be scalars'), end
if any(([m,n]~=floor([m,n]))|([m,n]<1))
	error('n & m must be natural numbers'), end
if sum(size(v)>1)>1, error('v must be a vector'), end

%-Computation
%-----------------------------------------------------------------------
si = ones(1,ndims(v)); si(find(size(v)>1))=n*m*length(v);
varargout = {reshape(repmat(v(:)',m,n),si)};



case {'iscon','allcon','conr','cono'}
%=======================================================================
% i = spm_DesUtil('isCon',X,c,tol)
if nargin<2, X=[]; else, X=varargin{2}; end
if isempty(X), varargout={[]}; return, end
if nargin<3, c=eye(size(X,2)); else, c=varargin{3}; end
if nargin<4, tol=max(size(X))*norm(X)*eps; else, tol=varargin{4}; end

switch lower(action), case 'iscon'
	varargout = {all(abs(c'-X'*pinv(X')*c')<tol,1)'};
case 'allcon'
	varargout = { rank(X)==rank([X;c]) };
case 'conr'
	w   = c*pinv(X);		%-contrast vectors for data
	r   = w*w';			%-inner products of data weight vectors
	tmp = diag(r); tmp = sqrt(tmp)*sqrt(tmp)';
	r   = r./tmp;			%-normalise r
	r(abs(r)<tol)=0;		%-set near-zeros to zero
	varargout = {r};		%-return r
case 'cono'
	%-This is the same as ~spm_DesUtil('ConR',X,c), and so returns
	% the contrast orthogonality (though not their corelations).
	varargout = {abs(c*X'*X*c')<tol};
end


case 'fcon'
%=======================================================================
% c = spm_DesUtil('FCon',X,i0)
if nargin<3, error('insufficient arguments'), end
X  = varargin{2};
i0 = varargin{3};

%-Argument checks
%-----------------------------------------------------------------------
p = size(X,2);
if any(floor(i0)~=i0) | any(i0<1) | any(i0>p)
	error('vector of column indices required')
end

%-Computation
%-----------------------------------------------------------------------
i1   = setdiff([1:p],i0);	%-Indices of columns hypothesised redundant
X0   = X(:,i0);			%-Design matrix of hypothesised null model
X1   = X(:,i1);			%-Design matrix columns hypothesised redundant

X1o  = X1 - X0*pinv(X0)*X1;	%-X1 columns orthogonalised wirit X0

oX1o = orth(X1o);		%-Orthonormal basis for columns of X1o
				% space spanned by this basis set is the
				% design space hypothesised redundant

c1   = orth(X1o');		%-Orthonormal basis set for rows of X1o
				%-These are the contrast coefficients for X1

%c1   = orth(X1o'*pinv(X1o')*eye(length(i1)));	%-JB's version
				%-Orthonormal basis set for identity matrix
				% projected onto rows of X1o
				%-These are the contrast coefficients for X1

c       = zeros(size(X,2),size(c1,2));
c(i1,:) = c1;			%-Form full contrast matrix
c       = c';			%-Transpose to usual SPM contrast orientation
				% with contrasts as row vectors

varargout = {c};


case 'cnullspace'
%=======================================================================
% X0 = spm_DesUtil('cNullSpace',X,c)
if nargin<3, error('Insufficient arguments'), end
X = varargin{2};
c = varargin{3};

%-Computation
%-----------------------------------------------------------------------
c0 = orth(X' -c'*pinv(c')*X')';		%-Compliment of current contrast
					% (Null contrast space)

X0 = X*c0';				%-Design space corresponding to null
					% contrast space

X1 = orth(X*c' -X0*pinv(X0)*(X*c'));	%-Contrast design space
					% (orthogonalised wirit X0)
varargout = {X0,X1};


otherwise
%=======================================================================
error('Unknown action string')



%=======================================================================
end
