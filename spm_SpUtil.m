function varargout = spm_DesUtil(varargin)
% Design matrix utilities
% FORMAT varargout = spm_DesUtil(action,varargin)
%
%_______________________________________________________________________
%
% spm_DesUtil is a multi-function function containing various utilities
% for Design matrix construction and manipulation
%
%=======================================================================
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
%=======================================================================
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
%=======================================================================
% FORMAT i = spm_DesUtil('allCon',X,c)
% Tests whether all weight vectors specify contrasts:
% Same as all(spm_DesUtil('isCon',X,c)), but works directly.
%
% Since a contrast must have weight vector contained in the space
% spanned by the rows of the design matrix X, [X;c] will have the same
% rank as X if all the weight vectors of c specify contrasts.
%
%=======================================================================
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
%_______________________________________________________________________
% %E% Andrew Holmes %W%

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



case {'iscon','allcon','conr'}
%=======================================================================
% i = spm_DesUtil('isCon',X,c,tol)
if nargin<2, X=[]; else, X=varargin{2}; end
if isempty(X), varargout={[]}; return, end
if nargin<3, c=eye(size(X,2)); else, c=varargin{3}; end
if nargin<4, tol=max(size(X))*norm(X)*eps; else, tol=varargin{4}; end

switch lower(action), case 'iscon'
	varargout = {all(abs(c'-X'*pinv(X')*c')<tol,1)};
case 'allcon'
	varargout = { rank(X)==rank([X;c]) };
case 'conr'
	w   = c*pinv(X);		%-contrast vectors for data
	r   = w*w';			%-inner products of data weight vectors
	tmp = diag(r); tmp = sqrt(tmp)*sqrt(tmp)';
	r   = r./tmp;			%-normalise r
	r(abs(r)<tol)=0;		%-set near-zeros to zero
	varargout = {r};		%-return r
end



otherwise
%=======================================================================
error('Unknown action string')



%=======================================================================
end
