function varargout = spm_sp(varargin)
% Orthogonal (design) matrix space setting & manipulation
% FORMAT varargout = spm_spc(action,varargin)
%
% This function computes the different projectors related to the row
% and column spaces X. It should be used to avoid redundant computation
% of svd on large X matrix.  It is divided into actions that set up the
% space, (Create,Set,...) and actions that compute projections (pinv,
% pinvXpX, pinvXXp, ...) This is motivated by the problem of rounding
% errors that can invalidate some computation and is a tool to work
% with spaces.
%
% The only thing that is not easily computed is the null space of
% the line of X (assuming size(X,1) > size(X,2)).
% To get this space (a basis of it or a projector on it) use spm_sp on X'.
%
% The only restriction on the use of the space structure is when X is
% so big that you can't fit X and its svd in memory at the same time.
% Otherwise, the use of spm_sp will generally speed up computations and
% optimise memory use.
%
% Note that since the design matrix is stored in the space structure,
% there is no need to keep a separate copy of it.
%
%                           ----------------
%
% The structure is:
%	x = struct(...
%		'X',	[],...		% Mtx
%		'tol',	[],...		% tolerance
%		'ds',	[],...		% vectors of singular values 
%		'u',	[],...		% u as in X = u*diag(ds)*v'
%		'v',	[],...		% v as in X = u*diag(ds)*v'
%		'rk',	[],...		% rank
%		'oP',	[],...		% orthogonal projector on X
%		'oPp',	[],...		% orthogonal projector on X'
%		'ups',	[],...		% space in which this one is embeded
%		'sus',	[]);		% subspace
%
% The basic required fields are X, tol, ds, u, v, rk.
%
% ======================================================================
%
% FORMAT x = spm_sp('Set',X)
% Set up space structure, storing matrix, singular values, rank & tolerance
% X - a (design) matrix (2D)
% x - the corresponding space structure, with basic fields filled in
%     The SVD is an "economy size" svd, using MatLab's svd(X,0)
%
%
% FORMAT r = spm_sp('oP',x[,Y])
% FORMAT r = spm_sp('oPp',x[,Y])
% Return orthogonal projectors, or orthogonal projection of data Y (if passed)
% x - space structure of matrix X
% r - ('oP' usage)  ortho. projection matrix projecting into column space of x.X
%   - ('oPp' usage) ortho. projection matrix projecting into row space of x.X
% Y - data (optional)
%   - If data are specified then the corresponding projection of data is
%     returned. This is usually more efficient that computing and applying
%     the projection matrix directly.
%
%
% FORMAT pX = spm_sp('pinv',x)
% Returns a pseudo-inverse of X - pinv(X) - computed efficiently
% x - space structure of matrix X
% pX - pseudo-inverse of X
% This is the same as MatLab's pinv - the Moore-Penrose pseudoinverse
% ( Note that because size(pinv(X)) == size(X'), it is not generally  )
% ( useful to compute pinv(X)*Data sequentially (as is the case for   )
% ( 'res' or 'oP')                                                    )
%
%
% FORMAT pXpX = spm_sp('pinvxpx',x)
% Returns a pseudo-inverse of X'X - pinv(X'*X) - computed efficiently
% x    - space structure of matrix X
% pXpX - pseudo-inverse of (X'X)
% ( Note that because size(pinv(X'*X)) == [size(X,2) size(X,2)],      )
% ( it is not useful to compute pinv(X'X)*Data sequentially unless    )
% ( size(X,1) < size(X,2)                                             )
%
%
% FORMAT XpX = spm_sp('xpx',x)
% Returns (X'X) - computed efficiently
% x    - space structure of matrix X
% XpX  - (X'X)
%
%
% FORMAT pXXp = spm_sp('pinvxxp',x)
% Returns a pseudo-inverse of XX' - pinv(X*X') - computed efficiently
% x    - space structure of matrix X
% pXXp - pseudo-inverse of (XX')
%
%
% FORMAT XXp = spm_sp('xxp',x)
% Returns (XX') - computed efficiently
% x    - space structure of matrix X
% XXp  - (XX')
%
%
% FORMAT b = spm_sp('isinsp',x,c[,tol])
% FORMAT b = spm_sp('isinspp',x,c[,tol])
% Check whether vectors c are in the column/row space of X
% x   - space structure of matrix X
% c   - vector(s) (Multiple vectors passed as a matrix)
% tol - (optional) tolerance (for rounding error)
%        [defaults to tolerance specified in space structure: x.tol]
% b   - ('isinsp'  usage) true if c is in the column space of X
%     - ('isinspp' usage) true if c is in the column space of X
% 
%
% FORMAT N = spm_sp('null',x)
% Simply returns the null space of matrix X (same as matlab NULL)
% (Null space = vectors associated with zero eigenvalues)
% x - space structure of matrix X
% N - null space
%
%
% FORMAT r = spm_sp('opnull',x[,Y])
% Orthogonal projector onto null space of X, or projection of data Y (if passed)
% x - space structure of matrix X
% Y - (optional) data
% r - (if no Y passed) orthogonal projection matrix  into the null space of X
%   - (if Y passed   ) orthogonal projection of data into the null space of X
% ( Note that if xp = spm_sp('set',x.X'), we have:                    )
% (       spm_sp('opnull',x) == spm_sp('res',xp)                      )
% ( or, equivalently:                                                 )
% (       spm_sp('opnull',x) + spm_sp('oP',xp) == eye(size(xp.X,1));  )
%
%
% FORMAT r = spm_sp('res',x[,Y])
% Returns residual formaing matrix wirit column space of X, or residuals (if Y)
% x - space structure of matrix X
% Y - (optional) data
% r - (if no Y passed) residual forming matrix for design matrix X
%   - (if Y passed   ) residuals, i.e. residual forming matrix times data
%                    ( This will be more efficient than
%                    ( spm_sp('res',x)*Data, when size(X,1) > size(X,2)
% Note that this can also be seen as the orthogonal projector onto the
% null space of x.X' (which is not generally computed in svd, unless
% size(X,1) < size(X,2)).
%
%
% FORMAT oX  = spm_sp('ox', x)
% FORMAT oXp = spm_sp('oxp',x)
% Returns an orthonormal basis for X ('ox' usage) or X' ('oxp' usage)
% x   - space structure of matrix X
% oX  - orthonormal basis for X - same as orth(x.X)
% xOp - *an* orthonormal for X' (but not the same as orth(x.X'))
%
%
% FORMAT b = spm_sp('isspc',x)
% Check a variable is a structure with the right fields for a space structure
% x - candidate variable
% b - true if x is a structure with fieldnames corresponding to spm_sp('create')
%
%
% FORMAT [b,e] = spm_sp('issetspc',x)
% Test whether a variable is a space structure with the basic fields set
% x - candidate variable
% b - true is x is a structure with fieldnames corresponding to
%     spm_sp('Create'), which has it's basic fields filled in.
% e - string describing why x fails the issetspc test (if it does)
% This is simply a gateway function combining spm_sp('isspc',x) with
% the internal subfunction sf_isset, which checks that the basic fields
% are not empty. See sf_isset (below).
%
%-----------------------------------------------------------------------
% SUBFUNCTIONS:
%
% FORMAT b = sf_isset(x)
% Checks that the basic fields are non-empty (doesn't check they're right!)
% x - space structure
% b - true if the basic fields are non-empty
%_______________________________________________________________________
% %W% Jean-Baptiste Poline %E%

if nargin==0
	error('Do what? no arguments given...')
else
	action = varargin{1};
end


switch lower(action), case 'create'             %-Create space structure
%=======================================================================
% x = spm_sp('Create')
varargout = {struct(...
	'X',	[],...		% Matrix
	'tol',	[],...		% tolerance
	'ds',	[],...		% vectors of singular values 
	'u',	[],...		% u as in X = u*diag(ds)*v'
	'v',	[], ...		% v as in X = u*diag(ds)*v'
	'rk',	[],...		% rank
	'oP', 	[],...		% orthogonal projector on X
	'oPp',	[],...		% orthogonal projector on X'
	'ups',	[],...		% space in which this one is embeded
	'sus',	[])};		% subspace 


case 'set'          %-Set singular values, space basis, rank & tolerance
%=======================================================================
% x = spm_sp('Set',X)

if nargin==1 error('No design matrix : can''t do much!'), end
if isempty(varargin{2}), varargout = {spm_sp('create')}; return, end

x   = spm_sp('create');
x.X = varargin{2};

%-- Compute the svd with svd(X,0) : find all the singular values of x.X
%-- SVD(FULL(A)) will usually perform better than SVDS(A,MIN(SIZE(A)))

[x.u, s, x.v] = svd(full(x.X),0);

if size(x.X,1)==1, x.ds = s; else, x.ds = diag(s); end
clear s

%-- compute the tolerance
x.tol =  max(size(x.X))*max(abs(x.ds))*eps;

%-- compute the rank
x.rk =  sum(x.ds > x.tol);

varargout = {x};	


case {'op', 'opp'}                               %-Orthogonal projectors
%=======================================================================
% r = spm_sp('oP', x[,Data])   
% r = spm_sp('oPp',x[,Data])   

%-Check arguments
%-----------------------------------------------------------------------
if nargin==1, error('No structure : can''t do much!'), end
[ok,str] = spm_sp('issetspc',varargin{2}); if ~ok, error(str), end

%-Compute
%-----------------------------------------------------------------------
rk = varargin{2}.rk;
if rk > 0 
  if nargin==2	
    switch lower(action), 
    case 'op'
      varargout = {varargin{2}.u(:,[1:rk])*varargin{2}.u(:,[1:rk])'};
    case 'opp'
      varargout = {varargin{2}.v(:,[1:rk])*varargin{2}.v(:,[1:rk])'};	
    end
  else % nargin==2
     switch lower(action), 
    case 'op'
      if size(varargin{3},1) ~= size(varargin{2}.X,1),
	  error('Dimension dont match'); end 
      varargout = {varargin{2}.u(:,[1:rk])* ...
		(varargin{2}.u(:,[1:rk])' * varargin{3})};

    case 'opp'
      if size(varargin{3},1) ~= size(varargin{2}.X,2),
	  error('Dimension dont match'); end 
      varargout = {varargin{2}.v(:,[1:rk])* ...
		(varargin{2}.v(:,[1:rk])' * varargin{3})};
    end
  end % nargin==2

else varargout = {0}; end; %if rk > 0


case 'pinv'                              %-Pseudo-inverse of X - pinv(X)
%=======================================================================
% pX = spm_sp('pinv',x) 
if nargin==1, error('No structure : can''t do much!'), end
[ok,str] = spm_sp('issetspc',varargin{2}); if ~ok, error(str), end

p   = max(size(varargin{2}.ds));
ds1 = zeros(size(varargin{2}.ds)); 
ds1([1:varargin{2}.rk]) = varargin{2}.ds([1:varargin{2}.rk]).^(-1);
if p == size(varargin{2}.v,1)
	varargout = { varargin{2}.v*diag(ds1)*varargin{2}.u' };
else
	%- then p < size(varargin{2}.v,1)
	varargout = { varargin{2}.v(:,1:p)*diag(ds1)*varargin{2}.u' };
end


case 'pinvxpx'                    %-Pseudo-inverse of (X'X) - pinv(X'*X)
%=======================================================================
% pXpX = spm_sp('pinvxpx',x) 
if nargin==1, error('No structure : can''t do much!'), end
[ok,str] = spm_sp('issetspc',varargin{2}); if ~ok, error(str), end

p   = max(size(varargin{2}.ds));
ds1 = zeros(size(varargin{2}.ds)); 
ds1([1:varargin{2}.rk]) = varargin{2}.ds([1:varargin{2}.rk]).^(-2);
if p == size(varargin{2}.v,1)
	varargout = { varargin{2}.v*diag(ds1)*varargin{2}.v' };
else
	%- then p < size(varargin{2}.v,1)
	varargout = { varargin{2}.v(:,1:p)*diag(ds1)*varargin{2}.v(:,1:p)' };
end


case 'xpx'                                       %-Computation of (X'*X)
%=======================================================================
% XpX = spm_sp('xpx',x)
if nargin==1, error('No structure : can''t do much!'), end
[ok,str] = spm_sp('issetspc',varargin{2}); if ~ok, error(str), end

p   = max(size(varargin{2}.ds));
ds1 = zeros(size(varargin{2}.ds)); 
ds1([1:varargin{2}.rk]) = varargin{2}.ds([1:varargin{2}.rk]).^2;
if p == size(varargin{2}.v,1)
	varargout = { varargin{2}.v*diag(ds1)*varargin{2}.v' };
else
	%- then p < size(varargin{2}.v,1)
	varargout = { varargin{2}.v(:,1:p)*diag(ds1)*varargin{2}.v(:,1:p)' };
end


case 'pinvxxp'                    %-Pseudo-inverse of (XX') - pinv(X*X')
%=======================================================================
% pXXp = spm_sp('pinvxxp',x)
if nargin==1, error('No structure : can''t do much!'), end
[ok,str] = spm_sp('issetspc',varargin{2}); if ~ok, error(str), end

ds1 = zeros(size(varargin{2}.ds)); 
ds1([1:varargin{2}.rk]) = varargin{2}.ds([1:varargin{2}.rk]).^(-2);
varargout = { varargin{2}.u*diag(ds1)*varargin{2}.u' }; clear ds1;
	

case 'xxp'                                       %-Computation of (X*X')
%=======================================================================
% XXp = spm_sp('xxp',x)
if nargin==1, error('No structure : can''t do much!'), end
[ok,str] = spm_sp('issetspc',varargin{2}); if ~ok, error(str), end

ds1 = zeros(size(varargin{2}.ds)); 
ds1([1:varargin{2}.rk]) = varargin{2}.ds([1:varargin{2}.rk]).^2;
varargout = { varargin{2}.u*diag(ds1)*varargin{2}.u' }; clear ds1;
	

case {'isinsp', 'isinspp'}
                    %-Check whether vectors are in row/column space of X
%=======================================================================
% b = spm_sp('isinsp',x,c[,tol])
% b = spm_sp('isinspp',x,c[,tol])

%-Check arguments
%-----------------------------------------------------------------------
if nargin<3, error('insufficient arguments - action,x,c required'), end
[ok,str] = spm_sp('issetspc',varargin{2}); if ~ok, error(str), end
if nargin<4, tol=varargin{2}.tol; else, tol = varargin{4}; end

%-Compute according to case
%-----------------------------------------------------------------------
switch lower(action)
case 'isinsp'
	%-Check dimensions
	if size(varargin{2}.X,1) ~= size(varargin{3},1) 
		error('Vector dimensions don''t match column dimension...'); end
	if ~isempty(varargin{2}.oP)
		varargout = {all((varargin{2}.oP*varargin{3} - ...
							varargin{3})<tol)};
	else
		varargout = {all((spm_sp('oP',varargin{2},varargin{3}) - ...
							varargin{3})<tol )};
	end
case 'isinspp'
	%- check dimensions
	if size(varargin{2}.X,2) ~= size(varargin{3},1) 
		error('Vector dimensions don''t match X row dimension...'); end
	if ~isempty(varargin{2}.oPp)
		varargout = {all((varargin{2}.oPp*varargin{3} - ...
							varargin{3})<tol)};
	else
		varargout = {all((spm_sp('oPp',varargin{2},varargin{3}) - ...
							varargin{3})<tol)};
	end
end


case {'null', 'opnull'}    %-Null space / projector(ion) into null space
%=======================================================================
% N = spm_sp('null',x)
% r = spm_sp('opnull',x[,Y])
if nargin==1, error('No structure : can''t do much!'), end
[ok,str] = spm_sp('issetspc',varargin{2}); if ~ok, error(str), end

rk = varargin{2}.rk;
p  = size(varargin{2}.X,2);
if p==rk | rk == 0, varargout = {[]}; return, end

switch lower(action)
case 'null'
	% null space = vectors associated with 0 eigenvalues
	if nargin==3 error('too many input arguments'), end
	varargout = {varargin{2}.v(:,rk+1:p)};

case 'opnull'
	if nargin==3	% Apply to arg proj. on the 0-space of x.X 
			% check consistency of third argument
		if size(varargin{3},1) ~= size(varargin{2}.X,2)
			error('Inconsistant arg size(Y,1)~=size(x.X,2)'), end
		varargout = {varargin{2}.v(:,rk+1:p)*varargin{2}.v(:,rk+1:p)'...
			*varargin{3} }; 
	else
		varargout = {varargin{2}.v(:,rk+1:p)*varargin{2}.v(:,rk+1:p)'};
		% or varargout = {eye(size(varargin{2}.oPp))-varargin{2}.oPp};
	end
end


case 'res'                        %-Residual formaing matrix / residuals
%=======================================================================
% r = spm_sp('res',x[,Y])
if nargin==1, error('No structure : can''t do much!'), end
[ok,str] = spm_sp('issetspc',varargin{2}); if ~ok, error(str), end

rk = varargin{2}.rk;
if rk > 0 
    switch nargin
    case 2
	if isempty(varargin{2}.oP)
		%-avoid storing varargin{2}.oP = spm_sp('oP',varargin{2});
		varargout = { eye(size(varargin{2}.X,1)) - ...
			varargin{2}.u(:,[1:rk])*varargin{2}.u(:,[1:rk])' };
	 else
	 	varargout = {eye(size(varargin{2}.oP))-varargin{2}.oP};
	 end
	
    case 3
	if size(varargin{3},1) ~= size(varargin{2}.X,1) 
		error('Data and space dim. dont match '); 
	 else 
		varargout = {varargin{3} - varargin{2}.u(:,[1:rk])* ...
			    	(varargin{2}.u(:,[1:rk])'*varargin{3}) };
	 end
	
    end % (switch nargin)

else
	varargout = { 0 };
end


case {'ox', 'oxp'}                    %-Orthonormal basis sets for X / X'
%=======================================================================
% oX  = spm_sp('ox', x)
% oXp = spm_sp('oxp',x)
if nargin==1, error('No structure : can''t do much!'), end
[ok,str] = spm_sp('issetspc',varargin{2}); if ~ok, error(str), end

rk = varargin{2}.rk;
if rk > 0 
	switch lower(action)
	case 'ox'
		varargout = {varargin{2}.u(:,[1:rk])};
	case 'oxp'
		varargout = {varargin{2}.v(:,[1:rk])};
	end
else
	varargout = {0};
end


case 'isspc'                                     %-Space structure check
%=======================================================================
% b = spm_sp('isspc',x)
if nargin~=2, error('too few/many input arguments - need 2'), end

%-Check we've been passed a structure
if ~isstruct(varargin{2}), varargout={0}; return, end

%-Go through required field names checking their existance
% (Get fieldnames once and compare: isfield doesn't work for multiple )
% (fields, and repeated calls to isfield would result in unnecessary  )
% (replicate calls to fieldnames(varargin{2}).                        )
b       = 1;
fnames  = fieldnames(varargin{2});
for str = fieldnames(spm_sp('Create'))'
	b = b & any(strcmp(str,fnames));
	if ~b, break, end
end
varargout = {b};


case 'issetspc'                   %-Is this a completed space structure?
%=======================================================================
% [b,e] = spm_sp('issetspc',x)
if nargin~=2, error('too few/many input arguments - need 2'), end
if ~spm_sp('isspc',varargin{2})
	varargout = {0,'not a space structure (wrong fieldnames)'};
elseif ~sf_isset(varargin{2})
	%-Basic fields aren't filled in
	varargout = {0,'space not defined (use ''set'')'};
else
	varargout = {1,'OK!'};
end


otherwise
%=======================================================================
error(['Invalid action (',action,')'])

%=======================================================================
end % (case lower(action))


%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================

function b = sf_isset(x)
%=======================================================================
b = ~(	isempty(x.X) 	|...
	isempty(x.u) 	|...
	isempty(x.v) 	|...
	isempty(x.ds)	|...
	isempty(x.tol)	|...
	isempty(x.rk)	);
