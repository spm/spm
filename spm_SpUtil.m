function varargout = spm_SpUtil(varargin)
% Space matrix utilities
% FORMAT varargout = spm_SpUtil(action,varargin)
%
%_______________________________________________________________________
%
% spm_SpUtil is a multi-function function containing various utilities
% for Design matrix construction and manipulation;
% In general, it accepts space structure or plain matrices. See 
% spm_sp for space structure.
%
% ======================================================================
%
% FORMAT i = spm_SpUtil('pds',v,m,n)
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
% spm_SpUtil('pds'... is a simple utility for patterned data setting,
% inspired by MINITAB's "SET" command. It is particularly useful for
% creating structured indicator vectors.
%
% The vector v has it's elements replicated m times, and the resulting
% vector is itself repeated n times, giving a resultant vector i of
% length n*m*length(v)
%
% Examples:
%     spm_SpUtil('pds',1:3)       % = [1,2,3]
%     spm_SpUtil('pds',1:3,2)     % = [1,1,2,2,3,3]
%     spm_SpUtil('pds',1:3,2,3)   % = [1,1,2,2,3,3,1,1,2,2,3,3,1,1,2,2,3,3]
% NB: MINITAB's "SET" command has syntax n(v)m:
%
% ======================================================================
%
% FORMAT i = spm_SpUtil('isCon',X,c,tol)
% Tests whether weight vectors specify contrasts
% X   - Space design matrix
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
% FORMAT i = spm_SpUtil('allCon',X,c)
% Tests whether all weight vectors specify contrasts:
% Same as all(spm_SpUtil('isCon',X,c)), but works directly.
%
% Since a contrast must have weight vector contained in the space
% spanned by the rows of the design matrix X, [X;c] will have the same
% rank as X if all the weight vectors of c specify contrasts.
%
% ======================================================================
%
% FORMAT r = spm_SpUtil('ConR',X,c,tol)
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
% FORMAT r = spm_SpUtil('ConO',X,c,tol)
% Assess orthogonality of contrasts (wirit the data)
% X   - design matrix
% c   - contrast matrix (I.e. matrix of contrast weights, contrasts in rows)
%       Must have row dimension matching that of X
%       [defaults to eye(size(X,2)) to test uniqueness of parameter estimates]
% tol - Tolerance for computation [default max(size(X))*norm(X)*eps, as in rank]
% r   - Contrast orthogonality matrix, of dimension the number of contrasts.
%
% This is the same as ~spm_SpUtil('ConR',X,c), but uses a quicker
% algorithm by looking at the orthogonality of the subspaces of the
% design space which are implied by the contrasts:
%       r = abs(c*X'*X*c')<tol
% 
% ======================================================================
%
% FORMAT c = spm_SpUtil('FCon',X,i0)
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
% FORMAT [X0,X1] = spm_SpUtil('cTestSp',X,c)
% Contrast test space for a design matrix and contrast
% X   - design matrix
% c   - contrast
% X1  - contrast space - design matrix corresponding according to contrast
%       (orthogonalised wirit X0)
% X0  - matrix reduced according to null hypothesis
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
% i = spm_SpUtil('pds',v,m,n)
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
% i = spm_SpUtil('isCon',x,c,tol)
if nargin==1, varargout={[]}; end;
if nargin==2, c=eye(size(varargin{2},2)); else, c=varargin{3}; end
if nargin>=4, tol=varargin{4}; else tol = []; end

switch lower(action), 
case 'iscon'
   if spm_sp('isspc',varargin{2})
	% if size(c,1) ~= size(varargin{2}.X,2) 
	%	error('Contrast not of the right size'); end; %-done in spm_sp
	varargout = { spm_sp('isinspp',varargin{2},c)};
   else
	if size(c,1) ~= size(varargin{2},2) 
		error('Contrast not of the right size'); end; 
	if isempty(tol)
	    tol = max(size(varargin{2})) * norm(varargin{2}) * eps; end;
	varargout = {  (varargin{2}*pinv(varargin{2})*c - c) < tol };
   end
case 'allcon'
   if spm_sp('isspc',varargin{2})
	% if size(c,1) ~= size(varargin{2}.X,2) 
	%	error('Contrast not of the right size'); end; %-done in spm_sp
	varargout = { all(spm_sp('isinspp',varargin{2},c))};
   else
	if size(c,1) ~= size(varargin{2},2) 
		error('Contrast not of the right size'); end; 
	if isempty(tol)
	    tol = max(size(varargin{2})) * norm(varargin{2}) * eps; end;
	varargout = { all( (varargin{2}*pinv(varargin{2})*c - c) < tol )};
   end

case 'conr'
   if spm_sp('isspc',varargin{2})
	if size(c,1) ~= size(varargin{2}.X,2) 
		error('Contrast not of the right size'); end; 
	r   = c'*spm_sp('pinvxpx',varargin{2})*c;	
					%-inner products of data weight vectors
	r   = r./(sqrt(diag(r))*sqrt(diag(r))');
					%-normalize by "cov(r)" to get
					%-correlations
	r(abs(r)<varargin{2}.tol)=0;	%-set near-zeros to zero
	varargout = {r};		%-return r
   else 
	if size(c,1) ~= size(varargin{2},2) 
		error('Contrast not of the right size'); end; 
	if isempty(tol)
	    tol = max(size(varargin{2})) * norm(varargin{2}) * eps; end;
	r = c'*pinv(varargin{2}'*varargin{2})*c;
	r   = r./(sqrt(diag(r))*sqrt(diag(r))');
					%-normalize by "cov(r)" to get 
					%-correlations
	r(abs(r)<tol)=0;		%-set near-zeros to zero
	varargout = {r};		%-return r
   end

case 'cono'
	%-This is the same as ~spm_SpUtil('ConR',X,c), and so returns
	% the contrast orthogonality (though not their corelations).
   if spm_sp('isspc',varargin{2})
	varargout = {abs(c'*varargin{2}.X'*varargin{2}.X*c)<tol};
   else 
	if isempty(tol)
	    tol = max(size(varargin{2})) * norm(varargin{2}) * eps; end;
	varargout = {abs(c'*varargin{2}'*varargin{2}*c)<tol};
   end;

end % switch lower(action), case {'iscon','allcon','conr','cono'}



case 'fcon'
%=======================================================================
% c = spm_SpUtil('FCon',x,i0)
if nargin<3, error('insufficient arguments'), end
i0 = varargin{3};

%-Argument checks and Computation
%-----------------------------------------------------------------------
if spm_sp('isspc',varargin{2}) 
   if isempty(i0) varargout = {varargin{2}.X'*varargin{2}.X}; return; end;
   p = size(varargin{2}.X,2);
   if any(floor(i0)~=i0) | any(i0<1) | any(i0>p)
	error('vector of column indices required'); end
   if max(size(i0))==p & all( ismember(i0,[0 1]) ) i0 = find(i0); end;
   i1   = setdiff([1:p],i0);	%-Indices of columns hypothesised redundant
   if isempty(i1),   varargout = {0}; 
      disp('Warning : dont you want to test for something ? '); return; 
   end;
   x0   = spm_sp('set',varargin{2}.X(:,i0)); 
				%- x0 = Space of hypothesised null model
   varargout = {varargin{2}.X'*spm_sp('res',x0,varargin{2}.X(:,i1))};
   % Note : the result can be orthogonolized if necessary
else
   if isempty(i0) varargout = {varargin{2}'*varargin{2}}; return; end;
   p = size(varargin{2},2);
   if any(floor(i0)~=i0) | any(i0<1) | any(i0>p)
	error('vector of column indices required'); end
   if max(size(i0))==p & all( ismember(i0,[0 1]) ) i0 = find(i0); end;
   i1   = setdiff([1:p],i0);	%-Indices of columns hypothesised redundant
   if isempty(i1),   varargout = {0}; 
      disp('Warning : dont you want to test for something ? '); return; 
   end;
   x0   = spm_sp('set',varargin{2}(:,i0)); 
				%- x0 = Space of hypothesised null model
   varargout = {varargin{2}'*spm_sp('res',x0,varargin{2}(:,i1))};
   % Note : the result can be orthogonolized if necessary
end

%-----------------------------------------------------------------------

case 'ctestsp' % 
%=======================================================================
% [X1,X0] = spm_SpUtil('cTestSp',X,c)
% return the space tested and not tested by the contrast
% Note that unless X0 is reduced to a set of linearely independant 
% vectors, c will only be contained in the null space of X0. 
% If X0 is "reduced", then the "parent" space of c must be reduced
% as well for c to be the actual null space of X0.

if nargin<3, error('Insufficient arguments'), end
if ~spm_sp('isspc',varargin{3})
	c = spm_sp('set',varargin{3});
else c = varargin{3}; end


if spm_sp('isspc',varargin{2})
   if size(c.X,1) ~= size(varargin{2}.X,2)
	error('dimensions argin 2 and 3 dont match'); end

   rk = varargin{2}.rk;
   if rk == 0 error('Rank is null '), end

   % X_0 = varargin{2}.X*orth(spm_sp('res',c,varargin{2}.X'));
   X0 = varargin{2}.X*spm_sp('res',c,varargin{2}.v(:,1:rk));
   X1 = varargin{2}.X*c.X - X0*pinv(X0)*varargin{2}.X*c.X;


   varargout = {X1, X0};
					% can be orthogonalised if necessary
else
   if size(c.X,1) ~= size(varargin{2},2)
	error('dimensions argin 2 and 3 dont match'); end
   
   X0 = varargin{2}*spm_sp('res',c,orth(varargin{2}'));
   X1 = varargin{2}*spm_sp('oP',c,orth(varargin{2}'));
   X1 = varargin{2}*c.X - X0*pinv(X0)*varargin{2}*c.X;
   varargout = {X1,X0};
					% can be orthogonalised if necessary
end;

case 'itestsp' % 
%=======================================================================
% [X1] = spm_SpUtil('iTestSp',X,i0)
% return the space tested while keeping X(i0) 
% Note that unless X0 is reduced to a set of linearely independant 
% vectors, c will only be contained in the null space of X0. 
% If X0 is "reduced", then the "parent" space of c must be reduced
% as well for c to be the actual null space of X0.

if nargin<3, error('Insufficient arguments'), end
i0 = varargin{3};

if spm_sp('isspc',varargin{2})
   p = size(varargin{2}.X,2);
   if max(size(i0))==p & all( ismember(i0,[0 1]) ) i0 = find(i0); end;
   if any(floor(i0)~=i0) | any(i0<1) | any(i0>p)
	error('vector of column indices required'); end
   i1   = setdiff([1:p],i0);	%-Indices of columns hypothesised redundant
   if isempty(i1) varargout = {spm_sp('oP',varargin{2})}; 
   else 
      varargout = {varargin{2}.X(:,i1) - ...
	varargin{2}.X(:,i0)*pinv(varargin{2}.X(:,i0))*varargin{2}.X(:,i1)};
   end
else 
   p = size(varargin{2},2);
   if max(size(i0))==p & all( ismember(i0,[0 1]) ) i0 = find(i0); end;
   if any(floor(i0)~=i0) | any(i0<1) | any(i0>p)
	error('vector of column indices required'); end;
   i1   = setdiff([1:p],i0);	%-Indices of columns hypothesised redundant
   if isempty(i1) varargout = {varargin{2}*pinv(varargin{2})}; 
   else 
      varargout = {varargin{2}(:,i1) -  ...
	      varargin{2}(:,i0)*pinv(varargin{2}(:,i0))*varargin{2}(:,i1)};
   end
end


case 'trrv'
%=======================================================================
% [trRV  [trRVRV]]= spm_SpUtil('trRV', x [, V])
% returns trace(RV) or trace(R) if only one input argument.
% returns [trace(RV) trace(RVRV)] if two ouput arguments.
% This uses the Karl's cunning understanding of the trace:
% (tr(A*B) = sum(sum(A'*B)).
% If the space of X is set, then it uses x.u to avoid 
% extra computation. 

if nargin == 1, error('insufficient arguments'), end
if nargin == 2, 			%- no V
   
   if spm_sp('isspc',varargin{2}) % assumes it's filled 
     	varargout = {size(varargin{2}.X,1) - varargin{2}.rk, ...
			size(varargin{2}.X,1) - varargin{2}.rk};
   else
	rk  = rank(varargin{2});
	varargout = {size(varargin{2},1) - rk, ...
		     size(varargin{2},1) - rk}; 
   end;
   
else 					%- V provided, assumed correct !
   if spm_sp('isspc',varargin{2}) % assumes varargin{2} is set
	
	rk   = varargin{2}.rk; 
        if rk == 0 error('Rank is null '), end

	if nargout==1 		% only trRV is needed then use : 
	   trMV	= sum(sum( varargin{2}.u(:,1:rk)' .* ...
			   (varargin{2}.u(:,1:rk)'*varargin{3}) ));
	   varargout = { trace(varargin{3}) - trMV};
	else 		 		%  trMV is needed then use : 
 	   % assume trRVRV is needed as well : 

	   MV	= varargin{2}.u(:,1:rk)*(varargin{2}.u(:,1:rk)'*varargin{3});
	   trMV	= trace(MV);		%- possible to avoid computing MV;
					%- Use MV = u*(u'*V) 
					
	   trRVRV	= sum(sum(varargin{3}.*varargin{3})) - ...
			2*sum(sum(varargin{3}.*MV)) + ...
			sum(sum(MV'.*MV));

	   varargout = {(trace(varargin{3}) - trMV), trRVRV};
	end
   else 
	if nargout==1 		% only trRV is needed then use : 
	   trMV	= sum(sum( varargin{2}' .* (pinv(varargin{2})*varargin{3}) ));
	   varargout = { trace(varargin{3}) - trMV};
	else 		 		%  trMV is needed then use : 
 	   % assume trRVRV is needed as well : 

	   MV	= varargin{2}*(pinv(varargin{2})*varargin{3});
	   trMV	= trace(MV);		%- possible to avoid computing MV;
					%- Use MV = u*(u'*V) 
					
	   trRVRV	= sum(sum(varargin{3}.*varargin{3})) - ...
			2*sum(sum(varargin{3}.*MV)) + ...
			sum(sum(MV'.*MV));

	   varargout = {(trace(varargin{3}) - trMV), trRVRV};
	end
   end
end

case 'trmv'
%=======================================================================
% [trMV [trMVMV]] = spm_SpUtil('trMV', x [, V])
% returns trace(MV) or trace(M) if only one input argument.
% returns [trace(MV) trace(MVMV)] if two ouput arguments.
% This uses the Karl's cunning understanding of the trace:
% (tr(A*B) = sum(sum(A'.*B)).
% If the space of X is set, then it uses x.u to avoid 
% extra computation. 

if nargin == 1, error('insufficient arguments'), end
if nargin == 2, 			%- no V
   
   if spm_sp('isspc',varargin{2}) % assumes it's filled 
     	varargout = {varargin{2}.rk, varargin{2}.rk};
   else rk = rank(varargin{2}); varargout = {rk,rk}; end;
   
else 					%- V provided, and assumed correct !

   if spm_sp('isspc',varargin{2})	% if varargin{2} is set as a space

	rk   = varargin{2}.rk; 
        if rk == 0 error('Rank is null '), end

	if nargout==1 		% only trMV is needed then use : 
	   trMV   = sum(sum(varargin{2}.u(:,1:rk)' .* ...
 			    (varargin{2}.u(:,1:rk)'*varargin{3}) ));
	   varargout = {trMV};
	else 
 	   % assume trMVMV is needed as well : 

	   MV	= varargin{2}.u(:,1:rk)*(varargin{2}.u(:,1:rk)'*varargin{3});
				%- See note on saving mem. for MV in 'trRV'
	   trMVMV	=  sum(sum(MV'.*MV));
	   disp([trace(MV)  trMVMV]);
	   varargout = {trace(MV), trMVMV};
	end

   else 
	if nargout==1 		% only trMV is needed then use : 
	   trMV = sum(sum(varargin{2}'.* (pinv(varargin{2})*varargin{3}) ))
	   varargout = {trMV};
	else
	   MV	=  varargin{2}*(pinv(varargin{2})*varargin{3});
	   trMVMV	=  sum(sum(MV'.*MV));
	   varargout = {trace(MV), trMVMV};
	end
   end
end



case 'betarc'
%=======================================================================
% ('BetaRc', x, c)
% input : space x of the design and {F|t}contrast  
% output : R such that beta'*R*beta = extra sum of square for contrast c


if nargin < 3, error('insufficient arguments'), end
c = varargin{3};
if spm_sp('isspc',varargin{2})
   if ~spm_sp('isinspp',varargin{2},varargin{3}), ...
		error('contrast not estimable'); end;
   varargout = { c*pinv(c'*spm_sp('pinvXpX',varargin{2})*c)*c' };
				% Note that there is a better way for 
				% that: first sets the space of c ..
else % no space provided, use 
   varargout = { c*pinv(c'*pinv(varargin{2}'*varargin{2})*c)*c' };
end

case 'mpc'
%=======================================================================
% ('Mpc', x, c)
% input : space x of the design and {F|t}contrast  
% output : Mp such that Y'*Mp*Y = extra sum of square

if nargin == 1, error('insufficient arguments'), end
c	= varargin{3};
if spm_sp('isspc',varargin{2})
   if ~spm_sp('isinspp',varargin{2},varargin{3}), 
				error('contrast not estimable'); end;
   varargout = ...
   {varargin{2}.X*c*pinv(c'*spm_sp('xpx', varargin{2})*c)*c'* varargin{2}.X'};
				% Note that there is a better way for 
				% that with first setting the space of c.
else % no space provided, use 
   varargout = ...
   {varargin{2}*c*pinv(c'*varargin{2}'*varargin{2}*c)*c'*varargin{2}'};
end

case 'mpx1'
%=======================================================================
% ('MpX1',x,i1)
% input : space x of the design and subspace of x you want to test for.
% output : MpX1 such that Y'*MpX1*Y = extra sum of square

if nargin<3, error('insufficient arguments'), end
i1 = varargin{3};

%-Argument checks
%-----------------------------------------------------------------------


if spm_sp('isspc',varargin{2})
   p = size(varargin{2}.X,2);
   if max(size(i1))==p & all( ismember(i1,[0 1]) ) i1 = find(i1); end;
   if any(floor(i1)~=i1) | any(i1<1) | any(i1>p)
	error('vector of column indices required'); end;
   P = varargin{2}.X(:,i1); 
   varargout = {varargin{2}.X(:,i1)*...
      pinv(varargin{2}.X(:,i1)'*varargin{2}.X(:,i1))*varargin{2}.X(:,i1)'};
else 
   p = size(varargin{2},2);
   if max(size(i1))==p & all( ismember(i1,[0 1]) ) i1 = find(i1); end;
   if any(floor(i1)~=i1) | any(i1<1) | any(i1>p)
	error('vector of column indices required'); end;
   %- P = varargin{2}(:,i1); 
   varargout = {varargin{2}(:,i1)*...
	pinv(varargin{2}(:,i1)'*varargin{2}(:,i1))*varargin{2}(:,i1)'};
end


case 'mpx0'
%=======================================================================
% MpX1 = spm_SpUtil('MpX0',x,i0)
% input : space x and i0: the subspace of x you dont want keep.
% output : MpX1 such that Y'*MpX1*Y = extra sum of square

%-Argument checks
%-----------------------------------------------------------------------
if nargin<3, error('insufficient arguments'), end
i0 = varargin{3};


if spm_sp('isspc',varargin{2})
   p = size(varargin{2}.X,2);
   if max(size(i0))==p & all( ismember(i0,[0 1]) ) i0 = find(i0); end;
   if any(floor(i0)~=i0) | any(i0<1) | any(i0>p)
	error('vector of column indices required'); end;
   i1   = setdiff([1:p],i0);	%-Indices of columns hypothesised redundant
   if isempty(i1) varargout = {spm_sp('oP',varargin{2})};
   else 
      P    = varargin{2}.X(:,i1) - ...
	  varargin{2}.X(:,i0)*pinv(varargin{2}.X(:,i0))*varargin{2}.X(:,i1);
      varargout = {P*pinv(P)};
   end
else 
   p = size(varargin{2},2);
   if max(size(i0))==p & all( ismember(i0,[0 1]) ) i0 = find(i0); end;
   if any(floor(i0)~=i0) | any(i0<1) | any(i0>p)
	error('vector of column indices required'); end;
   i1   = setdiff([1:p],i0);	%-Indices of columns hypothesised redundant
   if isempty(i1) varargout = {varargin{2}*pinv(varargin{2})}; 
   else 
      P = varargin{2}(:,i1) - ...
	  varargin{2}(:,i0)*pinv(varargin{2}(:,i0))*varargin{2}(:,i1); 
      varargout = {P*pinv(P)};
   end
end


case 'cxpequi'
%=======================================================================
% c = spm_SpUtil('cxpequi',x,c)
% input : space x of the design and c a contrast.
% output : 0 if x.X' and c dont span the same space

if nargin<3, error('insufficient arguments'), end
c = varargin{3};

if spm_sp('isspc',varargin{2})
	varargout = {all(spm_sp('isinspp',varargin{2},c)) & ...
	 	  all(spm_sp('isinsp',spm_sp('set',c), varargin{2}.X'))}; 
else 
	rkc	= rank(c);
	if size(c,1) ~= size(varargin{2},2)
	   error('dimensions dont match '), end
	varargout = { 	(rank([varargin{2}' c]) == rkc) & ...
			(rkc == rank(varargin{2}')) };
end

case 'edf'
%=======================================================================
%[df1 df2] = spm_SpUtil('edf', x, i0, V)
% returns df1 and df2 the residual df for the projector onto the null 
% space of x' (residual forming projector) and the numerator of 
% the F-test where i0 are the columns for the null hypothesis model.
    
   if spm_sp('isspc',varargin{2}), p = size(varargin{2}.X,2);
   else p = size(varargin{2},2); end;

   i0 = varargin{3};
   if max(size(i0))==p & all( ismember(i0,[0 1]) ) i0 = find(i0); end;
   if any(floor(i0)~=i0) | any(i0<1) | any(i0>p)
	error('vector of column indices required'); end;

   [trRV trRVRV] = spm_SpUtil('trRV', varargin{2}, varargin{4});
   [trMpV trMpVMpV] =  ...
	spm_SpUtil('trMV',spm_SpUtil('iTestSp',varargin{2}, i0), ...
		    varargin{4});
   varargout = {trMpV^2/trMpVMpV, trRV^2/trRVRV};

otherwise
%=======================================================================
error('Unknown action string in spm_SpUtil')



%=======================================================================
end
