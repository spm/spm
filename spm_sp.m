function varargout = spm_sp(varargin)
% Seting an orthogonal space of a (design) matrix
% FORMAT varargout = spm_spc(action,varargin)
%
% This function computes the different projectors related
% to the space of the line and the rows of X. It should be
% used to avoid redundant computation on large X matrix. 
% It's divided into action that set up the space (fillAll- Create)
% and actions that compute projections (pinv, pinvXpX, pinvXXp, ...)
% This is motivated by the problem of rounding errors that 
% can invalidate some computation and is a tool to work with spaces.
%
% The only thing that is not easely computed is the null space of
% the line of X (assuming size(X,1) > size(X,2)).
% To get this space (a basis of it or a projector on it) use spm_sp on X'.
%
% The structure looks like : 
%	x = struct( ...	
%		'X',[], ...		% Mtx
%		'tol',[], ...		% tolerance
%		'ds',[], ...		% vectors of singular values 
%		'u',[],	...		% u as in X = u*diag(ds)*v'
%		'v',[], ...		% v as in X = u*diag(ds)*v'
%		'rk',[], ...		% rank
%		'oP', [], ...		% orthogonal projector on X
%		'oPp',[], ...		% orthogonal projector on X'
%		'ups',[], ...		% space in which this one is embeded
%		'sus',[] ...		% subspace 
%	);
%
% action can be :
% case 'set',
%-----------------------------------------------------------------------
%- use svd(X,0) to set up singular values, space basis 
%- rank and tolerance	
%
% case {'op', 'opp'},
%-----------------------------------------------------------------------
%-  ('oP', spc_str[, Data])
%-  ('oPp', spc_str[, Data])
%- 'oP' stands for orthogonal projection on the space of the columns
%- 'oPp' stands for orthogonal projection on the space of the lines
%- 'oP' and 'oPp' are applied to the Data when provided.
%
% case 'pinv'
%-----------------------------------------------------------------------
%- 'pinv' stands for the standart pinv 
%- Note that because size(pinv(X)) == size(X'), it is not generally
%- useful to compute pinv(X)*Data sequentially (as it is the case
%- for 'res' or 'oP'
%
% case 'pinvxpx'
%-----------------------------------------------------------------------
%- 'pinvxpx' stands the standart pinv(X'X) 
%- Note that because size(pinv(X'*X)) == [size(X,2) size(X,2)],
%- it is not useful to compute pinv(X'X)*Data sequentially unless 
%- size(X,1) < size(X,2)
%
% case 'xpx'
%-----------------------------------------------------------------------
%- 'xpx' stands the standart (X'X) 
%
% case 'pinvxxp'
%-----------------------------------------------------------------------
%- 'pinvxxp' is pinv(XX') 
%
% case 'xxp'
%-----------------------------------------------------------------------
%- 'xxp' is (XX') 
%
% case {'isinsp', 'isinspp'}	
%-----------------------------------------------------------------------
%- ('isinsp', spc_str, vector(s), [tolerance])   
%- check wether vectors are in the space of sp_tructure
%
% case {'null', 'opnull'} 
%-----------------------------------------------------------------------
%- ('null', spc_str)
%- ('opnull', spc_str [, Data])
%- null simply returns the null space  (same as matlab NULL)
%- opnull returns the orth. proj. onto the null space of spc_str.
%- ('opnull', spc_str, Data) returns Data projected onto the null space of
%- the spc_str.
%- note that if y = spm_sp('set',x.X'), we have spm_sp('opnull',x) == 
%- spm_sp('res',y) or, equivalently:  
%- spm_sp('opnull',x) + spm_sp('oP',y) == eye(size(y.X,1));
%
%  case 'res'
%-----------------------------------------------------------------------
% ('res', spc_str [, Data]) 
%- returns the residual forming matrix wrt the space spc_str, or this 
%- matrix times the Data when Data are given. This will be more
%- efficient than spm_sp('res', spc_str)*Data when size(X,1) > size(X,2)
%- This can also be seen as the orthogonal projector onto the null space
%- of x.X' (which is not generally computed in svd, unless
%- size(X,1) < size(X,2)).
%
% case {'ox', 'oxp'}
%-----------------------------------------------------------------------
%- ('ox', spc_str)
%- ('oxp', spc_str)
%- oX gives orth(x.X)
%- oXp gives ONE orthonormal basis for x.X' but not the same as orth(x.X')
%
%  case 'isspc'	
%-----------------------------------------------------------------------
%- ('isspc', spc_str)	is space structure ?
%
%  case 'isset'	
%-----------------------------------------------------------------------
%- ('isset', spc_str)	is space set ?
%
% The only restriction on the use of space structure is when X is so big
% that you can't allow to have X and its svd in mem at the same time. 
% Otherwise, the use of spm_sp will generally speed up computations
% and better deal with the memory. 
%
%_______________________________________________________________________
% %W% Jean-Baptiste Poline %E%

if nargin==0, error('do what? no arguments given...')
	else, action = varargin{1}; end


switch lower(action), 


   case 'create', 

	x = struct( ...	
		'X',[], ...		% Mtx
		'tol',[], ...		% tolerance
		'ds',[], ...		% vectors of singular values 
		'u',[],	...		% u as in X = u*diag(ds)*v'
		'v',[], ...		% v as in X = u*diag(ds)*v'
		'rk',[], ...		% rank
		'oP', [], ...		% orthogonal projector on X
		'oPp',[], ...		% orthogonal projector on X'
		'ups',[], ...		% space in which this one is embeded
		'sus',[] ...		% subspace 
	);
		% 'pinvX', [], ...	% pinv of X
	
	varargout = {x};

   case 'set',
   %-----------------------------------------------------------------------
   %- use svd(X,0) to set up singular values, space basis 
   %- rank and tolerance	

	if nargin==1 error('no design matrix : can''t do much'); end
	if isempty(varargin{2}), 
		x = spm_sp('create'); 
		varargout = {x};
		return, 
	end;
 
	x = spm_sp('create');

	x.X = varargin{2};

	%-- compute the svd : find all the singular values of such x.X.
    	%-- SVD(FULL(A)) will usually perform better than SVDS(A,MIN(SIZE(A)))

	[x.u s x.v] = svd(full(x.X),0);

	if size(x.X,1) == 1, x.ds = s; else x.ds = diag(s); end;
	clear s;

	%-- compute the tolerance
	x.tol =  max(size(x.X))*max(abs(x.ds))*eps;

	%-- compute the rank
	x.rk =  sum(x.ds > x.tol);

	varargout = {x};	

   case {'op', 'opp'},
   %-----------------------------------------------------------------------
   %-  ('oP', spc_str[, Data])   
   %-  ('oPp', spc_str[, Data])   
   %- 'oP' stands for orthogonal projection on the space of the columns
   %- 'oPp' stands for orthogonal projection on the space of the lines
   %- 'oP' and 'oPp' are applied to the Data when provided. 

	if nargin==1 error('no structure : can''t do much'); end
	%-- check
	if ~isfield(varargin{2},fieldnames(spm_sp('create')));
	  error('Second argument is not of the appropriate structure'); end 
	if ~isset_(varargin{2})
	  error('Fill the structure first or use pinv'); end 
	if nargin>4 error('too many input arg'); end

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

   case 'pinv'
   %-----------------------------------------------------------------------
   %- 'pinv' stands for the standart pinv 
   %- Note that because size(pinv(X)) == size(X'), it is not generally
   %- useful to compute pinv(X)*Data sequentially (as it is the case
   %- for 'res' or 'oP'

	if nargin==1 error('no structure : can''t do much'); end
	%-- check
	if ~isfield(varargin{2},fieldnames(spm_sp('create')));
	  error('Second argument is not of the appropriate structure'); end 
	if ~isset_(varargin{2})
	  error('Fill the structure first or use pinv'); end 

	%-- compute the pseudo-inverse
	p   = max(size(varargin{2}.ds));
	ds1 = zeros(size(varargin{2}.ds)); 
	ds1([1:varargin{2}.rk]) = varargin{2}.ds([1:varargin{2}.rk]).^(-1);
	if p == size(varargin{2}.v,1)
	   varargout = { varargin{2}.v*diag(ds1)*varargin{2}.u' };
	else %- then p < size(varargin{2}.v,1)
	   varargout = { varargin{2}.v(:,1:p)*diag(ds1)*varargin{2}.u' };
	end

   case 'pinvxpx'
   %-----------------------------------------------------------------------
   %- 'pinvxpx' stands the standart pinv(X'X) 
   %- Note that because size(pinv(X'*X)) == [size(X,2) size(X,2)],
   %- it is not useful to compute pinv(X'X)*Data sequentially unless 
   %- size(X,1) < size(X,2)
	
	if nargin==1 error('no structure : can''t do much'); end
	%-- check
	if ~isfield(varargin{2},fieldnames(spm_sp('create')));
	  error('Second argument is not of the appropriate structure'); end 
	if  ~isset_(varargin{2})
	  error('Fill the structure first or use pinv'); end 

	%-- compute the pseudo-inverse of X'X
	p   = max(size(varargin{2}.ds));
	ds1 = zeros(size(varargin{2}.ds)); 
	ds1([1:varargin{2}.rk]) = varargin{2}.ds([1:varargin{2}.rk]).^(-2);
	if p == size(varargin{2}.v,1)
	  varargout = { varargin{2}.v*diag(ds1)*varargin{2}.v' };
	else %- then p < size(varargin{2}.v,1)
	  varargout = { varargin{2}.v(:,1:p)*diag(ds1)*varargin{2}.v(:,1:p)' };
	end

   case 'xpx'
   %-----------------------------------------------------------------------
   %- 'xpx' stands the standart (X'X) 
	
	if nargin==1 error('no structure : can''t do much'); end
	%-- check
	if ~isfield(varargin{2},fieldnames(spm_sp('create')));
	  error('Second argument is not of the appropriate structure'); end 
	if  ~isset_(varargin{2})
	  error('Fill the structure first or use pinv'); end 

	%-- compute  X'X
	p   = max(size(varargin{2}.ds));
	ds1 = zeros(size(varargin{2}.ds)); 
	ds1([1:varargin{2}.rk]) = varargin{2}.ds([1:varargin{2}.rk]).^2;
	if p == size(varargin{2}.v,1)
	  varargout = { varargin{2}.v*diag(ds1)*varargin{2}.v' };
	else %- then p < size(varargin{2}.v,1)
	  varargout = { varargin{2}.v(:,1:p)*diag(ds1)*varargin{2}.v(:,1:p)' };
	end

   case 'pinvxxp'
   %-----------------------------------------------------------------------
   %- 'pinvxxp' is pinv(XX') 
	
	if nargin==1 error('no structure : can''t do much'); end
	%-- check
	if ~isfield(varargin{2},fieldnames(spm_sp('create')));
	  error('Second argument is not of the appropriate structure'); end 
	if  ~isset_(varargin{2})
	  error('Fill the structure first or use pinv'); end 

	%-- compute the pseudo-inverse of X'X
	ds1 = zeros(size(varargin{2}.ds)); 
	ds1([1:varargin{2}.rk]) = varargin{2}.ds([1:varargin{2}.rk]).^(-2);
	varargout = { varargin{2}.u*diag(ds1)*varargin{2}.u' }; clear ds1;
	
   case 'xxp'
   %-----------------------------------------------------------------------
   %- 'xxp' is (XX') 
	
	if nargin==1 error('no structure : can''t do much'); end
	%-- check
	if ~isfield(varargin{2},fieldnames(spm_sp('create')));
	  error('Second argument is not of the appropriate structure'); end 
	if  ~isset_(varargin{2})
	  error('Fill the structure first or use pinv'); end 

	%-- compute the pseudo-inverse of X'X
	ds1 = zeros(size(varargin{2}.ds)); 
	ds1([1:varargin{2}.rk]) = varargin{2}.ds([1:varargin{2}.rk]).^2;
	varargout = { varargin{2}.u*diag(ds1)*varargin{2}.u' }; clear ds1;
	

   case {'isinsp', 'isinspp'}	
   %-----------------------------------------------------------------------
   %- ('isinsp', spc_str, vector(s), [tolerance])   
   %- check wether vectors are in the space of sp_tructure
	
	if nargin==1 error('no structure : can''t do much'); end
	%-- check structure
	if ~isfield(varargin{2},fieldnames(spm_sp('create')));
	  error('Second argument is not of the appropriate structure'); end 
	if  ~isset_(varargin{2})
	  error('Fill the structure first or use pinv'); end 
	if nargin==2 error('What is in the space of ... ?'); end
	if nargin > 3, tol = varargin{4}; else, tol =  varargin{2}.tol; end;

	switch lower(action)
	case 'isinsp'
	   %- check dimensions
	   if size(varargin{2}.X,1) ~= size(varargin{3},1) 
		error('Dimensions dont match...'); end
	   if ~isempty(varargin{2}.oP)
	      varargout = {all((varargin{2}.oP*varargin{3}-varargin{3})<tol)};
	   else
	      varargout = ...
		{all( (spm_sp('oP',varargin{2},varargin{3})-varargin{3})<tol )};
	   end
	case 'isinspp'
	   %- check dimensions
	   if size(varargin{2}.X,2) ~= size(varargin{3},1) 
		error('Dimensions dont match...'); end
	   if ~isempty(varargin{2}.oPp)
	      varargout = {all((varargin{2}.oPp*varargin{3}-varargin{3})<tol)};
	   else
	      varargout = ...
		{all( (spm_sp('oPp',varargin{2},varargin{3})-varargin{3})<tol)};
	   end

	end % 	switch lower(action)


   case {'null', 'opnull'} 
   %-----------------------------------------------------------------------
   %- ('null', spc_str)
   %- ('opnull', spc_str [, Data])
   %- null simply returns the null space  (same as matlab NULL)
   %- opnull returns the orth. proj. onto the null space of spc_str.
   %- ('opnull', spc_str, Data) returns Data projected onto the null space of
   %- the spc_str.
   %- note that if y = spm_sp('set',x.X'), we have spm_sp('opnull',x) == 
   %- spm_sp('res',y) or, equivalently:  
   %- spm_sp('opnull',x) + spm_sp('oP',y) == eye(size(y.X,1));

	if nargin==1 error('no structure : can''t do much'); end
	%-- check structure
	if ~isfield(varargin{2},fieldnames(spm_sp('create')));
	  error('Second argument is not of the appropriate structure'); end 
	if  ~isset_(varargin{2})
	  error('Fill the structure first or use pinv'); end 
	if nargin>3 error('too many input arg'); end

	rk = varargin{2}.rk;
	p  = size(varargin{2}.X,2);
	if p==rk | rk == 0, varargout = {[]}; return; end;

	switch lower(action)
	case 'null'
	    % null space = vectors associated with 0 eigenvalues
	    if nargin==3 error('too many input arg'); end
	    varargout = {varargin{2}.v(:,rk+1:p)};

	case 'opnull'
	     if nargin==3	% Apply to arg proj. on the 0-space of x.X 
				% check consistency of third argument
	        if size(varargin{3},1) ~= size(varargin{2}.X,2)
	  	   error('Inconsistant arg size(in{3},1)~size(in{2}.X,2)'); end;
	        varargout = {varargin{2}.v(:,rk+1:p)*varargin{2}.v(:,rk+1:p)'...
				*varargin{3} }; 
	     else
	   	varargout = {varargin{2}.v(:,rk+1:p)*varargin{2}.v(:,rk+1:p)'};
	   	% or varargout = {eye(size(varargin{2}.oPp))-varargin{2}.oPp};
	     end
	end % switch

   case 'res'
   %-----------------------------------------------------------------------
   % ('res', spc_str [, Data]) 
   %- returns the residual forming matrix wrt the space spc_str, or this 
   %- matrix times the Data when Data are given. This will be more
   %- efficient than spm_sp('res', spc_str)*Data when size(X,1) > size(X,2)
   %- This can also be seen as the orthogonal projector onto the null space
   %- of x.X' (which is not generally computed in svd, unless
   %- size(X,1) < size(X,2)).

   	if nargin==1 error('no structure : can''t do much'); end
	%-- check structure
	if ~isfield(varargin{2},fieldnames(spm_sp('create')));
	  error('Second argument is not of the appropriate structure'); end 
	if  ~isset_(varargin{2})
	  error('Fill the structure first or use pinv'); end 
	if nargin>4 error('too many input arg'); end

	rk = varargin{2}.rk;
	if rk > 0 
	  switch nargin
	  case 2
	     if isempty(varargin{2}.oP)
	        %---------- avoid storing varargin{2}.oP = ...
		%---------- spm_sp('oP',varargin{2});
	        varargout = { eye(size(varargin{2}.X,1)) - ...
			varargin{2}.u(:,[1:rk])*varargin{2}.u(:,[1:rk])' };
	     else
	        varargout = {eye(size(varargin{2}.oP))-varargin{2}.oP};
	     end

	  case 3
	     if size(varargin{3},1) ~= size(varargin{2}.X,1) 
		error('Data and space dim. dont match '); 
	     else 
		varargout = { 	varargin{3} - varargin{2}.u(:,[1:rk])* ...
				( varargin{2}.u(:,[1:rk])'*varargin{3} ) };
	     end

	  end %--- switch nargin

	else %--- if rk > 0
	  varargout = { 0 };
	end %---- if rk > 0


   case {'ox', 'oxp'}
   %-----------------------------------------------------------------------
   %- ('ox', spc_str)
   %- ('oxp', spc_str)
   %- oX gives orth(x.X)
   %- oXp gives ONE orthonormal basis for x.X' but not the same as orth(x.X')
	
	if nargin==1 error('no structure : can''t do much'); end
	%-- check structure
	if ~isfield(varargin{2},fieldnames(spm_sp('create')));
	  error('Second argument is not of the appropriate structure'); end 
	if  ~isset_(varargin{2})
	  error('Fill the structure first or use pinv'); end 

	rk = varargin{2}.rk;
	if rk > 0 
	   switch lower(action)
	      case 'ox'
	   	varargout = {varargin{2}.u(:,[1:rk])};
	      case 'oxp'
	   	varargout = {varargin{2}.v(:,[1:rk])};
	   end
	else varargout = {0};
	end


   case 'isspc'	
   %-----------------------------------------------------------------------
   %- ('isspc', spc_str)	is space structure ?
   if nargin~=2, error('Wrong number of argument'); end
   varargout = {isfield(varargin{2},fieldnames(spm_sp('create')))};

   case 'isset'	
   %-----------------------------------------------------------------------
   %- ('isset', spc_str)	is space set ?
   if nargin~=2, error('Wrong number of argument'); end
   varargout = {isset_(varargin{2})};

   otherwise
   %-----------------------------------------------------------------------
	error('Unknown action string')

end % switch


%--------------------------------------------
function b = isset_(x)
b = (	isempty(x.X) 	|...
	isempty(x.u) 	|...
	isempty(x.v) 	|...
	isempty(x.ds)	|...
	isempty(x.tol)	|...
	isempty(x.rk)	...
	 );
b = ~b;
return;
%--------------------------------------------



