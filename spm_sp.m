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
%   x = struct(...
%       'X',    [],...      % Mtx
%       'tol',  [],...      % tolerance
%       'ds',   [],...      % vectors of singular values 
%       'u',    [],...      % u as in X = u*diag(ds)*v'
%       'v',    [],...      % v as in X = u*diag(ds)*v'
%       'rk',   [],...      % rank
%       'oP',   [],...      % orthogonal projector on X
%       'oPp',  [],...      % orthogonal projector on X'
%       'ups',  [],...      % space in which this one is embeded
%       'sus',  []);        % subspace
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
% FORMAT b = spm_sp('eachinsp',x,c[,tol])
% FORMAT b = spm_sp('eachinspp',x,c[,tol])
% Same as 'isinsp' and 'isinspp' but returns a logical row vector of
% length size(c,2).
%
% FORMAT N = spm_sp('n',x)
% Simply returns the null space of matrix X (same as matlab NULL)
% (Null space = vectors associated with zero eigenvalues)
% x - space structure of matrix X
% N - null space
%
%
% FORMAT r = spm_sp('nop',x[,Y])
% Orthogonal projector onto null space of X, or projection of data Y (if passed)
% x - space structure of matrix X
% Y - (optional) data
% r - (if no Y passed) orthogonal projection matrix  into the null space of X
%   - (if Y passed   ) orthogonal projection of data into the null space of X
% ( Note that if xp = spm_sp('set',x.X'), we have:                    )
% (       spm_sp('nop',x) == spm_sp('res',xp)                      )
% ( or, equivalently:                                                 )
% (       spm_sp('nop',x) + spm_sp('oP',xp) == eye(size(xp.X,1));  )
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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jean-Baptiste Poline
% $Id: spm_sp.m 1143 2008-02-07 19:33:33Z spm $


if nargin==0
    error('Do what? no arguments given...')
else
    action = varargin{1};
end

%- check the very basics arguments

switch lower(action), 
case {'create','set','issetspc','isspc'}
    %- do nothing
otherwise,
    if nargin==1, error('No space : can''t do much!'), end
    [ok,str] = spm_sp('issetspc',varargin{2}); 
    if ~ok, error(str), else, sX = varargin{2}; end;
end;



switch lower(action), 

case 'create'             %-Create space structure
%=======================================================================
% x = spm_sp('Create')
varargout = {sf_create};

case 'set'          %-Set singular values, space basis, rank & tolerance
%=======================================================================
% x = spm_sp('Set',X)

if nargin==1 error('No design matrix : can''t do much!'), 
else X = varargin{2}; end
if isempty(X), varargout = {sf_create}; return, end

%- only sets plain matrices
%- check X has 2 dim only

if max(size(size(X))) > 2, error('Too many dim in the set'), end 
if ~isnumeric(X), error('only sets numeric matrices'), end

varargout = {sf_set(X)};


case {'p', 'transp'}   %-Transpose space of X
%=======================================================================
switch nargin
case 2  
    varargout = {sf_transp(sX)};
otherwise
    error('too many input argument in spm_sp');

end % switch nargin


case {'op', 'op:'}   %-Orthogonal projectors on space of X
%=======================================================================
% r = spm_sp('oP', sX[,Y])   
% r = spm_sp('oP:', sX[,Y])   %- set to 0 less than tolerence values
%
% if isempty(Y) returns as if Y not given
%-----------------------------------------------------------------------

switch nargin
case 2  
  switch lower(action),         
  case 'op'
     varargout = {sf_op(sX)};
  case 'op:'
     varargout = {sf_tol(sf_op(sX),sX.tol)};
  end %- switch lower(action), 

case 3
   Y = varargin{3};
   if isempty(Y), varargout = {spm_sp(action,sX)}; return, end
   if size(Y,1) ~= sf_s1(sX), error('Dim dont match'); end;

   switch lower(action), 
   case 'op'
     varargout = {sf_op(sX)*Y};
   case 'op:'
     varargout = {sf_tol(sf_op(sX)*Y,sX.tol)};
   end % switch lower(action)

otherwise
   error('too many input argument in spm_sp');

end % switch nargin

case {'opp', 'opp:'}   %-Orthogonal projectors on space of X'
%=======================================================================
% r = spm_sp('oPp',sX[,Y])   
% r = spm_sp('oPp:',sX[,Y])   %- set to 0 less than tolerence values
%
% if isempty(Y) returns as if Y not given
%-----------------------------------------------------------------------

switch nargin
case 2  
   switch lower(action),        
    case 'opp'
       varargout = {sf_opp(sX)};
    case 'opp:'
       varargout = {sf_tol(sf_opp(sX),sX.tol)};
   end %- switch lower(action), 

case 3
   Y = varargin{3};
   if isempty(Y), varargout = {spm_sp(action,sX)}; return, end
   if size(Y,1) ~= sf_s2(sX), error('Dim dont match'); end;

   switch lower(action), 
   case 'opp'
      varargout = {sf_opp(sX)*Y};
   case 'opp:'
      varargout = {sf_tol(sf_opp(sX)*Y,sX.tol)};
   end % switch lower(action)

otherwise
   error('too many input argument in spm_sp');

end % switch nargin


case {'x-','x-:'}                        %-Pseudo-inverse of X - pinv(X)
%=======================================================================
%  = spm_sp('x-',x) 

switch nargin
case 2  
   switch lower(action),        
      case {'x-'}
         varargout = { sf_pinv(sX) };
      case {'x-:'}
         varargout = {sf_tol( sf_pinv(sX), sf_t(sX) )};
   end

case 3
   %- check dimensions of Y 
   Y = varargin{3};
   if isempty(Y), varargout = {spm_sp(action,sX)}; return, end
   if size(Y,1) ~= sf_s1(sX), error(['Dim dont match ' action]); end

   switch lower(action),        
      case {'x-'}
         varargout = { sf_pinv(sX)*Y };
      case {'x-:'}
         varargout = {sf_tol( sf_pinv(sX)*Y, sf_t(sX) )};
   end

otherwise
   error(['too many input argument in spm_sp ' action]);

end % switch nargin


case {'xp-','xp-:','x-p','x-p:'}                 %- Pseudo-inverse of X'
%=======================================================================
% pX = spm_sp('xp-',x) 

switch nargin
case 2  
   switch lower(action),        
      case {'xp-','x-p'}
         varargout = { sf_pinvxp(sX) };
      case {'xp-:','x-p:'}
         varargout = {sf_tol( sf_pinvxp(sX), sf_t(sX) )};
   end

case 3
   %- check dimensions of Y 
   Y = varargin{3};
   if isempty(Y), varargout = {spm_sp(action,sX)}; return, end
   if size(Y,1) ~= sf_s2(sX), error(['Dim dont match ' action]); end

   switch lower(action),        
      case {'xp-','x-p'}
         varargout = { sf_pinvxp(sX)*Y };
      case {'xp-:','x-p:'}
         varargout = {sf_tol( sf_pinvxp(sX)*Y, sf_t(sX) )};
   end

otherwise
   error(['too many input argument in spm_sp ' action]);

end % switch nargin


case {'cukxp-','cukxp-:'}   %- Coordinates of pinv(X') in the base of uk
%=======================================================================
% pX = spm_sp('cukxp-',x) 

switch nargin
case 2  
   switch lower(action),        
      case {'cukxp-'}
         varargout = { sf_cukpinvxp(sX) };
      case {'cukxp-:'}
         varargout = {sf_tol(sf_cukpinvxp(sX),sX.tol)};
   end

case 3
   %- check dimensions of Y 
   Y = varargin{3};
   if isempty(Y), varargout = {spm_sp(action,sX)}; return, end
   if size(Y,1) ~= sf_s2(sX), error(['Dim dont match ' action]); end

   switch lower(action),        
      case {'cukxp-'}
         varargout = { sf_cukpinvxp(sX)*Y };
      case {'cukxp-:'}
         varargout = {sf_tol(sf_cukpinvxp(sX)*Y,sX.tol)};
   end

otherwise
   error(['too many input argument in spm_sp ' action]);

end % switch nargin




case {'cukx','cukx:'}                %- Coordinates of X in the base of uk 
%=======================================================================
% pX = spm_sp('cukx',x)

switch nargin
case 2  
   switch lower(action),        
      case {'cukx'}
         varargout = { sf_cukx(sX) };
      case {'cukx:'}
         varargout = {sf_tol(sf_cukx(sX),sX.tol)};
   end

case 3
   %- check dimensions of Y 
   Y = varargin{3};
   if isempty(Y), varargout = {spm_sp(action,sX)}; return, end
   if size(Y,1) ~= sf_s2(sX), error(['Dim dont match ' action]); end

   switch lower(action),        
      case {'cukx'}
         varargout = { sf_cukx(sX)*Y };
      case {'cukx:'}
         varargout = {sf_tol(sf_cukx(sX)*Y,sX.tol)};
   end

otherwise
   error(['too many input argument in spm_sp ' action]);

end % switch nargin


case {'rk'}                                            %- Returns rank
%=======================================================================
varargout = { sf_rk(sX) };


case {'ox', 'oxp'}                    %-Orthonormal basis sets for X / X'
%=======================================================================
% oX  = spm_sp('ox', x)
% oXp = spm_sp('oxp',x)

if sf_rk(sX) > 0 
   switch lower(action)
   case 'ox'
      varargout = {sf_uk(sX)};
   case 'oxp'
      varargout = {sf_vk(sX)};
   end
else
   switch lower(action)
   case 'ox'
      varargout = {zeros(sf_s1(sX),1)};
   case 'oxp'
      varargout = {zeros(sf_s2(sX),1)};
   end
end

case {'x', 'xp'}            %- X / X'           robust to spm_sp changes
%=======================================================================
% X  = spm_sp('x', x)
% X' = spm_sp('xp',x)

switch lower(action)
   case 'x', varargout = {sX.X};
   case 'xp', varargout = {sX.X'};
end

case {'xi', 'xpi'}   %- X(:,i) / X'(:,i)        robust to spm_sp changes
%=======================================================================
% X  = spm_sp('xi', x)
% X' = spm_sp('xpi',x)

i = varargin{3}; % NO CHECKING on i !!! assumes correct
switch lower(action)
   case 'xi',  varargout = {sX.X(:,i)};
   case 'xpi', varargout = {sX.X(i,:)'};
end

case {'uk','uk:'}                                    %- Returns u(:,1:r) 
%=======================================================================
% pX = spm_sp('uk',x) 
% Notice the difference with 'ox' : 'ox' always returns a basis of the
% proper siwe while this returns empty if rank is null

warning('can''t you use ox ?');

switch nargin

case 2  
   switch lower(action),        
      case {'uk'}
         varargout = { sf_uk(sX) };
      case {'uk:'}
         varargout = { sf_tol(sf_uk(sX),sX.tol) };
   end

case 3
   %- check dimensions of Y 
   Y = varargin{3};
   if isempty(Y), varargout = {spm_sp(action,sX)}; return, end
   if size(Y,1) ~= sf_rk(sX), error(['Dim dont match ' action]); end

   switch lower(action),        
      case {'uk'}
         varargout = { sf_uk(sX)*Y };
      case {'uk:'}
         varargout = {sf_tol(sf_uk(sX)*Y,sX.tol)};
   end

otherwise
   error(['too many input argument in spm_sp ' action]);

end % switch nargin




case {'pinvxpx', 'xpx-', 'pinvxpx:', 'xpx-:',}    %- Pseudo-inv of (X'X)
%=======================================================================
% pXpX = spm_sp('pinvxpx',x [,Y]) 

switch nargin
case 2  
   switch lower(action),        
    case {'xpx-','pinvxpx'}
        varargout = {sf_pinvxpx(sX)};
    case {'xpx-:','pinvxpx:'}
      varargout = {sf_tol(sf_pinvxpx(sX),sX.tol)};
   end %- 

case 3
   %- check dimensions of Y 
   Y = varargin{3};
   if isempty(Y), varargout = {spm_sp(action,sX)}; return, end
   if size(Y,1) ~= sf_s2(sX), error('Dim dont match'); end;

   switch lower(action),        
    case {'xpx-','pinvxpx'}
       varargout = {sf_pinvxpx(sX)*Y};
    case {'xpx-:','pinvxpx:'}
           varargout = {sf_tol(sf_pinvxpx(sX)*Y,sX.tol)};
   end %-

otherwise
   error('too many input argument in spm_sp');

end % switch nargin


case {'xpx','xpx:'}              %-Computation of (X'*X)
%=======================================================================
% XpX = spm_sp('xpx',x [,Y])

switch nargin
case 2  
   switch lower(action),        
    case {'xpx'}
        varargout = {sf_xpx(sX)};
    case {'xpx:'}
      varargout = {sf_tol(sf_xpx(sX),sX.tol)};
   end %- 

case 3
    %- check dimensions of Y    
    Y = varargin{3};
    if isempty(Y), varargout = {spm_sp(action,sX)}; return, end
    if size(Y,1) ~= sf_s2(sX), error('Dim dont match'); end;

   switch lower(action),        
    case {'xpx'}
        varargout = {sf_xpx(sX)*Y};
    case {'xpx:'}
      varargout = {sf_tol(sf_xpx(sX)*Y,sX.tol)};
   end %-

otherwise
   error('too many input argument in spm_sp');

end % switch nargin




case {'cx->cu','cx->cu:'}    %-coordinates in the basis of X -> basis u
%=======================================================================
% 
% returns cu such that sX.X*cx == sX.u*cu 

switch nargin
case 2  
   switch lower(action),        
    case {'cx->cu'}
        varargout = {sf_cxtwdcu(sX)};
    case {'cx->cu:'}
      varargout = {sf_tol(sf_cxtwdcu(sX),sX.tol)};
   end %- 

case 3
    %- check dimensions of Y    
    Y = varargin{3};
    if isempty(Y), varargout = {spm_sp(action,sX)}; return, end
    if size(Y,1) ~= sf_s2(sX), error('Dim dont match'); end;

   switch lower(action),        
    case {'cx->cu'}
        varargout = {sf_cxtwdcu(sX)*Y};
    case {'cx->cu:'}
      varargout = {sf_tol(sf_cxtwdcu(sX)*Y,sX.tol)};
   end %-

otherwise
   error('too many input argument in spm_sp');

end % switch nargin




case {'xxp-','xxp-:','pinvxxp','pinvxxp:'}    %-Pseudo-inverse of (XX')
%=======================================================================
% pXXp = spm_sp('pinvxxp',x [,Y])

switch nargin
case 2  
   switch lower(action),        
    case {'xxp-','pinvxxp'}
        varargout = {sf_pinvxxp(sX)};
    case {'xxp-:','pinvxxp:'}
      varargout = {sf_tol(sf_pinvxxp(sX),sX.tol)};
   end %- 

case 3
    %- check dimensions of Y    
    Y = varargin{3};
    if isempty(Y), varargout = {spm_sp(action,sX)}; return, end
    if size(Y,1) ~= sf_s1(sX), error('Dim dont match'); end;

   switch lower(action),        
    case {'xxp-','pinvxxp'}
        varargout = {sf_pinvxxp(sX)*Y};
    case {'xxp-:','pinvxxp:'}
      varargout = {sf_tol(sf_pinvxxp(sX)*Y,sX.tol)};
   end %-

otherwise
   error('too many input argument in spm_sp');

end % switch nargin
    

case {'xxp','xxp:'}                         %-Computation of (X*X')
%=======================================================================
% XXp = spm_sp('xxp',x)


switch nargin
case 2  
   switch lower(action),        
    case {'xxp'}
        varargout = {sf_xxp(sX)};
    case {'xxp:'}
      varargout = {sf_tol(sf_xxp(sX),sX.tol)};
   end %- 

case 3
    %- check dimensions of Y    
    Y = varargin{3};
    if isempty(Y), varargout = {spm_sp(action,sX)}; return, end
    if size(Y,1) ~= sf_s1(sX), error('Dim dont match'); end;

   switch lower(action),        
    case {'xxp'}
        varargout = {sf_xxpY(sX,Y)};
    case {'xxp:'}
      varargout = {sf_tol(sf_xxpY(sX,Y),sX.tol)};
   end %-

otherwise
   error('too many input argument in spm_sp');

end % switch nargin


case {'^p','^p:'}                %-Computation of v*(diag(s.^n))*v'
%=======================================================================

switch nargin
case {2,3}
    if nargin==2, n = 1; else n = varargin{3}; end;
    if ~isnumeric(n), error('~isnumeric(n)'), end;

   switch lower(action),        
    case {'^p'}
        varargout = {sf_jbp(sX,n)};
    case {'^p:'}
      varargout = {sf_tol(sf_jbp(sX,n),sX.tol)};
   end %- 

case 4
    n = varargin{3};
    if ~isnumeric(n), error('~isnumeric(n)'), end;
    Y = varargin{4};
    if isempty(Y), varargout = {spm_sp(action,sX,n)}; return, end
    if size(Y,1) ~= sf_s2(sX), error('Dim dont match'); end;

   switch lower(action),        
    case {'^p'}
        varargout = {sf_jbp(sX,n)*Y};
    case {'^p:'}
      varargout = {sf_tol(sf_jbp(sX,n)*Y,sX.tol)};
   end %- 

otherwise
   error('too many input argument in spm_sp');

end % switch nargin



case {'^','^:'}                %-Computation of v*(diag(s.^n))*v'
%=======================================================================

switch nargin
case {2,3}
    if nargin==2, n = 1; else n = varargin{3}; end;
    if ~isnumeric(n), error('~isnumeric(n)'), end;

   switch lower(action),        
    case {'^'}
        varargout = {sf_jb(sX,n)};
    case {'^:'}
      varargout = {sf_tol(sf_jb(sX,n),sX.tol)};
   end %- 

case 4
    n = varargin{3};
    if ~isnumeric(n), error('~isnumeric(n)'), end;
    Y = varargin{4};
    if isempty(Y), varargout = {spm_sp(action,sX,n)}; return, end
    if size(Y,1) ~= sf_s1(sX), error('Dim dont match'); end;

   switch lower(action),        
    case {'^'}
        varargout = {sf_jbY(sX,n,Y)};
    case {'^:'}
      varargout = {sf_tol(sf_jbY(sX,n,Y),sX.tol)};
   end %- 

otherwise
   error('too many input argument in spm_sp');

end % switch nargin



case {'n'}    %-Null space of sX
%=======================================================================

switch nargin
case 2  
    varargout = {sf_n(sX)};
otherwise
   error('too many input argument in spm_sp');
end % switch nargin


case {'np'}    %-Null space of sX'
%=======================================================================

switch nargin
case 2  
    varargout = {sf_n(sf_transp(sX))};
otherwise
   error('too many input argument in spm_sp');
end % switch nargin


case {'nop', 'nop:'}    %- project(or)(ion) into null space
%=======================================================================
%
%
% 

switch nargin

case 2  
   switch lower(action),        
    case {'nop'}
        n = sf_n(sX);
        varargout = {n*n'};
    case {'nop:'}
        n = sf_n(sX);
        varargout = {sf_tol(n*n',sX.tol)};
   end %-

case 3
    %- check dimensions of Y    
    Y = varargin{3};
    if isempty(Y), varargout = {spm_sp(action,sX)}; return, end
    if size(Y,1) ~= sf_s2(sX), error('Dim dont match'); end;

   switch lower(action),        
    case {'nop'}
        n = sf_n(sX);
        varargout = {n*(n'*Y)};
    case {'nop:'}
        n = sf_n(sX);
        varargout = {sf_tol(n*(n'*Y),sX.tol)};
   end %-

otherwise
   error('too many input argument in spm_sp');

end % switch nargin



case {'nopp', 'nopp:'}    %- projector(ion) into null space of X'
%=======================================================================
%
%

switch nargin

case 2  
    switch lower(action),       
    case {'nopp'}
        varargout = {spm_sp('nop',sf_transp(sX))};
    case {'nopp:'}
        varargout = {spm_sp('nop:',sf_transp(sX))};
    end %-
case 3  
    switch lower(action),       
    case {'nopp'}
        varargout = {spm_sp('nop',sf_transp(sX),varargin{3})};
    case {'nopp:'}
        varargout = {spm_sp('nop:',sf_transp(sX),varargin{3})};
    end %-
otherwise
   error('too many input argument in spm_sp');

end % switch nargin


case {'res', 'r','r:'}             %-Residual formaing matrix / residuals
%=======================================================================
% r = spm_sp('res',sX[,Y])
%
%- 'res' will become obsolete : use 'r' or 'r:' instead
%- At some stage, should be merged with 'nop'


switch nargin

case 2  
    switch lower(action) 
    case {'r','res'} 
        varargout = {sf_r(sX)};
    case {'r:','res:'} 
        varargout = {sf_tol(sf_r(sX),sX.tol)};
    end %-
case 3

    %- check dimensions of Y    
    Y = varargin{3};
    if isempty(Y), varargout = {spm_sp(action,sX)}; return, end
    if size(Y,1) ~= sf_s1(sX), error('Dim dont match'); end;
    
    switch lower(action) 
    case {'r','res'} 
        varargout = {sf_rY(sX,Y)};
    case {'r:','res:'} 
        varargout = {sf_tol(sf_rY(sX,Y),sX.tol)};
    end %-

otherwise
   error('too many input argument in spm_sp');

end % switch nargin


case {':'}
%=======================================================================
% spm_sp(':',sX [,Y [,tol]])

%- Sets Y and tol according to arguments

if nargin > 4
    error('too many input argument in spm_sp'); 

else
    if nargin > 3
        if isnumeric(varargin{4}), tol = varargin{4}; 
        else error('tol must be numeric'); 
        end
    else
        tol = sX.tol;
    end
    if nargin > 2
        Y = varargin{3}; %- if isempty, returns empty, end
    else
        Y = sX.X;   
    end
end

varargout = {sf_tol(Y,tol)};




case {'isinsp', 'isinspp'}    %- is in space or is in dual space
%=======================================================================
% b = spm_sp('isinsp',x,c[,tol])
% b = spm_sp('isinspp',x,c[,tol])
%-Check whether vectors are in row/column space of X

%-Check arguments
%-----------------------------------------------------------------------
if nargin<3, error('insufficient arguments - action,x,c required'), end
c = varargin{3}; %- if isempty(c), dim wont match exept for empty sp.
if nargin<4, tol=sX.tol; else, tol = varargin{4}; end

%-Compute according to case
%-----------------------------------------------------------------------
switch lower(action)

case 'isinsp'
    %-Check dimensions
    if size(sX.X,1) ~= size(c,1) 
        warning('Vector dim don''t match col. dim : not in space !'); 
        varargout = { 0 }; return;
    end
    varargout = {all(all( abs(sf_op(sX)*c - c) <= tol ))};

case 'isinspp'
    %- check dimensions
    if size(sX.X,2) ~= size(c,1) 
        warning('Vector dim don''t match row dim : not in space !'); 
        varargout = { 0 }; return;
    end
    varargout = {all(all( abs(sf_opp(sX)*c - c) <= tol ))};
end




case {'eachinsp', 'eachinspp'}  %- each column of c in space or in dual space
%=======================================================================
% b = spm_sp('eachinsp',x,c[,tol])
% b = spm_sp('eachinspp',x,c[,tol])
%-Check whether vectors are in row/column space of X

%-Check arguments
%-----------------------------------------------------------------------
if nargin<3, error('insufficient arguments - action,x,c required'), end
c = varargin{3}; %- if isempty(c), dim wont match exept for empty sp.
if nargin<4, tol=sX.tol; else, tol = varargin{4}; end

%-Compute according to case
%-----------------------------------------------------------------------
switch lower(action)

case 'eachinsp'
    %-Check dimensions
    if size(sX.X,1) ~= size(c,1) 
        warning('Vector dim don''t match col. dim : not in space !'); 
        varargout = { 0 }; return;
    end
    varargout = {all( abs(sf_op(sX)*c - c) <= tol )};

case 'eachinspp'
    %- check dimensions
    if size(sX.X,2) ~= size(c,1) 
        warning('Vector dim don''t match row dim : not in space !'); 
        varargout = { 0 }; return;
    end
    varargout = {all( abs(sf_opp(sX)*c - c) <= tol )};
end





case '=='       % test wether two spaces are the same
%=======================================================================
% b = spm_sp('==',x1,X2)
if nargin~=3, error('too few/many input arguments - need 3');
else X2 = varargin{3}; end;

if isempty(sX.X)
   if isempty(X2), 
      warning('Both spaces empty');
        varargout = { 1 };
   else 
      warning('one space empty');
        varargout = { 0 }; 
    end;

else 
    x2 = spm_sp('Set',X2);
    maxtol = max(sX.tol,x2.tol);
    varargout = { all( spm_sp('isinsp',sX,X2,maxtol)) & ...
    all( spm_sp('isinsp',x2,sX.X,maxtol) ) };

    %- I have encountered one case where the max of tol was needed.

end;

case 'isspc'                                     %-Space structure check
%=======================================================================
% [b,str] = spm_sp('isspc',x)
if nargin~=2, error('too few/many input arguments - need 2'), end

%-Check we've been passed a structure
if ~isstruct(varargin{2}), varargout={0}; return, end

%-Go through required field names checking their existance
% (Get fieldnames once and compare: isfield doesn't work for multiple )
% (fields, and repeated calls to isfield would result in unnecessary  )
% (replicate calls to fieldnames(varargin{2}).                        )

b       = 1;
fnames  = fieldnames(varargin{2});
for str = fieldnames(sf_create)'
    b = b & any(strcmp(str,fnames));
    if ~b, break, end
end
if nargout > 1, 
    if b, str = 'ok'; else, str = 'not a space'; end;
    varargout = {b,str};
else, varargout = {b}; end;


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

case 'size'                                      %- gives the size of sX
%=======================================================================
% size = spm_sp('size',x,dim)
%
if nargin > 3, error('too many input arguments'), end
if nargin == 2, dim = []; else dim = varargin{3}; end

if ~isempty(dim)
   switch dim
      case 1, varargout = { sf_s1(sX) };
      case 2, varargout = { sf_s2(sX) };
      otherwise, error(['unknown dimension in ' action]);
   end
else %- assumes want both dim
switch nargout
   case {0,1}
      varargout = { sf_si(sX) };
   case 2
      varargout = { sf_s1(sX), sf_s2(sX) };
   otherwise
      error(['too many output arg in ' mfilename ' ' action]);
   end
end


otherwise
%=======================================================================
error(['Invalid action (',action,')'])

%=======================================================================
end % (case lower(action))


%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================
%
% The rule I tried to follow is that the space structure is accessed
% only in this sub function part : any sX.whatever should be 
% prohibited elsewhere .... still a lot to clean !!! 


function x = sf_create
%=======================================================================

x = struct(...
    'X',    [],...         % Matrix
    'tol',  [],...      % tolerance
    'ds',   [],...         % vectors of singular values 
    'u',    [],...         % u as in X = u*diag(ds)*v'
    'v',    [], ...        % v as in X = u*diag(ds)*v'
    'rk',   [],...         % rank
    'oP',   [],...      % orthogonal projector on X
    'oPp',  [],...      % orthogonal projector on X'
    'ups',  [],...      % space in which this one is embeded
    'sus',  []);           % subspace 



function x = sf_set(X)
%=======================================================================

x   = sf_create;
x.X = X;

%-- Compute the svd with svd(X,0) : find all the singular values of x.X
%-- SVD(FULL(A)) will usually perform better than SVDS(A,MIN(SIZE(A)))

%- if more column that lines, performs on X'

if size(X,1) < size(X,2)
    [x.v, s, x.u] = svd(full(X'),0);
else  
    [x.u, s, x.v] = svd(full(X),0);
end

x.ds = diag(s); clear s;

%-- compute the tolerance
x.tol =  max(size(x.X))*max(abs(x.ds))*eps;

%-- compute the rank
x.rk =  sum(x.ds > x.tol);


function x = sf_transp(x)
%=======================================================================
%
%- Tranpspose the space : note that tmp is not touched, therefore
%- only contains the address, no duplication of data is performed.

x.X     = x.X';

tmp     = x.v;
x.v     = x.u;
x.u     = tmp;

tmp     = x.oP;
x.oP    = x.oPp;
x.oPp   = tmp;
clear   tmp;


function b = sf_isset(x)
%=======================================================================
b = ~(  isempty(x.X)    |...
    isempty(x.u)    |...
    isempty(x.v)    |...
    isempty(x.ds)   |...
    isempty(x.tol)  |...
    isempty(x.rk)   );



function s1 = sf_s1(x)
%=======================================================================
s1 = size(x.X,1);

function s2 = sf_s2(x)
%=======================================================================
s2 = size(x.X,2);

function si = sf_si(x)
%=======================================================================
si = size(x.X);

function r = sf_rk(x)
%=======================================================================
r = x.rk;

function uk = sf_uk(x)
%=======================================================================
uk = x.u(:,1:sf_rk(x));

function vk = sf_vk(x)
%=======================================================================
vk = x.v(:,1:sf_rk(x));

function sk = sf_sk(x)
%=======================================================================
sk = x.ds(1:sf_rk(x));

function t = sf_t(x)
%=======================================================================
t = x.tol;

function x = sf_tol(x,t)
%=======================================================================
x(abs(x) < t) = 0;


function op = sf_op(sX)
%=======================================================================
if sX.rk > 0 
    op = sX.u(:,[1:sX.rk])*sX.u(:,[1:sX.rk])';
else  
    op = zeros( size(sX.X,1) ); 
end;

%!!!! to implement : a clever version of sf_opY (see sf_rY)


function opp = sf_opp(sX)
%=======================================================================
if sX.rk > 0 
    opp = sX.v(:,[1:sX.rk])*sX.v(:,[1:sX.rk])';
else  
    opp = zeros( size(sX.X,2) ); 
end;

%!!!! to implement : a clever version of sf_oppY (see sf_rY)


function px = sf_pinv(sX)
%=======================================================================
r = sX.rk;
if r > 0 
    px = sX.v(:,1:r)*diag( ones(r,1)./sX.ds(1:r) )*sX.u(:,1:r)';
else 
    px = zeros(size(sX.X,2),size(sX.X,1));
end

function px = sf_pinvxp(sX)
%=======================================================================
r = sX.rk;
if r > 0 
    px = sX.u(:,1:r)*diag( ones(r,1)./sX.ds(1:r) )*sX.v(:,1:r)';
else 
    px = zeros(size(sX.X));
end

function px = sf_pinvxpx(sX)
%=======================================================================
r = sX.rk;
if r > 0 
    px = sX.v(:,1:r)*diag( sX.ds(1:r).^(-2) )*sX.v(:,1:r)';
else
    px = zeros(size(sX.X,2));
end

function px = sf_jbp(sX,n)
%=======================================================================
r = sX.rk;
if r > 0 
    px = sX.v(:,1:r)*diag( sX.ds(1:r).^(n) )*sX.v(:,1:r)';
else
    px = zeros(size(sX.X,2));
end

function x = sf_jb(sX,n)
%=======================================================================
r = sX.rk;
if r > 0 
    x = sX.u(:,1:r)*diag( sX.ds(1:r).^(n) )*sX.u(:,1:r)';
else
    x = zeros(size(sX.X,1));
end

function y = sf_jbY(sX,n,Y)
%=======================================================================
r = sX.rk;
if r > 0 
    y = ( sX.u(:,1:r)*diag(sX.ds(1:r).^n) )*(sX.u(:,1:r)'*Y);
else
    y = zeros(size(sX.X,1),size(Y,2));
end
%!!!! to implement : a clever version of sf_jbY (see sf_rY)



function x = sf_cxtwdcu(sX) 
%=======================================================================
%- coordinates in sX.X -> coordinates in sX.u(:,1:rk)

x = diag(sX.ds)*sX.v';


function x = sf_cukpinvxp(sX) 
%=======================================================================
%- coordinates of pinv(sX.X') in the basis sX.u(:,1:rk)

r = sX.rk;
if r > 0 
    x = diag( ones(r,1)./sX.ds(1:r) )*sX.v(:,1:r)';
else 
    x = zeros( size(sX.X,2) );
end

function x = sf_cukx(sX) 
%=======================================================================
%- coordinates of sX.X in the basis sX.u(:,1:rk)

r = sX.rk;
if r > 0 
    x = diag( sX.ds(1:r) )*sX.v(:,1:r)';
else 
    x = zeros( size(sX.X,2) );
end


function x = sf_xpx(sX)
%=======================================================================
r = sX.rk;
if r > 0 
    x = sX.v(:,1:r)*diag( sX.ds(1:r).^2 )*sX.v(:,1:r)';
else
    x = zeros(size(sX.X,2));
end

function x = sf_xxp(sX)
%=======================================================================
r = sX.rk;
if r > 0 
    x = sX.u(:,1:r)*diag( sX.ds(1:r).^2 )*sX.u(:,1:r)';
else
    x = zeros(size(sX.X,1));
end

function x = sf_xxpY(sX,Y)
%=======================================================================
r = sX.rk;
if r > 0 
    x = sX.u(:,1:r)*diag( sX.ds(1:r).^2 )*(sX.u(:,1:r)'*Y);
else
    x = zeros(size(sX.X,1),size(Y,2));
end

function x = sf_pinvxxp(sX)
%=======================================================================
r = sX.rk;
if r > 0 
    x = sX.u(:,1:r)*diag( sX.ds(1:r).^(-2) )*sX.u(:,1:r)';
else
    x = zeros(size(sX.X,1));
end

function n = sf_n(sX)
%=======================================================================
% if the null space is in sX.v, returns it
% otherwise, performs Gramm Schmidt orthogonalisation.
% 
%
r = sX.rk;
[q p]= size(sX.X);
if r > 0
    if q >= p  %- the null space is entirely in v
        if r == p, n = zeros(p,1); else, n = sX.v(:,r+1:p); end
    else %- only part of n is in v: same as computing the null sp of sX'

        n = null(sX.X); 
%----- BUG !!!! in that part ----------------------------------------
%-      v = zeros(p,p-q); j = 1; i = 1; z = zeros(p,1);
%-      while i <= p
%-          o = z; o(i) = 1; vpoi = [sX.v(i,:) v(i,1:j-1)]';
%-          o = sf_tol(o - [sX.v v(:,1:j-1)]*vpoi,sX.tol);
%-          if any(o), v(:,j) = o/((o'*o)^(1/2)); j = j + 1; end;
%-          i = i + 1; %- if i>p, error('gramm shmidt error'); end;
%-      end
%-      n = [sX.v(:,r+1:q) v];
%--------------------------------------------------------------------
    end
else 
    n = eye(p);
end



function r = sf_r(sX)
%=======================================================================
%-
%- returns the residual forming matrix for the space sX
%- for internal use. doesn't Check whether oP exist.

r = eye(size(sX.X,1)) - sf_op(sX) ;


function Y = sf_rY(sX,Y)
%=======================================================================
% r = spm_sp('r',sX[,Y])
% 
%- tries to minimise the computation by looking whether we should do
%- I - u*(u'*Y) or n*(n'*Y) as in 'nop'

r = sX.rk;
[q p]= size(sX.X);

if r > 0 %- else returns the input;
    
    if r < q-r %- we better do I - u*u' 
        Y = Y - sX.u(:,[1:r])*(sX.u(:,[1:r])'*Y); % warning('route1');
    else
        %- is it worth computing the n ortho basis ? 
        if size(Y,2) < 5*q
            Y = sf_r(sX)*Y;            % warning('route2');
        else 
            n = sf_n(sf_transp(sX));   % warning('route3');
            Y = n*(n'*Y);
        end
    end

end


