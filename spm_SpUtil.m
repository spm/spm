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
% Quick Reference    :
%---------------------
% ('isCon',x,c)      :
% ('allCon',x,c)     :
% ('ConR',x,c)       :
% ('ConO',x,c)       :
% ('size',x,dim)     :
% ('iX0check',i0,sL) :
%---------------------
% ('i0->c',x,i0)     : Out : c
% ('c->Tsp',x,c)     : Out : [X1o [X0]]
% ('+c->Tsp',x,c)    : Out : [ukX1o [ukX0]]
% ('i0->x1o',x,i0)   : Use ('i0->c',x,i0) and ('c->Tsp',X,c)
% ('+i0->x1o',x,i0)  : Use ('i0->c',x,i0) and ('+c->Tsp',X,c)
% ('X0->c',x,X0)     :~ 
% ('+X0->c',x,cukX0) :~ 
%---------------------
% ('trRV',x[,V])     :
% ('trMV',x[,V])     :
% ('i0->edf',x,i0,V) :
%
%---------------------
%
% Improvement compared to the spm99 beta version :
%
% Improvements in df computation using spm_SpUtil('trRV',x[,V]) and
% spm_SpUtil('trMV',sX [,V]). The degrees of freedom computation requires
% in general that the trace of RV and of RVRV be computed, where R is a
% projector onto either a sub space of the design space or the residual
% space, namely the space that is orthogonal to the design space. V is
% the (estimated or assumed) variance covariance matrix and is a number
% of scans by number of scans matrix which can be huge in some cases. We
% have (thanks to S Rouquette and JB) speed up this computation 
% by using matlab built in functions of the frobenius norm and some theorems
% on trace computations. 
%
% ======================================================================
%
% FORMAT i = spm_SpUtil('isCon',x,c)
% Tests whether weight vectors specify contrasts
% x   - Design matrix X, or space structure of X
% c   - contrast matrix (I.e. matrix of contrast weights, contrasts in columns)
%       Must have column dimension matching that of X
%       [defaults to eye(size(X,2)) to test uniqueness of parameter estimates]
% i   - logical row vector indiciating estimability of contrasts in c
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
%   "Plane Answers to Complex Questions"
%    2nd Ed. Springer-Verlag, New York
%
% Andrade A, Paradis AL, Rouquette S and Poline JB, NeuroImage 9, 1999
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
%       defaults to eye(size(X,2)) to test independence of parameter estimates
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
% FORMAT c = spm_SpUtil('i0->c',x,i0)
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
% case {'x0->c'}                %- 
% FORMAT c = spm_SpUtil('X0->c',sX,X0)
%                           ----------------
%
% FORMAT [X1,X0] = spm_SpUtil('c->TSp',X,c)
% Orthogonalised partitioning of design space implied by F-contrast
% x   - Design matrix X, or space structure of X
% c   - contrast matrix (I.e. matrix of contrast weights, contrasts in columns)
%       Must have column dimension matching that of X
% X1o - contrast space - design matrix corresponding according to contrast
%       (orthogonalised wirit X0)
% X0  - matrix reduced according to null hypothesis
%       (of same size as X but rank deficient)
% FORMAT [uX1,uX0] = spm_SpUtil('c->TSp+',X,c)
%        + version to deal with the X1o and X0 partitions in the "uk basis"
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
% FORMAT X1 = spm_SpUtil('i0->x1o',X,i0)
% x   - Design matrix X, or space structure of X
% i0  - Columns of X that make up X0 - the reduced model (Ho:B1=0)
% X1  - Hypothesised null design space, i.e. that part of X orthogonal to X0
% This offers the same functionality as the 'c->TSp' option, but for
% simple reduced models formed from the columns of X.
%
% FORMAT X1 = spm_SpUtil('i0->x1o+',X,i0)
%        + version to deal with the X1o and X0 partitions in the "uk basis"
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
% OBSOLETE use FcUtil('H') for spm_SpUtil('c->H',x,c) 
% Extra sum of squares matrix O for beta's from contrast
% x   - Design matrix X, or space structure of X
% c   - contrast matrix (I.e. matrix of contrast weights, contrasts in columns)
%       Must have column dimension matching that of X
% O   - Matrix such that b'*O*b = extra sum of squares for F-test of contrast c
%
%                           ----------------
%
% OBSOLETE use spm_sp('=='...) for spm_SpUtil('c==X1o',x,c) {or 'cxpequi'}
% x   - Design matrix X, or space structure of X
% c   - contrast matrix (I.e. matrix of contrast weights, contrasts in columns)
%       Must have column dimension matching that of X
% b   - True is c is a spanning set for space of X
%       (I.e. if contrast and space test the same thing)
%
%                           ----------------
%
% FORMAT [df1,df2] = spm_SpUtil('i0->edf',x,i0,V) {or 'edf'}
% (effective) df1 and df2 the residual df for the projector onto the
% null space of x' (residual forming projector) and the numerator of
% the F-test where i0 are the columns for the null hypothesis model.
% x   - Design matrix X, or space structure of X
% i0  - Columns of X corresponding to X0 partition X = [X1,X0] & with
%       parameters B = [B1;B0]. Ho:B1=0
% V   - V matrix
%
%                           ----------------
%
% FORMAT  sz           = spm_SpUtil('size',x,dim)
% FORMAT [sz1,sz2,...] = spm_SpUtil('size',x)
% Returns size of design matrix
% (Like MatLab's `size`, but copes with design matrices inside structures.)
% x   - Design matrix X, or structure containing design matrix in field X
%       (Structure needn't be a space structure.)
% dim - dimension which to size
% sz  - size
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes Jean-Baptiste Poline
% $Id: spm_SpUtil.m 4137 2010-12-15 17:18:32Z guillaume $

% (frobenius norm trick by S. Rouquette)

%-Format arguments
%-----------------------------------------------------------------------
if nargin==0, error('do what? no arguments given...')
else action = varargin{1}; end


switch lower(action), 

case {'iscon','allcon','conr','cono'}
%=======================================================================
% i = spm_SpUtil('isCon',x,c)
if nargin==0, varargout={[]}; error('isCon : no argument specified'), end;
if nargin==1, 
   varargout={[]}; warning('isCon : no contrast specified'); return; 
end;
if ~spm_sp('isspc',varargin{2})
     sX = spm_sp('Set',varargin{2});
else sX = varargin{2}; end
if nargin==2, c=eye(spm_sp('size',sX,2)); else c=varargin{3}; end;
if isempty(c), varargout={[]}; return, end

switch lower(action)
    case 'iscon'
        varargout = { spm_sp('eachinspp',sX,c) };
    case 'allcon'
        varargout = {spm_sp('isinspp',sX,c)};
    case 'conr'
        if size(c,1) ~= spm_sp('size',sX,2) 
            error('Contrast not of the right size'), end
        %-Compute inner products of data weight vectors
        % (c'b = c'*pinv(X)*Y = w'*Y
        % (=> w*w' = c'*pinv(X)*pinv(X)'*c == c'*pinv(X'*X)*c
        r   = c'*spm_sp('XpX-',sX)*c;
        %-normalize by "cov(r)" to get correlations
        r   = r./(sqrt(diag(r))*sqrt(diag(r))');
        r(abs(r) < sX.tol)=0;       %-set near-zeros to zero
        varargout = {r};                %-return r
    case 'cono'
        %-This is the same as ~spm_SpUtil('ConR',x,c), and so returns
        % the contrast orthogonality (though not their corelations).
        varargout = { abs(c'* spm_sp('XpX',sX) *c) < sX.tol};
end


case {'+c->tsp','c->tsp'}    %- Ortho. partitioning implied by F-contrast
%=======================================================================
% spm_SpUtil('c->Tsp',sX,c)
% + version of 'c->tsp'. 
% The + version returns the same in the base u(:,1:r).

%--------- begin argument check ------------------------------
if nargin ~= 3, error(['Wrong number of arguments in ' action])
else sX = varargin{2}; c = varargin{3}; end;
if nargout > 2, error(['Too many output arguments in ' action]), end;
if ~spm_sp('isspc',sX),  sX = spm_sp('set',varargin{2}); end;
if sX.rk == 0, error('c->Tsp null rank sX == 0'); end;
if ~isempty(c) && spm_sp('size',sX,2) ~= size(c,1), 
   error(' c->TSp matrix & contrast dimensions don''t match');   
end
%--------- end argument check --------------------------------- 

%- project c onto the space of  X' if needed 
%-------------------------------------------
if ~isempty(c) && ~spm_sp('isinspp',sX,c),
   warning([sprintf('\n') 'c is not a proper contrast in ' action ...
            ' in ' mfilename sprintf('\n') '!!! projecting...' ]);
   disp('from'), c, disp('to'), c = spm_sp('oPp:',sX,c)
end

cukFlag = strcmp(lower(action),'+c->tsp');

switch nargout
% case 0
%   warning(['no output demanded in ' mfilename ' ' action])
case {0,1}
   if ~isempty(c) && any(any(c))               %- c not empty & not null
      if cukFlag, varargout = { spm_sp('cukxp-:',sX,c) };
      else        varargout = { spm_sp('xp-:',sX,c) };
      end
   else if isempty(c), varargout = { [] };     %- c empty
       else                                    %- c null
           if cukFlag, varargout = { spm_sp('cukx',sX,c) }; 
           else        varargout = { spm_sp('x',sX)*c };
           end
       end
   end

case 2
   if ~isempty(c) && any(any(c))         %- not empty and not null
      if cukFlag,
         varargout = {  
        spm_sp('cukxp-:',sX,c), ...      %- X1o
        spm_sp('cukx',sX,spm_sp('r',spm_sp('set',c))) };      %- X0
      else
         varargout = {  
        spm_sp(':',sX, spm_sp('xp-:',sX,c)), ...              %- X1o
        spm_sp(':',sX, ...
               spm_sp('x',sX)*spm_sp('r',spm_sp('set',c))) }; %- X0
      end
   else            
      if isempty(c),                     %- empty
         if cukFlag, varargout = { [], spm_sp('cukx',sX) }; 
         else        varargout = { [], spm_sp('x',sX) };
         end
      else                               %- null
         if cukFlag, 
            varargout = { spm_sp(':',sX,spm_sp('cukx',sX,c)), ...
                          spm_sp(':',sX,spm_sp('cukx',sX)) }; 
         else
            varargout = { spm_sp('x',sX)*c, spm_sp('x',sX)};
         end
      end;
   end

otherwise
   error(['wrong number of output argument in ' action]);
end


case {'i0->x1o','+i0->x1o'} %- Space tested whilst keeping size of X(i0)
%=======================================================================
% X1o = spm_SpUtil('i0->X1o',sX,i0)

% arguments are checked in calls to spm_Util
%--------------------------------------------
if nargin<3, error('Insufficient arguments'),  
else sX = varargin{2}; i0 = varargin{3}; end;

cukFlag = strcmp(lower(action),'+i0->x1o');

c  = spm_SpUtil('i0->c',sX,i0);

if cukFlag, 
   varargout = { spm_SpUtil('+c->TSp',sX,c) };
else 
   varargout = { spm_SpUtil('c->TSp',sX,c) };
end


case {'i0->c'}                 %- 
%=======================================================================
% c = spm_SpUtil('i0->c',sX,i0)
%
% if i0 == [] : returns a proper contrast
% if i0 == [1:size(sX.X,2)] : returns [];
%
%- Algorithm : avoids the pinv(X0) and insure consistency
%- Get the estimable parts of c0 and c1
%- remove from c1_estimable the estimable part of c0.
%- Use the rotation making X1o orthog. to X0.
%- i0 is checked when Fc is created
%- If i0 defines a space that is the space of X but with
%- fewer vectors, c is null. 

%--------- begin argument check --------------------------------
if nargin<3, error('Insufficient arguments'),  
else sX = varargin{2}; i0 = varargin{3}; end;
if ~spm_sp('isspc',sX),  sX = spm_sp('set',varargin{2}); end;
if spm_sp('rk',sX) == 0, error('i0->c null rank sX == 0'); end;
sL  = spm_sp('size',sX,2);
i0  = sf_check_i0(i0,sL);
%--------- end argument check ---------------------------------- 

c0  = eye(sL); c0 = c0(:,i0);
c1  = eye(sL); c1 = c1(:,setdiff(1:sL,i0));

%- try to avoid the matlab error when doing svd of matrices with
%- high multiplicities. (svd convergence pb)

if ~ spm_sp('isinspp',sX,c0), c0 = spm_sp('oPp:',sX,c0); end;
if ~ spm_sp('isinspp',sX,c1), c1 = spm_sp('oPp:',sX,c1); end;

if ~isempty(c1)
   if ~isempty(c0) 
      %- varargout = { spm_sp('res',spm_sp('set',opp*c0),opp*c1) };
      %- varargout = { c1 - c0*pinv(c0)*c1 }; NB: matlab pinv uses
      %- svd: will fail if spm_sp('set')  fails. 

      varargout = { spm_sp('r:',spm_sp('set',c0),c1) };

   else varargout = { spm_sp('xpx',sX) }; end;
else
   varargout = { [] };  %- not zeros(sL,1) : this is return when  
            %- appropriate
end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
case {'+x0->c','x0->c'}             %- 
%=======================================================================
% c = spm_SpUtil('X0->c',sX,X0)
% c = spm_SpUtil('+X0->c',sX,cukX0)
% + version of 'x0->c'. 
% The + version returns the same in the base u(:,1:r).

warning('Not tested for release - provided for completeness');
cukFlag = strcmp(lower(action),'+x0->c');

%--------- begin argument check --------- 
if nargin<3, error('Insufficient arguments'),  
else
   sX = varargin{2}; 
   if cukFlag, cukX0 = varargin{3}; else X0 = varargin{3}; end 
end
if ~spm_sp('isspc',sX),  sX = spm_sp('set',varargin{2}); end
if spm_sp('rk',sX) == 0, error(['null rank sX == 0 in ' action]); end
if cukFlag
   if ~isempty(cukX0) && spm_sp('rk',sX) ~= size(cukX0,1),
      cukX0, spm_sp('rk',sX),
      error(['cukX0 of wrong size ' mfilename ' ' action]), end
else
   if ~isempty(X0) && spm_sp('size',sX,1) ~= size(X0,1),
      X0, spm_sp('size',sX,1),
      error(['X0 of wrong size ' mfilename ' ' action]),X0, end
end
%--------- end argument check --------- 

if cukFlag
   if isempty(cukX0), X0 = []; else X0 = spm_sp('ox',sX)*cukX0; end
end

varargout = { sf_X0_2_c(X0,sX) };



case {'c->h','betarc'}  %-Extra sum of squares matrix for beta's from 
                        %- contrast : use F-contrast if possible
%=======================================================================
% H = spm_SpUtil('c->H',sX,c)
error(' Obsolete : Use F-contrast utilities  ''H'' or ''Hsqr''... ');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



%=======================================================================
%=======================================================================
%       trace part
%=======================================================================
%=======================================================================
%

case 'trrv'                      %-Traces for (effective) df calculation
%=======================================================================
% [trRV,trRVRV]= spm_SpUtil('trRV',x[,V])

if nargin == 1, error('insufficient arguments');
else sX = varargin{2};  end;

if ~spm_sp('isspc',sX), sX = spm_sp('Set',sX); end;

rk = spm_sp('rk',sX);
sL = spm_sp('size',sX,1);

if  sL == 0,
    warning('space with no dimension ');    
    if nargout==1, varargout = {[]};
    else varargout = {[], []}; end
else

   if nargin > 2 && ~isempty(varargin{3})    

    V = varargin{3};
    u = sX.u(:,1:rk);
    clear sX;
    if nargout==1
        %-only trRV needed
        if rk==0 || isempty(rk),  trMV = 0; 
        else trMV = sum(sum( u .* (V*u) ));
        end
        varargout = { trace(V) - trMV};
    else
        %-trRVRV is needed as well
        if rk==0 || isempty(rk),  
            trMV = 0; 
            trRVRV = (norm(V,'fro'))^2;
         trV = trace(V);
         clear V u
        else 
           Vu = V*u;
           trV = trace(V);
           trRVRV = (norm(V,'fro'))^2;
                   clear V; 
           trRVRV = trRVRV - 2*(norm(Vu,'fro'))^2;
           trRVRV = trRVRV + (norm(u'*Vu,'fro'))^2;
           trMV = sum(sum( u .* Vu ));
                   clear u Vu
        end
        varargout = {(trV - trMV), trRVRV};
    end
           
   else  %- nargin == 2 | isempty(varargin{3})

    if nargout==1
        if rk==0 || isempty(rk), varargout = {sL}; 
        else varargout = {sL - rk}; 
        end
    else
        if rk==0 || isempty(rk),  varargout = {sL,sL};
        else varargout = {sL - rk, sL - rk};
        end    
    end

    end
end 


case 'trmv'                     %-Traces for (effective) Fdf calculation
%=======================================================================
% [trMV, trMVMV]] = spm_SpUtil('trMV',sX [,V])
%
% NB : When V is given empty, the routine asssumes it's identity
% This is used in spm_FcUtil.

if nargin == 1, error('insufficient arguments');
else sX = varargin{2}; end;
if ~spm_sp('isspc',sX), sX = spm_sp('Set',sX); end;
rk = spm_sp('rk',sX);

if isempty(rk)
    warning('Rank is empty');   
    if nargout==1, varargout = {[]};
    else varargout = {[], []}; end
    return; 
elseif  rk==0, warning('Rank is null in spm_SpUtil trMV ');
    if nargout==1, varargout = {0};
    else varargout = {0, 0}; end
    return; 
end;

if nargin > 2 && ~isempty(varargin{3}) %- V provided, and assumed correct !

    V = varargin{3};
    u = sX.u(:,1:rk);
    clear sX;

    if nargout==1
        %-only trMV needed
        trMV = sum(sum(u' .* (u'*V) ));
        varargout = {trMV};
    else 
        %-trMVMV is needed as well
        Vu = V*u;
                clear V
        trMV = sum(sum( u .* Vu ));
        trMVMV = (norm(u'*Vu,'fro'))^2;
            clear u Vu
        varargout = {trMV, trMVMV};
    end

else  % nargin == 2 | isempty(varargin{3}) %-no V specified: trMV == trMVMV
    if nargout==1
        varargout = {rk};
    else
        varargout = {rk, rk};
    end
end
 


case {'i0->edf','edf'}                  %-Effective F degrees of freedom
%=======================================================================
% [df1,df2] = spm_SpUtil('i0->edf',sX,i0,V)
%-----------------------------------------------------------------------
%--------- begin argument check ---------------------------------------- 
if nargin<3, error('insufficient arguments'),
else i0 = varargin{3}; sX = varargin{2}; end
if ~spm_sp('isspc',sX), sX = spm_sp('Set',sX); end;
i0  = sf_check_i0(i0,spm_sp('size',sX,2));
if nargin == 4, V=varargin{4}; else V = eye(spm_sp('size',sX,1)); end;
if nargin>4, error('Too many input arguments'), end;
%--------- end argument check ------------------------------------------ 

warning(' Use F-contrast utilities if possible ... ');

[trRV,trRVRV]    = spm_SpUtil('trRV', sX, V);
[trMpV,trMpVMpV] = spm_SpUtil('trMV',spm_SpUtil('i0->x1o',sX, i0),V);
varargout = {trMpV^2/trMpVMpV, trRV^2/trRVRV};


%=======================================================================
%=======================================================================
%       Utilities
%=======================================================================
%=======================================================================


case 'size'                                      %-Size of design matrix
%=======================================================================
% sz = spm_SpUtil('size',x,dim)

if nargin<3, dim=[]; else dim = varargin{3}; end
if nargin<2, error('insufficient arguments'), end

if isstruct(varargin{2})
   if isfield(varargin{2},'X')
    sz = size(varargin{2}.X);
   else error('no X field'); end;  
else
    sz = size(varargin{2});
end

if ~isempty(dim)
    if dim>length(sz), sz = 1; else sz = sz(dim); end
    varargout = {sz};
elseif nargout>1
    varargout = cell(1,min(nargout,length(sz)));
    for i=1:min(nargout,length(sz)), varargout{i} = sz(i); end
else
    varargout = {sz};
end


case 'ix0check'                                      %-
%=======================================================================
% i0c = spm_SpUtil('iX0check',i0,sL)

if nargin<3, error('insufficient arguments'),
else i0 = varargin{2}; sL = varargin{3}; end;

varargout = {sf_check_i0(i0,sL)};


otherwise
%=======================================================================
error('Unknown action string in spm_SpUtil')

%=======================================================================
end




%=======================================================================
function i0c = sf_check_i0(i0,sL)
% NB : [] = sf_check_i0([],SL);
%

if all(ismember(i0,[0,1])) && length(i0(:))==sL, i0c=find(i0); 
elseif ~isempty(i0) && any(floor(i0)~=i0) || any(i0<1) || any(i0>sL)
    error('logical mask or vector of column indices required')
else i0c = i0; end

%=======================================================================
function c = sf_X0_2_c(X0,sX)  
%
%- Algorithm to avoids the pinv(X0) and insure consistency
%- Get a contrast that span the space of X0 and is estimable
%- Get the orthogonal complement and project onto the estimable space
%- Strip zeros columns and use the rotation making X1o orthog. to X0
% !!! tolerance dealing ?

if ~isempty(X0)

   sc0 = spm_sp('set',spm_sp('x-',sX,X0));
   if  sc0.rk
      c = spm_sp('oPp:',sX,spm_sp('r',sc0));
   else 
      c = spm_sp('oPp',sX);
   end;
   c  = c(:,any(c));
   sL = spm_sp('size',sX,2);

   %- why the   "& size(X0,2) ~= sL"   !!!?
   if isempty(c) && size(X0,2) ~= sL
      c = zeros(sL,1);
   end

else 
   c = spm_sp('xpx',sX);
end

%- c = spm_sp('r',sc0,spm_sp('oxp',sX)); would also works.
