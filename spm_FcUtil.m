function varargout = spm_FcUtil(varargin)
% Contrast utilities
% FORMAT varargout = spm_FcUtil(action,varargin)
%_______________________________________________________________________
%
% spm_FcUtil is a multi-function function containing various utilities
% for contrast construction and manipulation. In general, it accepts
% design matrices as plain matrices or as space structures setup by
% spm_sp (that is preferable in general).
% 
% The use of spm_FcUtil should help with robustness issues and
% maintainability of SPM.  % Note that when space structures are passed
% as arguments is is assummed that their basic fields are filled in.
% See spm_sp for details of (design) space structures and their
% manipulation.
%
%
% ======================================================================
% case 'fconfields'             %- fields of F contrast
% Fc = spm_FcUtil('FconFields')
%
%- simply returns the fields of a contrast structure.
%
%=======================================================================
% case 'set'                    %- Create an F contrast
% Fc = spm_FcUtil('Set',name, STAT, set_action, value, sX)
%
%- Set will fill in the contrast structure, in particular 
%- c (in the contrast space), X1o (the space actually tested) and
%- X0 (the space left untested), such that space([X1o X0]) == sX.
%- STAT is either 'F' or 'T';
%- name is a string descibing the contrast.
%
%- There are three ways to set a contrast :
%- set_action is 'c','c+'   : value can then be zeros.
%-                dimensions are in X', 
%-                f c+ is used, value is projected onto sX';
%-                                iX0 is set to 'c' or 'c+';
%- set_action is 'iX0'      : defines the indices of the columns 
%-                that will not be tested. Can be empty.
%- set_action is 'X0'       : defines the space that will remain 
%-                unchanged. The orthogonal complement is
%-                tested; iX0 is set to 'X0';
%-                                  
%=======================================================================
% case 'isfcon'                 %- Is it an F contrast ?
% b = spm_FcUtil('IsFcon',Fc)
%
%=======================================================================
% case 'fconedf'                    %- F contrast edf
% [edf_tsp edf_Xsp] = spm_FcUtil('FconEdf', Fc, sX [, V])
%
%- compute the effective degrees of freedom of the numerator edf_tsp
%- and (optionally) the denominator edf_Xsp of the contrast.
%
%=======================================================================
% case 'hsqr' %-Extra sum of squares sqr matrix for beta's from contrast
% hsqr = spm_FcUtil('Hsqr',Fc, sX)
%
%- This computes the matrix hsqr such that a the numerator of an F test
%- will be beta'*hsqr'*hsqr*beta
%
%=======================================================================
% case 'h'     %-Extra sum of squares matrix for beta's from contrast
% H = spm_FcUtil('H',Fc, sX)
%
%- This computes the matrix H such that a the numerator of an F test
%- will be beta'*H*beta
%-                                  
%=======================================================================
% case 'yc' %- Fitted data corrected for confounds defined by Fc 
% Yc = spm_FcUtil('Yc',Fc, sX, b)
%
%- Input : b : the betas                            
%- Returns the corrected data Yc for given contrast. Y = Yc + Y0 + error
%
%=======================================================================
% case 'y0' %-  Confounds data defined by Fc 
% Y0 = spm_FcUtil('Y0',Fc, sX, b)
%
%- Input : b : the betas                            
%- Returns the confound data Y0 for a given contrast. Y = Yc + Y0 + error
%
%=======================================================================
% case {'|_'}       %-  Fc orthogonalisation 
% Fc = spm_FcUtil('|_',Fc1, sX, Fc2)
%
%- Orthogonolise a (list of) contrasts Fc1 wrt a (list of) contrast Fc2
%- such that the space these contrasts test are orthogonal.
%- If contrasts are not estimable contrasts, works with the estimable 
%- part. In any case, returns estimable contrasts.  
%
%=======================================================================
% case {'|_?'}      %-  Are contrasts orthogonals 
% b = spm_FcUtil('|_?',Fc1, sX [, Fc2])
%
%- Tests whether a (list of) contrast is orthogonal. Works with the
%- estimable part if they are not estimable. With only one argument,
%- tests whether the list is made of orthogonal contrasts. With Fc2
%- provided, tests whether the two (list of) contrast are orthogonal. 
%
%=======================================================================
% case 'in'    %-  Fc1 is in list of  contrasts Fc2
% [iFc2 iFc1] = spm_FcUtil('In', Fc1, sX, Fc2)
%
%- Tests wether a (list of) contrast Fc1 is in a list of contrast Fc2.
%- returns the indices iFc2 where element of Fc1 have been found
%- in Fc2 and the indices iFc1 of the element of Fc1 found in Fc2.
%- These indices are not necessarily unique.
%
%=======================================================================
% case '~unique'     %-  Fc list unique 
% idx = spm_FcUtil('~unique', Fc, sX)
%
%- returns indices ofredundant contrasts in Fc
%- such that Fc(idx) = [] makes Fc unique.
%
%=======================================================================
% case {'0|[]','[]|0'}     %-  Fc is null or empty 
% b = spm_FcUtil('0|[]', Fc, sX)
%
%- NB : for the "null" part, checks if the contrast is in the null space 
%- of sX (completely non estimable !)
%=======================================================================
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jean-Baptiste Poline
% $Id: spm_FcUtil.m 4137 2010-12-15 17:18:32Z guillaume $


%-Format arguments
%-----------------------------------------------------------------------
if nargin==0, error('do what? no arguments given...')
else action = varargin{1}; end


switch lower(action),

case 'fconfields'               %- fields of F contrast
%=======================================================================
% Fc = spm_FcUtil('FconFields')

if nargout > 1, error('Too many output arguments: FconFields'), end
if nargin  > 1, error('Too many input arguments:  FconFields'), end

varargout = {sf_FconFields};

case {'set','v1set'}                %- Create an F contrast
%=======================================================================
% Fc = spm_FcUtil('Set',name, STAT, set_action, value, sX)
%
% Sets the contrast structure with set_action either 'c', 'X0' or 'iX0'
% resp. for a contrast, the null hyp. space or the indices of which.
% STAT can be 'T' or 'F'.
%
% If not set by iX0 (in which case field .iX0 containes the indices),
% field .iX0 is set as a string containing the set_action: {'X0','c','c+','ukX0'}
%
% if STAT is T, then set_action should be 'c' or 'c+' 
% (at the moment, just a warning...)
% if STAT is T and set_action is 'c' or 'c+', then 
% checks whether it is a real T. 
%
% 'v1set' is NOT provided for backward compatibility so far ...

%-check # arguments...
%--------------------------------------------------------------------------
if nargin<6, error('insufficient arguments'), end
if nargout > 1, error('Too many output arguments Set'), end

%-check arguments...
%--------------------------------------------------------------------------
if ~ischar(varargin{2}), error('~ischar(name)'), end
if ~(varargin{3}=='F'||varargin{3}=='T'||varargin{3}=='P'), 
    error('~(STAT==F|STAT==T|STAT==P)'), end
if ~ischar(varargin{4}), error('~ischar(varargin{4})'); 
else set_action = varargin{4}; end

sX = varargin{6};
if ~spm_sp('isspc',sX), sX = spm_sp('set',sX);  end
if isempty(sX.X), error('Empty space X in Set'); end

Fc = sf_FconFields;
%- use the name as a flag to insure that F-contrast has been 
%- properly created;
Fc.name = varargin{2};

Fc.STAT = varargin{3};
if Fc.STAT=='T' &&  ~(any(strcmp(set_action,{'c+','c'}))) 
   warning('enter T stat with contrast - here no check rank == 1');
end

[sC sL] = spm_sp('size',sX);

%- allow to define the contrast the old (version 1) way ?
%- NO. v1 = strcmp(action,'v1set');

switch set_action,
    case {'c','c+'}
        Fc.iX0  = set_action;
        c = spm_sp(':', sX, varargin{5});
        if isempty(c)
            [Fc.X1o.ukX1o Fc.X0.ukX0] = spm_SpUtil('+c->Tsp',sX,[]);
            %- v1 [Fc.X1o Fc.X0] = spm_SpUtil('c->Tsp',sX,[]);
            Fc.c    = c;
        elseif size(c,1) ~= sL,
            error(['not contrast dim. in ' mfilename ' ' set_action]);
        else
            if strcmp(set_action,'c+')
                if ~spm_sp('isinspp',sX,c), c = spm_sp('oPp:',sX,c); end
            end;
            if Fc.STAT=='T' &&  ~sf_is_T(sX,c)
                %- Could be make more self-correcting by giving back an F
                error('trying to define a t that looks like an F');
            end
            Fc.c   = c;
            [Fc.X1o.ukX1o Fc.X0.ukX0] = spm_SpUtil('+c->Tsp',sX,c);
            %- v1 [Fc.X1o Fc.X0] = spm_SpUtil('c->Tsp',sX,c);
        end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 'option given for completeness - not for SPM use'

    case {'X0'}
        warning('option given for completeness - not for SPM use');
        Fc.iX0  = set_action;
        X0 = spm_sp(':', sX, varargin{5});
        if isempty(X0),
            Fc.c         = spm_sp('xpx',sX);
            Fc.X1o.ukX1o = spm_sp('cukx',sX);
            Fc.X0.ukX0   = [];
        elseif size(X0,1) ~= sC,
            error('dimension of X0 wrong in Set');
        else
            Fc.c         = spm_SpUtil('X0->c',sX,X0);
            Fc.X0.ukX0   = spm_sp('ox',sX)'*X0;
            Fc.X1o.ukX1o = spm_SpUtil('+c->Tsp',sX,Fc.c);
        end
        
    case 'ukX0'
        warning('option given for completeness - not for SPM use');
        Fc.iX0  = set_action;
        if isempty(ukX0),
            Fc.c         = spm_sp('xpx',sX);
            Fc.X1o.ukX1o = spm_sp('cukx',sX);
            Fc.X0.ukX0   = [];
        elseif size(ukX0,1) ~= spm_sp('rk',sX),
            error('dimension of cukX0 wrong in Set');
        else
            Fc.c         = spm_SpUtil('+X0->c',sX,ukX0);
            Fc.X0.ukX0   = ukX0;
            Fc.X1o.ukX1o = spm_SpUtil('+c->Tsp',sX,Fc.c);
        end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    case 'iX0'
        iX0        = varargin{5};
        iX0        = spm_SpUtil('iX0check',iX0,sL);
        Fc.iX0     = iX0;
        Fc.X0.ukX0 = spm_sp('ox',sX)' * spm_sp('Xi',sX,iX0);
        if isempty(iX0),
            Fc.c          = spm_sp('xpx',sX);
            Fc.X1o.ukX1o  = spm_sp('cukx',sX);
        else
            Fc.c          = spm_SpUtil('i0->c',sX,iX0);
            Fc.X1o.ukX1o  = spm_SpUtil('+c->Tsp',sX,Fc.c);
        end
        
    otherwise
        error('wrong action in Set ');
        
end
varargout = {Fc};

case 'x0' % spm_FcUtil('X0',Fc,sX)
%=======================================================================
if nargin ~= 3, 
   error('too few/many input arguments - need 2');
else 
   Fc = varargin{2}; sX = varargin{3};
end
if nargout ~= 1, error('too few/many output arguments - need 1'), end
if ~sf_IsFcon(Fc), error('argument is not a contrast struct'), end
if ~spm_sp('isspc',sX), sX = spm_sp('set',sX); end

varargout = {sf_X0(Fc,sX)};

case 'x1o' % spm_FcUtil('X1o',Fc,sX)
%=======================================================================
if nargin ~= 3, 
   error('too few/many input arguments - need 2');
else 
   Fc = varargin{2}; sX = varargin{3};
end
if nargout ~= 1, error('too few/many output arguments - need 1'), end
if ~sf_IsFcon(Fc), error('argument is not a contrast struct'), end
if ~spm_sp('isspc',sX), sX = spm_sp('set',sX); end

varargout = {sf_X1o(Fc,sX)};


case 'isfcon'               %- Is it an F contrast ?
%=======================================================================
% yes_no = spm_FcUtil('IsFcon',Fc)

if nargin~=2, error('too few/many input arguments - need 2'), end
if ~isstruct(varargin{2}), varargout={0}; 
else varargout = {sf_IsFcon(varargin{2})};
end


case 'fconedf'              %- F contrast edf
%=======================================================================
% [edf_tsp edf_Xsp] = spm_FcUtil('FconEdf', Fc, sX [, V])

if nargin<3, error('Insufficient arguments'), end
if nargout >= 3, error('Too many output argument.'), end
Fc = varargin{2};
sX = varargin{3};
if nargin == 4, V = varargin{4}; else V = []; end

if ~sf_IsFcon(Fc), error('Fc must be Fcon'), end
if ~spm_sp('isspc',sX)
    sX = spm_sp('set',sX);  end

if ~sf_isempty_X1o(Fc)
   [trMV, trMVMV] = spm_SpUtil('trMV',sf_X1o(Fc,sX),V);
else
   trMV   = 0;
   trMVMV = 0;
end 
if ~trMVMV, edf_tsp = 0; warning('edf_tsp = 0'), 
else  edf_tsp = trMV^2/trMVMV; end;    

if nargout == 2

   [trRV, trRVRV] = spm_SpUtil('trRV',sX,V);
   if ~trRVRV, edf_Xsp = 0; warning('edf_Xsp = 0'),
   else  edf_Xsp = trRV^2/trRVRV; end;

   varargout = {edf_tsp, edf_Xsp};
else    
   varargout = {edf_tsp};
end



%=======================================================================
%=======================================================================
%       parts that use F contrast
%=======================================================================
%=======================================================================
%
% Quick reference        : L : can be lists of ...
%-------------------------
% ('Hsqr',Fc, sX)        : Out: Hsqr / ESS = b' * Hsqr' * Hsqr * b
% ('H',Fc, sX)           : Out: H / ESS = b' * H * b
% ('Yc',Fc, sX, b)       : Out: Y corrected = X*b - X0*X0- *Y 
% ('Y0',Fc, sX, b)       : Out: Y0 = X0*X0- *Y
% ('|_',LFc1, sX, LFc2)  : Out: Fc1 orthog. wrt Fc2 
% ('|_?',LFc1,sX [,LFc2]): Out: is Fc2 ortho to Fc1 or is Fc1 ortho ? 
% ('In', LFc1, sX, LFc2) : Out: indices of Fc2 if "in", 0 otherwise
% ('~unique', LFc, sX)   : Out: indices of redundant contrasts 
% ('0|[]', Fc, sX)       : Out: 1 if Fc is zero or empty, 0 otherwise


case 'hsqr'     %-Extra sum of squares matrix for beta's from contrast
%=======================================================================
% hsqr = spm_FcUtil('Hsqr',Fc, sX)

if nargin<3, error('Insufficient arguments'), end
if nargout>1, error('Too many output argument.'), end
Fc = varargin{2};
sX = varargin{3};

if ~sf_IsFcon(Fc), error('Fc must be F-contrast'), end
if ~sf_IsSet(Fc), error('Fcon must be set'); end; %-
if ~spm_sp('isspc',sX), sX = spm_sp('set',sX);  end;

if sf_isempty_X1o(Fc)
    if ~sf_isempty_X0(Fc)
        %- assumes that X0 is sX.X
        %- warning(' Empty X1o in spm_FcUtil(''Hsqr'',Fc,sX) ');
        varargout = { zeros(1,spm_sp('size',sX,2)) };
    else    
        error(' Fc must be set ');
    end
else
    varargout = { sf_Hsqr(Fc,sX) };
end


case 'h'     %-Extra sum of squares matrix for beta's from contrast
%=======================================================================
% H = spm_FcUtil('H',Fc, sX)

% Empty and zeros dealing : 
% This routine never returns an empty matrix. 
% If sf_isempty_X1o(Fc) | isempty(Fc.c) it explicitly 
% returns a zeros projection matrix.

if nargin<2, error('Insufficient arguments'), end
if nargout>1, error('Too many output argument.'), end
Fc = varargin{2};
sX = varargin{3};

if ~sf_IsFcon(Fc), error('Fc must be F-contrast'), end
if ~sf_IsSet(Fc), error('Fcon must be set'); end; %-
if ~spm_sp('isspc',sX), sX = spm_sp('set',sX);  end;

if sf_isempty_X1o(Fc)
    if ~sf_isempty_X0(Fc)
        %- assumes that X0 is sX.X
        %- warning(' Empty X1o in spm_FcUtil(''H'',Fc,sX) ');
        varargout = { zeros(spm_sp('size',sX,2)) };
    else    
        error(' Fc must be set ');
    end
else
    varargout = { sf_H(Fc,sX) };  
end




case 'yc'     %- Fitted data corrected for confounds defined by Fc 
%=======================================================================
% Yc = spm_FcUtil('Yc',Fc, sX, b)

if nargin < 4, error('Insufficient arguments'), end
if nargout > 1, error('Too many output argument.'), end

Fc = varargin{2}; sX = varargin{3}; b = varargin{4};

if ~sf_IsFcon(Fc), error('Fc must be F-contrast'), end
if ~sf_IsSet(Fc), error('Fcon must be set'); end; 
if ~spm_sp('isspc',sX), sX = spm_sp('set',sX);  end;
% if ~spm_FcUtil('Rcompatible',Fc,sX), ...
%   error('sX and Fc must be compatible'), end;
if spm_sp('size',sX,2) ~= size(b,1), 
    error('sX and b must be compatible'), end;

if sf_isempty_X1o(Fc)
    if ~sf_isempty_X0(Fc)
       %- if space of interest empty or null, returns zeros !
       varargout = { zeros(spm_sp('size',sX,1),size(b,2)) };
    else    
       error(' Fc must be set ');
    end
else
    varargout = { sf_Yc(Fc,sX,b) };  
end



case 'y0'     %- Fitted data corrected for confounds defined by Fc 
%=======================================================================
% Y0 = spm_FcUtil('Y0',Fc, sX, b)

if nargin < 4, error('Insufficient arguments'), end
if nargout > 1, error('Too many output argument.'), end

Fc = varargin{2}; sX = varargin{3}; b = varargin{4};

if ~sf_IsFcon(Fc), error('Fc must be F-contrast'), end
if ~sf_IsSet(Fc), error('Fcon must be set'); end; 
if ~spm_sp('isspc',sX), sX = spm_sp('set',sX);  end;
if spm_sp('size',sX,2) ~= size(b,1), 
    error('sX and b must be compatible'), end;

if sf_isempty_X1o(Fc)
    if ~sf_isempty_X0(Fc)
       %- if space of interest empty or null, returns zeros !
       varargout = { sX.X*b };
    else    
       error(' Fc must be set ');
    end
else
    varargout = { sf_Y0(Fc,sX,b) };  
end


case {'|_'}     %-  Fc orthogonalisation 
%=======================================================================
% Fc = spm_FcUtil('|_',Fc1, sX, Fc2)

% returns Fc1 orthogonolised wrt Fc2 

if nargin < 4, error('Insufficient arguments'), end
if nargout > 1, error('Too many output argument.'), end
Fc1 = varargin{2}; sX = varargin{3}; Fc2 = varargin{4}; 

%-check arguments
%-----------------------------------------------------------------------
L1 = length(Fc1);
if ~L1, warning('no contrast given to |_'); varargout = {[]}; return; end
for i=1:L1
    if ~sf_IsFcon(Fc1(i)), error('Fc1(i) must be a contrast'), end
end
L2 = length(Fc2);
if ~L2, error('must have at least a contrast in Fc2'); end
for i=1:L2
    if ~sf_IsFcon(Fc2(i)), error('Fc2(i) must be a contrast'), end
end
if ~spm_sp('isspc',sX), sX = spm_sp('set',sX);  end;

%-create an F-contrast for all the Fc2
%--------------------------------------------------------------------------
str  = Fc2(1).name; for i=2:L2 str = [str ' ' Fc2(i).name]; end;
Fc2  = spm_FcUtil('Set',str,'F','c+',cat(2,Fc2(:).c),sX);

if sf_isempty_X1o(Fc2) || sf_isnull(Fc2,sX)
    varargout = {Fc1};
else
    for i=1:L1
        if sf_isempty_X1o(Fc1(i)) || sf_isnull(Fc1(i),sX)
            %- Fc1(i) is an [] or 0 contrast : ortho to anything; 
            out(i) = Fc1(i);
        else
            out(i) = sf_fcortho(Fc1(i), sX, Fc2);
        end
    end
    varargout = {out};  
end

case {'|_?'}        %-  Are contrasts orthogonals 
%=======================================================================
% b = spm_FcUtil('|_?',Fc1, sX [, Fc2])

if nargin < 3, error('Insufficient arguments'), end
Fc1 = varargin{2}; sX = varargin{3};
if nargin > 3, Fc2 = varargin{4}; else Fc2 = []; end;
if isempty(Fc1), error('give at least one non empty contrast'), end;

if ~spm_sp('isspc',sX), sX = spm_sp('set',sX);  end;
for i=1:length(Fc1)
    if ~sf_IsFcon(Fc1(i)), error('Fc1(i) must be a contrast'), end
end
for i=1:length(Fc2)
    if ~sf_IsFcon(Fc2(i)), error('Fc2(i) must be a contrast'), end
end
varargout = { sf_Rortho(Fc1,sX,Fc2) };


case 'in'     %-  Fc1 is in list of  contrasts Fc2
%=======================================================================
% [iFc2 iFc1] = spm_FcUtil('In', Fc1, sX, Fc2)

% returns indice of Fc2 if "in", 0 otherwise 
% NB : If T- stat, the routine checks whether Fc.c is of
% size one. This is ensure if contrast is set 
% or manipulated (ortho ..) with spm_FcUtil
% note that the algorithmn works \emph{only because} Fc2(?).c 
% and Fc1.c are in space(X')

if nargin < 4, error('Insufficient arguments'), end
if nargout > 2, error('Too many output argument.'), end

Fc1 = varargin{2}; Fc2 = varargin{4}; sX = varargin{3};

L1  = length(Fc1);
if  ~L1, warning('no contrast given to in'); 
     if nargout == 2, varargout = {[] []}; 
     else varargout = {[]}; end;
     return; 
end
for i=1:L1
    if ~sf_IsFcon(Fc1(i)), error('Fc1(i) must be a contrast'), end
end
L2 = length(Fc2);
if ~L2, error('must have at least a contrast in Fc2'); end
for i=1:L2
    if ~sf_IsFcon(Fc2(i)), error('Fc2(i) must be F-contrast'), end
end
if ~spm_sp('isspc',sX), sX = spm_sp('set',sX);  end;


[idxFc2 idxFc1] =  sf_in(Fc1, sX, Fc2);
if isempty(idxFc2), idxFc2 = 0; end
if isempty(idxFc1), idxFc1 = 0; end

switch nargout 
    case {0,1}
        varargout = { idxFc2 };
    case 2
        varargout = { idxFc2 idxFc1 };
    otherwise 
        error('Too many or not enough output arguments');
end

case '~unique'     %-  Fc list unique 
%=======================================================================
% idx = spm_FcUtil('~unique', Fc, sX)

%- returns indices of redundant contrasts in Fc
%- such that Fc(idx) = [] makes Fc unique. 
%- if already unique returns [] 

if nargin ~= 3, error('Insufficient/too many arguments'), end
Fc = varargin{2}; sX = varargin{3}; 

%----------------------------
L1 = length(Fc);
if ~L1, warning('no contrast given '); varargout = {[]}; return; end
for i=1:L1
    if ~sf_IsFcon(Fc(i)), error('Fc(i) must be a contrast'), end
end
if ~spm_sp('isspc',sX), sX = spm_sp('set',sX);  end;
%----------------------------

varargout = { unique(sf_notunique(Fc, sX))};



case {'0|[]','[]|0'}     %-  Fc is null or empty 
%=======================================================================
% b = spm_FcUtil('0|[]', Fc, sX)

% returns 1 if F-contrast is empty or null; assumes the contrast is set.

if nargin ~= 3, error('Insufficient/too many arguments'), end
Fc = varargin{2}; sX = varargin{3}; 

%----------------------------
L1 = length(Fc);
if ~L1, warning('no contrast given to |_'); varargout = {[]}; return; end
for i=1:L1
    if ~sf_IsFcon(Fc(i)), error('Fc(i) must be a contrast'), end
end
if ~spm_sp('isspc',sX), sX = spm_sp('set',sX);  end;
%----------------------------

idx = [];
for i=1:L1
    if sf_isempty_X1o(Fc(i)) || sf_isnull(Fc(i),sX), idx = [idx i]; end
end
if isempty(idx) 
    varargout = {0};
else 
    varargout = {idx};
end

%=======================================================================

otherwise
%=======================================================================
error('Unknown action string in spm_FcUtil')


end; %---- switch lower(action),

%=======================================================================
%=======================================================================
% Sub Functions
%=======================================================================
%=======================================================================

%=======================================================================
% Fcon = spm_FcUtil('FconFields')

function Fc = sf_FconFields

Fc = struct(...
    'name',     '',...
    'STAT',     '',...
    'c',        [],...
    'X0',       struct('ukX0',[]),...
    'iX0',      [],...
    'X1o',      struct('ukX1o',[]),...
    'eidf',     [],...
    'Vcon',     [],...
    'Vspm',     []);
 
%=======================================================================
% used internally. Minimum contrast structure

function minFc = sf_MinFcFields

minFc = struct(...
    'name',     '',...
    'STAT',     '',...
    'c',        [],...
    'X0',       [],...
    'X1o',      []);

%=======================================================================
% yes_no = spm_FcUtil('IsFcon',Fc)

function b = sf_IsFcon(Fc)

%- check that minimum fields of a contrast are in Fc 
b    = 1;
minnames = fieldnames(sf_MinFcFields);
FCnames  = fieldnames(Fc);
for  str = minnames'
   b = b & any(strcmp(str,FCnames));
   if ~b, break, end
end

%=======================================================================
% used internally; To be set, a contrast structure should have
% either X1o or X0 non empty. X1o can be non empty because c
% is non empty.

function b = sf_IsSet(Fc)

b = ~sf_isempty_X0(Fc) | ~sf_isempty_X1o(Fc);

%=======================================================================
% used internally
%
function v = sf_ver(Fc)

if isstruct(Fc.X0), v = 2; else v = 1; end

%=======================================================================
% used internally

function b = sf_isempty_X1o(Fc)

if sf_ver(Fc) > 1, 
   b = isempty(Fc.X1o.ukX1o); 
   %- consistency check
   if b ~= isempty(Fc.c), 
      Fc.c, Fc.X1o.ukX1o, error('Contrast internally not consistent');
   end
else
   b = isempty(Fc.X1o);
   %- consistency check
   if b ~= isempty(Fc.c), 
      Fc.c, Fc.X1o, error('Contrast internally not consistent');
   end
end
%=======================================================================
% used internally

function b = sf_X1o(Fc,sX)

if sf_ver(Fc) > 1, 
   b = spm_sp('ox',sX)*Fc.X1o.ukX1o; 
else
   b = Fc.X1o;
end

%=======================================================================
% used internally

function b = sf_X0(Fc,sX)

if sf_ver(Fc) > 1, 
   b = spm_sp('ox',sX)*Fc.X0.ukX0; 
else
   b = Fc.X0;
end

%=======================================================================
% used internally

function b = sf_isempty_X0(Fc)

if sf_ver(Fc) > 1, 
   b = isempty(Fc.X0.ukX0); 
else
   b = isempty(Fc.X0);
end

%=======================================================================
% Hsqr = spm_Fcutil('Hsqr',Fc,sX)

function hsqr = sf_Hsqr(Fc,sX)
%
% Notations : H equiv to X1o, H = uk*a1, X = uk*ax, r = rk(sX), 
%             sX.X is (n,p), H is (n,q), a1 is (r,q), ax is (r,p)
%             oxa1 is an orthonormal basis for a1, oxa1 is (r,qo<=q)
% Algorithm : 
% v1 : Y'*H*(H'*H)-*H'*Y = b'*X'*H*(H'*H)-*H'*X*b
%                        = b'*X'*oxH*oxH'*X*b
% so hsrq is set to oxH'*X, a (q,n)*(n,p) op. + computation of oxH 
% v2 : X'*H*(H'*H)-*H'*X = ax'*uk'*uk*a1*(a1'*uk'*uk*a1)-*a1'*uk'*uk*ax
%                        = ax'*a1*(a1'*a1)-*a1'*ax
%                        = ax'*oxa1*oxa1'*ax
%                        
% so hsrq is set to oxa1'*ax : a (qo,r)*(r,p) operation!  -:)) 
% + computation of oxa1.

%-**** fprintf('v%d\n',sf_ver(Fc));
if sf_ver(Fc) > 1,    
   hsqr = spm_sp('ox',spm_sp('set',Fc.X1o.ukX1o))' * spm_sp('cukx',sX);
else
   hsqr = spm_sp('ox',spm_sp('set',Fc.X1o))'*spm_sp('x',sX);
end
%=======================================================================
% H = spm_FcUtil('H',Fc)

function H = sf_H(Fc,sX)
% 
% Notations : H equiv to X1o, H = uk*a1, X = uk*ax
% Algorithm : 
% v1 : Y'*H*(H'*H)-*H'*Y = b'*X'*H*(H'*H)-*H'*X*b
%                        = b'*c*(H'*H)-*c'*b
%                        = b'*c*(c'*(X'*X)-*c)-*c'*b
%- v1 : Note that pinv(Fc.X1o' * Fc.X1o) is not too bad
%- because dimensions are only (q,q). See sf_hsqr for notations.
%- Fc.c and Fc.X1o should match. This is ensure by using FcUtil.

%-**** fprintf('v%d\n',sf_ver(Fc));
if sf_ver(Fc) > 1, 
   hsqr = sf_Hsqr(Fc,sX);
   H    = hsqr' * hsqr;
else
   H = Fc.c * pinv(Fc.X1o' * Fc.X1o) * Fc.c';
   % H = {c*spm_sp('x-',spm_sp('Set',c'*spm_sp('xpx-',sX)*c) )*c'}
end

%=======================================================================
% Yc = spm_FcUtil('Yc',Fc,sX,b)

function Yc = sf_Yc(Fc,sX,b)

Yc =  sX.X*spm_sp('xpx-',sX)*sf_H(Fc,sX)*b;

%=======================================================================
% Y0 = spm_FcUtil('Y0',Fc,sX,b)
function Y0 = sf_Y0(Fc,sX,b)

Y0 = sX.X*(eye(spm_sp('size',sX,2)) - spm_sp('xpx-',sX)*sf_H(Fc,sX))*b;

%=======================================================================
% Fc = spm_FcUtil('|_',Fc1, sX, Fc2)

function Fc1o = sf_fcortho(Fc1, sX, Fc2)

%--- use the space facility to ensure the proper tolerance dealing...

c1_2  = Fc1.c - sf_H(Fc2,sX)*spm_sp('xpx-:',sX,Fc1.c);
Fc1o  = spm_FcUtil('Set',['(' Fc1.name ' |_ (' Fc2.name '))'], ...
            Fc1.STAT, 'c+',c1_2,sX);

%- In the large (scans) dimension :
%- c = sX.X'*spm_sp('r:',spm_sp('set',Fc2.X1o),Fc1.X1o);
%- Fc1o  = spm_FcUtil('Set',['(' Fc1.name ' |_ (' Fc2.name '))'], ...
%-                          Fc1.STAT, 'c',c,sX);

%=======================================================================
function b = sf_Rortho(Fc1,sX,Fc2)

if isempty(Fc2)
   if length(Fc1) <= 1, b = 0; 
   else
      c1 = cat(2,Fc1(:).c); 
      b = ~any(any( abs(triu(c1'*spm_sp('xpx-:',sX,c1), 1)) > sX.tol));
   end
else 
   c1 = cat(2,Fc1(:).c); c2 = cat(2,Fc2(:).c); 
   b = ~any(any( abs(c1'*spm_sp('xpx-:',sX,c2)) > sX.tol ));
end

%=======================================================================
% b = spm_FcUtil('0|[]', Fc, sX)

%- returns 1 if F-contrast is empty or null; assumes the contrast is set.
%- Assumes that if Fc.c contains only zeros, so does Fc.X1o. 
%- this is ensured if spm_FcUtil is used

function boul = sf_isnull(Fc,sX)
%
boul = ~any(any(spm_sp('oPp:',sX,Fc.c)));

%=======================================================================
% Fc = spm_FcUtil('Set',name, STAT, set_action, value, sX)

function boul = sf_is_T(sX,c)
%- assumes that the dimensions are OK
%- assumes c is not empty
%- Does NOT assumes that c is space of sX'
%- A rank of zero can be defined
%- if the rank == 1, checks whether same directions

boul = 1;
if ~spm_sp('isinspp',sX,c), c = spm_sp('oPp:',sX,c); end;
if rank(c) > 1 || any(any(c'*c < 0)), boul = 0; end;

%=======================================================================
function [idxFc2, idxFc1] =  sf_in(Fc1, sX, Fc2)

L2 = length(Fc2);
L1 = length(Fc1);

idxFc1 = []; idxFc2 = [];
for j=1:L1
   %- project Fc1(j).c if not estimable
   if ~spm_sp('isinspp',sX,Fc1(j).c), %- warning ? 
      c1  = spm_sp('oPp:',sX,Fc1(j).c); 
   else
      c1 = Fc1(j).c; 
   end
   sc1 = spm_sp('Set',c1); 
   S   = Fc1(j).STAT;

   boul = 0;
   for i =1:L2
      if Fc2(i).STAT == S
         %- the same statistics. else just go on to the next contrast
         boul = spm_sp('==',sc1,spm_sp('oPp',sX,Fc2(i).c));

         %- if they are the same space and T stat (same direction),
         %- then check wether they are in the same ORIENTATION
         %- works because size(X1o,2) == 1, else .* with (Fc1(j).c'*Fc2(i).c)
         if boul && S == 'T'
            atmp = sf_X1o(Fc1(j),sX);
            btmp = sf_X1o(Fc2(i),sX);
            boul = ~any(any( (atmp' * btmp) < 0 ));
         end
         %- note the indices
         if boul, idxFc1 = [idxFc1 j]; idxFc2 = [idxFc2 i]; end
      end
   end
end %- for j=1:L1

%=======================================================================
function idx = sf_notunique(Fc, sX)

%- works recursively ...
%- and use the fact that [] + i == []
%- quite long for large sets ... 

l = length(Fc);
if l == 1, idx = []; 
else
   idx = [ (1+sf_in(Fc(1),sX,Fc(2:l))) (1+sf_notunique(Fc(2:l), sX))];
end
