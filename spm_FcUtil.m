function varargout = spm_FcUtil(varargin)
%-----------------------------------------------------------------------
% F contrast utilities
% FORMAT varargout = spm_FcUtil(action,varargin)
%=======================================================================
%
%_______________________________________________________________________
%
% spm_FcUtil is a multi-function function containing various utilities
% for F contrast construction and manipulation. In
% general, it accepts design matrices as plain matrices or as space
% structures setup by spm_sp.
% 
% The use of spm_FcUtil should help with robustness issues and
% maintainability of SPM. 
%
% Note that when space structures are passed as arguments is is
% assummed that their basic fields are filled in. See spm_sp for
% details of (design) space structures and their manipulation.
%
% ======================================================================
%
%                           ----------------
% FORMAT spm_FcUtil('IsFcon',Fc), 
%                           ----------------
% FORMAT spm_FcUtil('FconFields'))'
%                           ----------------
% FORMAT Fcon = spm_FcUtil('Set',name, STAT, action, value, sX)
%                           ----------------
% ======================================================================

%-Format arguments
%-----------------------------------------------------------------------
if nargin==0, error('do what? no arguments given...')
	else, action = varargin{1}; end


switch lower(action),

case 'fconfields'				%- fields of F contrast
%=======================================================================
% Fcon = spm_FcUtil('FconFields')
%
%

if nargout > 1, error('Too many output arguments FconFields'), end;
if nargin > 1, error('Too many input arguments FconFields'), end;

varargout = {struct(...
	'name',		'',...
	'STAT',		'',...
	'c',			[],...
	'X0',			[],...
	'iX0',		[],...
	'X1o',		[],...
	'eidf',		[],...
	'Vcon',		[],...
	'Vspm',		[]	)};


case 'set'				%- Create an F contrast
%=======================================================================
% Fc = spm_FcUtil('Set',name, STAT, action, value, sX)
%
%
%

%--- check # arguments...
if nargin<6, error('insufficient arguments'), end;
if nargout > 1, error('Too many output arguments Set'), end;

%--- check arguments...
if ~isstr(varargin{2}), error('~isstr(name)'), end;
if ~(varargin{3}=='F'|varargin{3}=='T'), 
	error('~(STAT==F|STAT==T)'), end;
if ~isstr(varargin{4}), error('~isstr(varargin{4})'), end;

sX = varargin{6};
if ~spm_sp('isspc',sX), sX = spm_sp('set',sX);	end;
if isempty(sX.X), error('Empty space X in Set'); end;

Fc = sf_FconFields;
%- use the name as a flag to insure that F-contrast has been 
%- properly created;
Fc.name = varargin{2};
if isempty(Fc.name), Fc.name = ' '; end;
Fc.STAT = varargin{3};



[sC sL] = spm_sp('size',sX);

switch varargin{4},
	case 'c'
		c = varargin{5};
		if isempty(c)
				Fc.c 		= [];
				Fc.X1o 	= [];
				Fc.X0		= sX.X;
				Fc.iX0	= [];
		elseif size(c,1) ~= sL, 
	      	error('not contrast dimension in Set'); 
		else 
			   Fc.c 		= c;
   			[Fc.X1o Fc.X0] = spm_SpUtil('c->Tsp',sX,c);
				Fc.iX0	= [];
		end;

	case 'X0'
		X0 = varargin{5};
		if isempty(X0), 
				Fc.c 		= spm_sp('xpx',sX); 
				Fc.X1o 	= sX.X;
				Fc.X0 	= [];	
   			Fc.iX0	= 0;
	   elseif size(X0,1) ~= sC, 
				error('dimension of X0 wrong in Set');
	   else 
				Fc.c 		= spm_SpUtil('X0->c',sX,X0);
				Fc.X0 	= X0;	
				Fc.X1o 	= spm_SpUtil('c->Tsp',sX,Fc.c);
   			Fc.iX0	= 0;
		end  

	case 'iX0'
		iX0 = varargin{5};
	   Fc.iX0 	= iX0;
	   Fc.X0 	= sX.X(:,iX0);
		if isempty(iX0), 
				Fc.c 		= spm_sp('xpx',sX); 
				Fc.X1o 	= sX.X;  
		else 	
				iX0 		= spm_SpUtil('iX0check',iX0,sL);
   			Fc.c		= spm_SpUtil('i0->c',sX,iX0);
   			Fc.X1o	= spm_SpUtil('c->Tsp',sX,Fc.c);
		end;

	otherwise 
	   error('wrong action in Set '); 

end;
varargout = {Fc};



case 'isfcon'				%- Is it an F contrast ?
%=======================================================================
% yes_no = spm_FcUtil('IsFcon',Fc)

if nargin~=2, error('too few/many input arguments - need 2'), end
if ~isstruct(varargin{2}), varargout={0}; 
else, varargout = {sf_IsFcon(varargin{2})};
end



case 'fconedf'				%- F contrast edf
%=======================================================================
% [edf_testsp edf_Xsp] = spm_FcUtil('FconEdf', Fc, sX [, V])

if nargin<3, error('Insufficient arguments'), end
if nargout >= 3, error('Too many output argument.'), end
Fc = varargin{2};
sX = varargin{3};
if nargin == 4, V = varargin{4}; V_flag =1; else V_flag = 0; end;

if ~sf_IsFcon(Fc), error('Fc must be Fcon'), end
if ~spm_sp('isspc',sX)
	sX = spm_sp('set',sX);	end;
if ~spm_FcUtil('Rcompatible',Fc,sX), ...
 	error('sX and Fc must be compatible'), end

if ~V_flag
	if nargout == 2
   	[trMV, trMVMV] = spm_SpUtil('trMV',Fc.X1o);
   	[trRV, trRVRV] = spm_SpUtil('trRV',sX);
		varargout = {trMV^2/trMVMV, trRV^2/trRVRV};
	else 
   	[trMV, trMVMV] = spm_SpUtil('trMV',Fc.X1o);
		varargout = {trMV^2/trMVMV};
	end;
else 
	if nargout == 2
   	[trMV, trMVMV] = spm_SpUtil('trMV',Fc.X1o,V);
   	[trRV, trRVRV] = spm_SpUtil('trRV',sX,V);
		varargout = {trMV^2/trMVMV, trRV^2/trRVRV};
	else 
   	[trMV, trMVMV] = spm_SpUtil('trMV',Fc.X1o,V);
		varargout = {trMV^2/trMVMV};
	end;
end

case 'rcompatible'				%- Rcompatible F and sX
%=======================================================================
% Y_N = spm_FcUtil('Rcompatible',Fcon, sX)

% !!!! TO IMPLEMENT

varargout = {1};



%=======================================================================
%=======================================================================
%		part that use F contrast
%=======================================================================
%=======================================================================


case 'hsqr'     %-Extra sum of squares matrix for beta's from contrast
%=======================================================================
% hsqr = spm_FcUtil('Hsqr',Fc, sX)

if nargin<3, error('Insufficient arguments'), end
if nargout>1, error('Too many output argument.'), end
Fc = varargin{2};
sX = varargin{3};

if ~sf_IsFcon(Fc), error('Fc must be F-contrast'), end
if isempty(Fc.name), error('Fcon must be set'); end; %-
if ~spm_sp('isspc',sX), sX = spm_sp('set',sX);	end;
if ~spm_FcUtil('Rcompatible',Fc,sX), ...
	error('sX and Fc must be compatible'), end;

if isempty(Fc.X1o) | isempty(Fc.c)
	if ~empty(Fc.X0)
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

if nargin<2, error('Insufficient arguments'), end
if nargout>1, error('Too many output argument.'), end
Fc = varargin{2};
sX = varargin{3};

if ~sf_IsFcon(Fc), error('Fc must be F-contrast'), end
if isempty(Fc.name), error('Fcon must be set'); end; %-
if ~spm_sp('isspc',sX), sX = spm_sp('set',sX);	end;
if ~spm_FcUtil('Rcompatible',Fc,sX), ...
	error('sX and Fc must be compatible'), end;

if isempty(Fc.X1o) | isempty(Fc.c)
	if ~isempty(Fc.X0)
		%- assumes that X0 is sX.X
		%- warning(' Empty X1o in spm_FcUtil(''H'',Fc,sX) ');
		varargout = { zeros(spm_sp('size',sX,2)) };
	else 	
		error(' Fc must be set ');
	end
else
	varargout = { sf_H(Fc) };  
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
%
function Fc = sf_FconFields

Fc = struct(...
	'name',		'',...
	'STAT',		'',...
	'c',			[],...
	'X0',			[],...
	'iX0',		[],...
	'X1o',		[],...
	'eidf',		[],...
	'Vcon',		[],...
	'Vspm',		[]	);

%=======================================================================
% edf  = spm_FcUtil('FconEdf', Fc, sX, V)
%
% function Fc = sf_FconEdf(Fc,sX,V)



% warning('not to be used')
% Hb = spm_SpUtil('H',sX,Fc);
% if ~V_flag
%   Vb = (spm_sp('x-',sX)*V)*spm_sp('x-',sX)';
%   [trMV, trMVMV] = spm_SpUtil('trMV',Hb,Vb);
% else 
%   [trMV, trMVMV] = spm_SpUtil('trMV',Hb,spm_sp('xpx-',sX));
% end;


%=======================================================================
% Fcon = spm_FcUtil('fillFcon',Fcin, sX)
%
function Fc = sf_fillFcon(Fcin, sX)

Fc = Fcin;

if ~isempty(Fc.c) %- construct X1o and X0;
   [Fc.X1o Fc.X0] = spm_SpUtil('c->Tsp',sX,Fc.c);
   Fc.iX0 = [];
elseif ~isempty(Fc.iX0) 
   Fc.X0 	= sX.X(:,Fc.iX0);
   Fc.c		= spm_SpUtil('i0->c',sX,Fc.iX0);
   Fc.X1o	= spm_SpUtil('c->Tsp',sX,Fc.c);
else %- X0 not empty
   Fc.c		= spm_SpUtil('X0->c',sX,Fc.X0);
   Fc.X1o	= spm_SpUtil('c->Tsp',sX,Fc.c);
   Fc.iX0	= 0;
end


%=======================================================================
% yes_no = spm_FcUtil('IsFcon',Fc)
%
function b = sf_IsFcon(Fc)

b	= 1;
fnames  	= fieldnames(Fc);
for str 	= fieldnames(sf_FconFields)'
	b = b & any(strcmp(str,fnames));
	if ~b, break, end
end

%=======================================================================
% hsqr = spm_FcUtil('Hsqr',Fc,sX)
%
function hsqr = sf_Hsqr(Fc,sX)

hsqr = spm_sp('ox',spm_sp('set',Fc.X1o))'*sX.X;

%=======================================================================
% H = spm_FcUtil('H',Fc)
%
function H = sf_H(Fc)

%- Before a more efficient and robust way is implemented,
%- get the pinv(X1o'*X1o) ... JB :to be looked into !!!!
%- Note that pinv(Fc.X1o' * Fc.X1o) is not too bad
%- because the dimension should be small.

H = Fc.c * pinv(Fc.X1o' * Fc.X1o) * Fc.c';
% varargout = { c*spm_sp('x-',spm_sp('Set',c'*spm_sp('xpx-',sX)*c) )*c' }
