function varargout = spm_FcUtil(varargin)
% Contrast utilities
% FORMAT varargout = spm_FcUtil(action,varargin)
%_______________________________________________________________________
%
% spm_FcUtil is a multi-function function containing various utilities
% for contrast construction and manipulation. In general, it accepts
% design matrices as plain matrices or as space structures setup by
% spm_sp.
% 
% The use of spm_FcUtil should help with robustness issues and
% maintainability of SPM.  % Note that when space structures are passed
% as arguments is is assummed that their basic fields are filled in.
% See spm_sp for details of (design) space structures and their
% manipulation.
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
%_______________________________________________________________________
% %W% Jean-Baptiste Poline %E%

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
	'c',		[],...
	'X0',		[],...
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
if Fc.STAT=='T'  &  ~(action == 'c' | action == 'c+') 
	   warning('enter T stat with contrast');
end;


[sC sL] = spm_sp('size',sX);

switch varargin{4},
	case {'c','c+'}
		c = varargin{5};
		if (varargin{4} == 'c+')  
		   if ~spm_sp('isinspp',sX,c), c = spm_sp('oPp',sX,c); end;
		end;
		if isempty(c)
			Fc.c 	= [];
			Fc.X1o 	= [];
			Fc.X0	= sX.X;
			Fc.iX0	= [];
		elseif size(c,1) ~= sL, 
	      		error('not contrast dimension in Set'); 
		else 
			Fc.c 	= c;
   			[Fc.X1o Fc.X0] = spm_SpUtil('c->Tsp',sX,c);
			Fc.iX0	= [];

			%-- check the rank if T stat
			if Fc.STAT=='T'
			   if abs(rank(Fc.X1o) - 1) > eps,
			      error('Set tries to define a t that is a F');
			   end
			end
		end;

	case 'X0'
		X0 = varargin{5};
		if isempty(X0), 
			Fc.c 	= spm_sp('xpx',sX); 
			Fc.X1o 	= sX.X;
			Fc.X0 	= [];	
   			Fc.iX0	= 0;
	   	elseif size(X0,1) ~= sC, 
			error('dimension of X0 wrong in Set');
	   	else 
			Fc.c 	= spm_SpUtil('X0->c',sX,X0);
			Fc.X0 	= X0;	
			Fc.X1o 	= spm_SpUtil('c->Tsp',sX,Fc.c);
   			Fc.iX0	= 0;
		end  

	case 'iX0'
		iX0 	= varargin{5};
		iX0 	= spm_SpUtil('iX0check',iX0,sL);
	   	Fc.iX0 	= iX0;
	   	Fc.X0 	= sX.X(:,iX0);
		if isempty(iX0), 
			Fc.c 	= spm_sp('xpx',sX); 
			Fc.X1o 	= sX.X;  
		else 			
   			Fc.c	= spm_SpUtil('i0->c',sX,iX0);
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
% [edf_tsp edf_Xsp] = spm_FcUtil('FconEdf', Fc, sX [, V])

if nargin<3, error('Insufficient arguments'), end
if nargout >= 3, error('Too many output argument.'), end
Fc = varargin{2};
sX = varargin{3};
if nargin == 4, V = varargin{4}; else V = []; end;

if ~sf_IsFcon(Fc), error('Fc must be Fcon'), end
if ~spm_sp('isspc',sX)
	sX = spm_sp('set',sX);	end;
if ~spm_FcUtil('Rcompatible',Fc,sX), ...
 	error('sX and Fc must be compatible'), end

 
[trMV, trMVMV] = spm_SpUtil('trMV',Fc.X1o,V);
	
if ~trMVMV, edf_tsp = 0; warning('edf_tsp = 0'), 
else,  edf_tsp = trMV^2/trMVMV; end;	

if nargout == 2

   [trRV, trRVRV] = spm_SpUtil('trRV',sX,V);
   if ~trRVRV, edf_Xsp = 0, warning('edf_Xsp = 0'),
   else,  edf_Xsp = trRV^2/trRVRV; end;

   varargout = {edf_tsp, edf_Xsp};

else 	
   varargout = {edf_tsp};
end;



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
	if ~isempty(Fc.X0)
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
%
% Empty and zeros dealing : 
% This routine never returns an empty matrix. 
% If isempty(Fc.X1o) | isempty(Fc.c) it explicitly 
% returns a zeros projection matrix.

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




case 'yc'     %- Fitted data corrected for confounds defined by Fc 
%=======================================================================
% Yc = spm_FcUtil('Yc',Fc, sX, b)

if nargin < 4, error('Insufficient arguments'), end
if nargout > 1, error('Too many output argument.'), end

Fc = varargin{2}; sX = varargin{3}; b = varargin{4};

if ~sf_IsFcon(Fc), error('Fc must be F-contrast'), end
if isempty(Fc.name), error('Fcon must be set'); end; 
if ~spm_sp('isspc',sX), sX = spm_sp('set',sX);	end;
if ~spm_FcUtil('Rcompatible',Fc,sX), ...
	error('sX and Fc must be compatible'), end;
if spm_sp('size',sX,2) ~= size(b,1), 
	error('sX and b must be compatible'), end;

if isempty(Fc.X1o) | isempty(Fc.c)
	if ~isempty(Fc.X0)
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
if isempty(Fc.name), error('Fcon must be set'); end; 
if ~spm_sp('isspc',sX), sX = spm_sp('set',sX);	end;
if ~spm_FcUtil('Rcompatible',Fc,sX), ...
	error('sX and Fc must be compatible'), end;
if spm_sp('size',sX,2) ~= size(b,1), 
	error('sX and b must be compatible'), end;

if isempty(Fc.X1o) | isempty(Fc.c)
	if ~isempty(Fc.X0)
	   %- if space of interest empty or null, returns zeros !
	   varargout = { sX.X*b };
	else 	
	   error(' Fc must be set ');
	end
else
	varargout = { sf_Y0(Fc,sX,b) };  
end





case {'|_','fcortho'}     %-  Fc orthogonalisation 
%=======================================================================
% Fc = spm_FcUtil('|_',Fc1, Fc2, sX)
%
% returns Fc1 orthogonolised wrt Fc2 

if nargin < 4, error('Insufficient arguments'), end
if nargout > 1, error('Too many output argument.'), end

Fc1 = varargin{2}; Fc2 = varargin{3}; sX = varargin{4};

if ~sf_IsFcon(Fc1), error('Fc1 must be F-contrast'), end
L = length(Fc2);
for i=1:L
    if ~sf_IsFcon(Fc2(i)), error('Fc2(i) must be F-contrast'), end
end
if ~spm_sp('isspc',sX), sX = spm_sp('set',sX);	end;


if sf_isempty(Fc1) | sf_isnull(Fc1)
	%- Fc1 is an empty or null F-contrast : orthognonal to anything
	varargout = { Fc1 };
else
	%- create an F-contrast for all the Fc2

	str  = Fc2(1).name;
	for i=2:L str = [str ' ' Fc2(i).name]; end;

	Fc2  = spm_FcUtil('Set',str,Fc1.STAT,'c',cat(2,Fc2(:).c),sX);

	if sf_isempty(Fc2) | sf_isnull(Fc2)
		varargout = { Fc1 };
	else 
		varargout = { sf_fcortho(Fc1, Fc2, sX) };  
	end
end





case 'in'     %-  Fc1 is in list of  contrasts Fc2
%=======================================================================
% iFc2 = spm_FcUtil('In',Fc1, Fc2, sX)
%
% returns indice of Fc2 if in, 0 otherwise 

%----------------------------
if nargin < 4, error('Insufficient arguments'), end
if nargout > 1, error('Too many output argument.'), end

Fc1 = varargin{2}; Fc2 = varargin{3}; sX = varargin{4};

if ~sf_IsFcon(Fc1), error('Fc1 must be F-contrast'), end
L = length(Fc2);
for i=1:L
    if ~sf_IsFcon(Fc2(i)), error('Fc2(i) must be F-contrast'), end
end
if ~spm_sp('isspc',sX), sX = spm_sp('set',sX);	end;
%----------------------------

%- project Fc1.c if not estimable
if ~spm_sp('isinspp',sX,Fc1.c), c1 = spm_sp('oPp',sX,Fc1.c);
else, c1 = Fc1.c; end

sc1 = spm_sp('Set',c1);
S   = Fc1.STAT;

boul = 0; i = 1;
while ~boul & i <= L,
   if Fc2(i).STAT==S
	boul = spm_sp('==',sc1,spm_sp('oPp',sX,Fc2(i).c));
	%- if they are the same space and T stat (same direction),
	%- then check wether they are in the same orientation
	if boul & S == 'T'
	   %- assumes size(Fc1.c,2) == 1
	   if size(Fc1.c,2) == 1 & size(Fc2(i).c,2) == 1,
	   	boul = ((sX.X*Fc1.c)' * sX.X*Fc2(i).c >0 );
	   else 
		error('not implemented');
	   end;

	   %if size(Fc1.X1o,1) ~= 1, X1o1 = orth(Fc1.X1o); 
	   %else, X1o1=Fc1.X1o; end
	   %if size(Fc2(i).X1o,1) ~= 1, X1o2 = orth(Fc2(i).X1o);
	   %else, X1o2=Fc2(i).X1o; end
	   %boul = ( X1o1'*X1o2 >= 0)  % same orientation
	end
   end;
   i = i+1;
end
if boul, varargout = { i-1 }; else varargout = { 0 }; end;





case {'0|[]','[]|0','isemptyornull'}     %-  Fc orthogonalisation 
%=======================================================================
% Fc = spm_FcUtil('0|[]',Fc)
%
% returns 1 if F-contrast is empty or null; assumes the contrast is set.

if nargin < 2, error('Insufficient arguments'), end
Fc = varargin{2}; 

if ~sf_IsFcon(Fc), error('Fc must be F-contrast'), end

varargout = { sf_isempty(Fc) | sf_isnull(Fc) };

 


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
	'c',		[],...
	'X0',		[],...
	'iX0',		[],...
	'X1o',		[],...
	'eidf',		[],...
	'Vcon',		[],...
	'Vspm',		[]	);

%=======================================================================
% edf  = spm_FcUtil('FconEdf', Fc, sX, V)
%
% function Fc = sf_FconEdf(Fc,sX,V)



% To investigate
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
fnames  = fieldnames(Fc);
for str = fieldnames(sf_FconFields)'
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
% varargout = {c*spm_sp('x-',spm_sp('Set',c'*spm_sp('xpx-',sX)*c) )*c'}

%=======================================================================
% Yc = spm_FcUtil('Yc',Fc,sX,b)
%
function Yc = sf_Yc(Fc,sX,b)

Yc =  sX.X*spm_sp('xpx-',sX)*sf_H(Fc)*b;

%=======================================================================
% Y0 = spm_FcUtil('Y0',Fc,sX,b)
%
function Y0 = sf_Y0(Fc,sX,b)

Y0 =  sX.X*( eye(spm_sp('size',sX,2)) - spm_sp('xpx-',sX)*sf_H(Fc) )*b;

%=======================================================================
% Fc = spm_FcUtil('|_',Fc1, Fc2, sX)
%
function Fc1o = sf_fcortho(Fc1, Fc2, sX)

%--- use the space facility to ensure the proper tolerance dealing...
sC2   = spm_sp('set',Fc2.c);
c1o   = spm_sp('xpx',sX)*spm_sp('res',sC2,Fc1.c);
%--  to put in spm_sp ?
%--  NB : usually, the tol of sX will be much greater than the tol of sC2
c1o( abs(c1o) < max(sX.tol,sC2.tol) ) = 0;
Fc1o  = spm_FcUtil('Set',['(' Fc1.name ' |_ (' Fc2.name '))'], ...
		    Fc1.STAT, 'c',c1o,sX);


%=======================================================================
% Fc = spm_FcUtil('0|[]',Fc)
%
% returns 1 if F-contrast is empty or null; assumes the contrast is set.
%

function boul = sf_isnull(Fc)
%
%- Assumes that if Fc.c contains only zeros, so does Fc.X1o. 
%- this is ensured if spm_FcUtil is used
boul = ~any(any(Fc.c));

function boul = sf_isempty(Fc)
%
%- Assumes that if Fc.c is empty, so is Fc.X1o. 
%- this is ensured if spm_FcUtil is used
boul = isempty(Fc.c);

