function varargout=spm_platform(varargin)
% Multifunction.  Platform specific routines for various spm functions
% FORMAT ans = spm_platform(arg)
% arg  - optional string argument, can be
%        - 'bigend'  - return whether this architecture is bigendian
%                      - Inf - is not IEEE floating point
%                      - 0   - is little end
%                      - 1   - big end
%        - 'filesys' - type of filesystem
%                      - 'unx' - UNIX
%                      - 'win' - DOS
%                      - 'mac' - Macintosh
%                      - 'vms' - VMS
%        - 'sepchar' - returns directory separator
%        - 'user'    - returns username
%        - 'tempdir' - returns name of temp directory
%
% Without input arguments, initialises the platform specific functions
% For calls with arguments, see body of function
%_______________________________________________________________________
% %W% Matthew Brett %E%

% Calls to the function work in one of two ways.  Calls which may be made
% often, and return values that are always the same over the spm session
% use the global variable spm_plat_vars as a static structure, initialised at
% first call to spm_platform, and return values from this structure.
% Other calls are standard, and process / return arguments at time of call

global spm_plat_vars

if isempty(spm_plat_vars)
	init_platform; end

if (nargin == 0)
	varargout={[]};
	return;
else
	Action = varargin{1};
end

switch lower(Action)
	
case 'bigend'
% return whether this architecture is big or little endian, with errors
varargout{1} = spm_plat_vars.bigend;
if ~finite(spm_plat_vars.bigend),
	if isnan(spm_plat_vars.bigend),
		error(['I don''t know if "' computer '" is big-endian.']);
	else,
		error(['I don''t think that "' computer '" uses IEEE floating point ops.']);
	end;
end

case 'filesys'
% return file system
varargout{1} = spm_plat_vars.filesys;

case 'sepchar'
% return file separator character
varargout{1} = spm_plat_vars.sepchar;

case 'user'
% return user string
varargout{1} = spm_plat_vars.user;

case 'shoes'
% warn about potential soling pitfalls
disp('Footwear fashion catastrophe');
error('Please retire to place of safety and destroy offensive walkwear');

case 'tempdir'
% return temporary directory
twd = getenv('SPMTMP');
if isempty(twd)
	switch spm_plat_vars.filesys
		case 'unx'
			twd = '/tmp';
		case 'win'
			twd = getenv(TEMP);
		otherwise
			error('Do not know how to set temp directory');
	end
end	
varargout{1} = twd;

otherwise
error('Illegal Action string')

end

%====================================================
% subfunctions
%====================================================

function init_platform
% Initialise variables on basis of architecture
global spm_plat_vars

% this and following arrays in matching order
computers = str2mat('PCWIN','MAC2','SUN4','SOL2','HP700','SGI','SGI64','IBM_RS',...
		'ALPHA','AXP_VMSG','AXP_VMSIEEE','LNX86','VAX_VMSG','VAX_VMSD');
% file systems
filesyses = str2mat('win','mac','unx','unx','vms','unx','unx','unx',...
			'unx','vms','vms','unx','vms','vms');
% endian codes
% NaN is don't know, Inf = not IEEE floating point, 0 is little end, 1 big end
endians = [0 1 1 1 1 1 1 1 0 Inf 0 0 Inf Inf];

% index into arrays
c = computer;
ci = NaN;
for i=1:size(computers,1),
	if strcmp(c,deblank(computers(i,:))),
		ci = i;
		break
	end
end
if isnan(ci)
	error(['Do not recognise architecture ' c]), end

% flag for big-endian.
spm_plat_vars.bigend =endians(ci);

% Last check for absence of IEEE floating point maths
if ~isieee
	spm_plat_vars.bigend =Inf; end

% assigns filesystem
spm_plat_vars.filesys = filesyses(ci, :);
switch (spm_plat_vars.filesys)
	case 'unx'
		spm_plat_vars.sepchar = '/';
	case 'win'
		spm_plat_vars.sepchar = '\';
	case 'mac'
		spm_plat_vars.sepchar = ':';
	otherwise
		error(['Do not know filesystem ' spm_plat_vars.filesys])
end

% gets user name
User = [];
switch (spm_plat_vars.filesys)
	case 'unx'
		User = getenv('USER');
	case 'win'
		User = getenv('USER');
end
if isempty(User), User='user'; end
spm_plat_vars.user = User;

% may want to set fonts etc here


% end of initialisation
return
