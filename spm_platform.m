function varargout=spm_platform(varargin)
% Platform specific configuration parameters for SPM
%
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
%        - 'rootlen' - returns number of chars in root directory name
%        - 'user'    - returns username
%        - 'tempdir' - returns name of temp directory
%
% FORMAT PlatFontNames = spm_platform('fonts')
% Returns structure with fields named after the generic (UNIX) fonts, the
% field containing the name of the platform specific font.
%
% FORMAT PlatFontName = spm_platform('font',GenFontName)
% Maps generic (UNIX) FontNames to platform specific FontNames
%
% FORMAT SPM_PLATFORM = spm_platform('init',comp)
% Initialises platform specific parameters in global SPM_PLATFORM
% (External gateway to init_platform(comp) subfunction)
% comp         - computer to use [defaults to MatLab's `computer`]
% SPM_PLATFORM - copy of global SPM_PLATFORM
%
% FORMAT spm_platform
% Initialises platform specific parameters in global SPM_PLATFORM
% (External gateway to init_platform(computer) subfunction)
%
% FORMAT spm_platform('clear')
% Clears global SPM_PLATFORM containing platform specific parameters 
%
%                           ----------------
% SUBFUNCTIONS:
%
% FORMAT init_platform(comp)
% Initialise platform specific parameters in global SPM_PLATFORM
% comp         - computer to use [defaults to MatLab's `computer`]
%
%-----------------------------------------------------------------------
%
% Since calls to spm_platform will be made frequently, most platform
% specific parameters are stored as a structure in the global variable
% SPM_PLATFORM. Subsequent calls use the information from this global
% variable, if it exists.
%
% Platform specific difinitions are contained in the data structures at
% the beginning of the init_platform subfunction at the end of this
% file.
%_______________________________________________________________________
% %W% Matthew Brett %E%


%-Initialise
%-----------------------------------------------------------------------
global SPM_PLATFORM
if isempty(SPM_PLATFORM), init_platform, end

if nargin==0, return, end


switch lower(varargin{1}), case 'init'                  %-(re)initialise
%=======================================================================
init_platform(varargin{2:end})
varargout = {SPM_PLATFORM};
   
case 'clear'                                       %-Clear SPM_PLATFORM
%=======================================================================
clear global SPM_PLATFORM

case 'bigend'                      %-Return endian for this architecture
%=======================================================================
varargout = {SPM_PLATFORM.bigend};
if ~finite(SPM_PLATFORM.bigend),
	if isnan(SPM_PLATFORM.bigend)
		error(['I don''t know if "',computer,'" is big-endian.'])
	else
		error(['I don''t think that "',computer,...
			'" uses IEEE floating point ops.'])
	end
end

case 'filesys'                                      %-Return file system
%=======================================================================
varargout = {SPM_PLATFORM.filesys};

case 'sepchar'                         %-Return file separator character
%=======================================================================
warning('use filesep instead (supported by MathWorks)')
varargout = {SPM_PLATFORM.sepchar};

case 'rootlen'           %-Return length in chars of root directory name 
%=======================================================================
varargout = {SPM_PLATFORM.rootlen};

case 'user'                                         %-Return user string
%=======================================================================
varargout = {SPM_PLATFORM.user};

case 'tempdir'                              %-Return temporary directory
%=======================================================================
twd = getenv('SPMTMP');
if isempty(twd)
    switch SPM_PLATFORM.filesys
    case 'unx'
        twd = '/tmp';
    case 'win'
        twd = getenv('TEMP');
        if isempty(twd)
            for tmp = {'c:\temp',[getenv('WINDIR'),'\Temp']}
                if mkdir('',tmp), twd=tmp; break, end
            end
            if isempty(twd)
                error('Could not find or create temporary directory')
            end
        end
    otherwise
        error('Do not know how to set temp directory');
    end
end 
varargout = {twd};


case {'font','fonts'}    %-Map default font names to platform font names
%=======================================================================
if nargin<2, varargout={SPM_PLATFORM.font}; return, end
switch lower(varargin{2})
case 'times'
	varargout = {SPM_PLATFORM.font.times};
case 'courier'
	varargout = {SPM_PLATFORM.font.courier};
case 'helvetica'
	varargout = {SPM_PLATFORM.font.helvetica};
case 'symbol'
	varargout = {SPM_PLATFORM.font.symbol};
otherwise
	warning(['Unknown font ',varargin{2},', using default'])
	varargout = {SPM_PLATFORM.font.helvetica};
end

otherwise                                        %-Unknown Action string
%=======================================================================
error('Unknown Action string')

%=======================================================================
end



%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================


function init_platform(comp)             %-Initialise platform variables
%=======================================================================
if nargin<1, comp=computer; end
global SPM_PLATFORM

%-Platform definitions
%-----------------------------------------------------------------------
PDefs = {	'PCWIN',	'win',	0;...
		'MAC2',		'mac',	1;...
		'SUN4',		'unx',	1;...
		'SOL2',		'unx',	1;...
		'HP700',	'vms',	1;...
		'SGI',		'unx',	1;...
		'SGI64',	'unx',	1;...
		'IBM_RS',	'unx',	1;...
		'ALPHA',	'unx',	0;...
		'AXP_VMSG',	'vms',	Inf;...
		'AXP_VMSIEEE',	'vms',	0;...
		'LNX86',	'unx',	0;...
		'VAX_VMSG',	'vms',	Inf;...
		'VAX_VMSD',	'vms',	Inf	};

PDefs = cell2struct(PDefs,{'computer','filesys','endian'},2);


%-Which computer?
%-----------------------------------------------------------------------
ci = find(strcmp({PDefs.computer},comp));
if isempty(ci), error([comp,' not supported architecture for SPM']), end


%-Set bigend
%-----------------------------------------------------------------------
SPM_PLATFORM.bigend = PDefs(ci).endian;
%-Last check for absence of IEEE floating point maths
if ~isieee, SPM_PLATFORM.bigend = Inf; end	%-Last check for IEEE math


%-Set filesys
%-----------------------------------------------------------------------
SPM_PLATFORM.filesys = PDefs(ci).filesys;


%-Set filesystem dependent stuff
%-----------------------------------------------------------------------
%-File separators character
%-Length of root directory strings
%-User name finding
%-(mouse button labels?)
switch (SPM_PLATFORM.filesys)
case 'unx'
	SPM_PLATFORM.sepchar = '/';
	SPM_PLATFORM.rootlen = 1;
	SPM_PLATFORM.user    = getenv('USER');
case 'win'
	SPM_PLATFORM.sepchar = '\';
	SPM_PLATFORM.rootlen = 3;
	SPM_PLATFORM.user    = getenv('USERNAME');
	if isempty(SPM_PLATFORM.user)
		SPM_PLATFORM.user = spm_win32utils('username'); end
case 'mac'
	SPM_PLATFORM.sepchar = ':';
	SPM_PLATFORM.rootlen = 1;			%-** Not sure!?
	SPM_PLATFORM.user    = '';			%-** Dunno!
otherwise
	error(['Don''t know filesystem ',SPM_PLATFORM.filesys])
end

%-Fonts
%-----------------------------------------------------------------------
switch comp
case {'SUN4','SOL2','HP700','SGI','SGI64','IBM_RS','ALPHA','LNX86'}
	SPM_PLATFORM.font.helvetica = 'Helvetica';
	SPM_PLATFORM.font.times     = 'Times';
	SPM_PLATFORM.font.courier   = 'Courier';
	SPM_PLATFORM.font.symbol    = 'Symbol';
case {'PCWIN'}
	SPM_PLATFORM.font.helvetica = 'Arial Narrow';
	SPM_PLATFORM.font.times     = 'Times New Roman';
	SPM_PLATFORM.font.courier   = 'Courier New';
	SPM_PLATFORM.font.symbol    = 'Symbol';
end
