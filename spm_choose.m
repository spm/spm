function spm_choose(varargin)
% SPM: Version chooser function
% FORMAT spm_choose(SPMverCode,varargin)
% SPMverCode - [optional] short string identifier for SPM version to run
%_______________________________________________________________________
%  ___  ____  __  __
% / __)(  _ \(  \/  )  Statistical Parametric Mapping
% \__ \ )___/ )    (   The Wellcome Department of Cognitive Neurology
% (___/(__)  (_/\/\_)  
%
%_______________________________________________________________________
%
% spm_choose is a standalone SPM function which allows the user to
% choose between multiple versions of SPM. If no particular version is
% specified, an SPM chooser window is displayed with pushbuttons
% leading to the defined versions.
%
% Versions are defined in the body of the program via the parameter
% cell array SPMs. This has a row for each SPM version (up to 4), and
% has three columns defining the respective version descriptions, short
% names, and paths. Multiple short names & paths are allowed for a
% single version, provided they are embedded in cell arrays of strings
% themselves. See the parameters section in the code for further details.
%
% Once a version is chosen, the appropriate directory is prepended to
% the MatLab path, all functions cleared, and `spm` called. Note that
% since the SPM path is prepended, standard SPM functions will take
% precedence over functions in other directories on the path.
%
% Versions can be selected from the command line by specifying their
% short names as parameter. MatLab's command-function duality means
% that `spm_choose('99')` is equivalent to `spm_choose 99`.
%
% Systems administrators should code the appropriate definitions into
% the Parameters section. Note that options without a readable spm.m
% function are removed from the list of choices. If only one (valid)
% version is defined then this is run.
%
%-----------------------------------------------------------------------
%
% NOTE to SPM administrators on using this program:
%
%      * To use this program you must customise the hard-coded path and
%        version information in the SPMs structure in the main program
%        body (below).
%
%      * This program needs to be self contained.
%        It cannot rely on other SPM routines!
%_______________________________________________________________________
% %W% Andrew Holmes %E%

%-Parameters
%=======================================================================

%-SPM versions - up to 4 versions may be defined
%-----------------------------------------------------------------------
%         Description         | Short name(s) |   Path(s)
%-------+---------------------+---------------+-------------------------
SPMs ={	'SPMdevel',		'devel',	'/local/spm/spm_devel';...
	'SPM99  (FIL)',		'99',		'/local/spm/spm99_fil';...
	'SPM99  (public)',	'99p',		'/local/spm/spm99_pub';...
	'SPM99b (public)',	'99b',		'/local/spm/spm99b_pub';...
	'invalid',		'inv',		'/local/nowhere';
	'worse still'		'worse',	{'/tmp','/var/tmp'}};


%-Check parameters
%-----------------------------------------------------------------------
if isempty(SPMs), return, end


%-Version specified as parameter by SPMverCode - delete other options
%=======================================================================
if nargin>0
	iSPM  = [];
	SPMvs = cell(0,0);
	for i = 1:size(SPMs,1)
		if iscell(SPMs{i,2})
			iSPM  = [iSPM, i*ones(1,length(SPMs{i,2}(:)))];
			SPMvs = {SPMvs{:}, SPMs{i,2}{:}};
		else
			iSPM  = [iSPM, i];
			SPMvs = {SPMvs{:}, SPMs{i,2}};
		end
	end
	i = min(find(strcmp(varargin{1},SPMvs)));
	if ~isempty(i)
		SPMs = SPMs(iSPM(i),:);
	else
		error('Unrecognised SPMverCode')
	end
end


%-Check path validity: prune options without a readable spm.m,
% Make all path subcells into cell arrays themselves.
%=======================================================================
ok = zeros(size(SPMs,1),1);
for i = 1:size(SPMs,1)
	if iscell(SPMs{i,3})
		for j=1:length(SPMs{i,3})
			ok(i) = ok(i) | exist([SPMs{i,3}{j},'/spm.m'],'file')==2;
		end
	else
		ok(i) = exist([SPMs{i,3},'/spm.m'],'file')==2;
		SPMs(i,3)={SPMs(i,3)};
	end
end
SPMs(find(~ok),:)=[];

if isempty(SPMs)
	fprintf(['%s\nNo valid paths to choose an SPM version from!\n',...
		'\t- set version names and paths correctly in spm.m\n',...
		'\t  or set your MATLABPATH appropriately\n\n'],7)
	return
end


%-If only one version (remaining) specified then launch it and return
%=======================================================================
if size(SPMs,1)==1
	addpath(SPMs{1,3}{:});
	eval('spm(varargin{2:end})')
	return
end


%-Open chooser window
%=======================================================================
if size(SPMs,1) > 4, error('Not enough space for >4 buttons!'), end
S = get(0,'ScreenSize');

F = figure('IntegerHandle','off',...
	'Name','',...
	'NumberTitle','off',...
	'Tag','SPMchooser',...
	'Position',[S(3)/2-200,S(4)/2-140,300,280],...
	'Resize','off',...
	'Pointer','Watch',...
	'Color',[1 1 1]*.8,...
	'MenuBar','none',...
	'DefaultUicontrolFontSize',12,...
	'HandleVisibility','off',...
	'Visible','off');

%-Fontnames for this platform
%-----------------------------------------------------------------------
switch computer
case {'SUN4','SOL2','HP700','SGI','SGI64','IBM_RS','ALPHA','LNX86'}
	PF = struct('helvetica','Helvetica',...			%-UNIX
		'times','Times','courier','Courier');
case {'PCWIN'}
	PF = struct('helvetica','Arial Narrow',...		%-PCWIN
		'times','Times New Roman','courier','Courier New');
otherwise
	error([mfilename,' not setup for platform: ',computer])
end

%-Frames and text
%-----------------------------------------------------------------------
hA = axes('Parent',F,'Position',[0 0 100/300 280/280],'Visible','Off');
text(0.5,0.5,'SPM',...
	'Parent',hA,...
	'FontName',PF.times,'FontSize',96,...
	'FontAngle','Italic','FontWeight','Bold',...
	'Rotation',90,...
	'VerticalAlignment','Middle','HorizontalAlignment','center',...
	'Color',[1 1 1]*.6)

text(0.3,0.95,'Statistical Parametric Mapping',...
	'Parent',hA,...
	'FontName',PF.times,'FontSize',16,...
	'FontAngle','Italic','FontWeight','Bold',...
	'Color',[1 1 1]*.6)

uicontrol(F,'Style','Frame','Position',[110 020 180 230])

%-Buttons to launch SPM versions
%-----------------------------------------------------------------------
uicontrol(F,'Style','Text',...
	'String','Select version...',...
	'FontName',PF.times,'FontSize',14,...
	'HorizontalAlignment','Center',...
	'Position',[120 210 160 030],...
	'ForegroundColor','k')
for i = 1:min(size(SPMs,1),4)
	tmp = ['''',SPMs{i,3}{1},''''];
	for j=2:length(SPMs{i,3})
		tmp = [tmp,',''',SPMs{i,3}{j},''''];
	end
	uicontrol(F,'String',SPMs{i,1},...
		'Position',[120 180-(i-1)*35 160 030],...
		'CallBack',[...
			'addpath(',tmp,'),',...
			'delete(gcbf),',...
			'clear functions,',...
			'spm'],...
		'Interruptible','on',...
		'ForegroundColor',[0 1 1])
end

uicontrol(F,'String','Quit',...
	'Position',[120 030 160 030],...
	'CallBack','close(gcbf)',...
	'ForegroundColor','r');

set(F,'Pointer','Arrow','Visible','on')
