function spm_choose(SPMverCode)
% SPM: Version chooser function
% FORMAT spm_choose(SPMverCode)
% SPMverCode - [optional] short string identifier for SPM version to run
%_______________________________________________________________________
%  ___  ____  __  __
% / __)(  _ \(  \/  )  Statistical Parametric Mapping
% \__ \ )___/ )    (   The Wellcome Department of Cognitive Neurology
% (___/(__)  (_/\/\_)  University College London
%
%_______________________________________________________________________
%
% spm_choose is a standalone SPM function which allows the user to
% choose between multiple versions of SPM. If no particular version is
% specified, an SPM chooser window is displayed with pushbuttons
% leading to the defined versions.
%
% Versions are defined in terms of a triple of string matrices, holding
% the version name (the button label), a short name (for command line
% use), and the directory holding that version. These hardcoded into
% the program, in the "Parameters" section.
%
% Once a version is chosen, the appropriate directory is prepended to
% the MatLab path, all functions cleared, and the `spm` called.
%
% Versions can be selected from the command line by specifying their
% short names as parameter. MatLab's command-function duality means
% that `spm_choose('96')` is equivalent to `spm_choose 96`.
%
% Systems administrators should code the appropriate definitions into
% the Parameters section. Note that unaccessible directories are
% removed from the list of choices. If only one (valid) version is
% defined then this is run.
%
%_______________________________________________________________________
% %W% Andrew Holmes %E%

%-Parameters - up to 4 versions may be defined
%=======================================================================

%-Long names of the SPM versions
SPMverNames = str2mat(...
		'SPM''94  (MRC-CU)',...
		'SPM''95  (FIL)',...
		'SPM''96  (FIL)',...
		'SPM''97d (devel)');

%-Shortcut keys for each version.
SPMverCodes = ...
	  str2mat('94','95','96','97d');

%-Paths for each version
SPMverPaths = str2mat(...
		'/local/spm/spm94_mrccu',...
		'/local/spm/spm95_fil',...
		'/local/spm/spm96_fil',...
		'/local/spm/spm_devel');



%-Check parameters
%=======================================================================
if isempty(SPMverNames), return, end
if size(SPMverNames,1)~=size(SPMverPaths,1) | ...
		size(SPMverNames,1)~=size(SPMverCodes,1)
	error('Versions - ShortCodes - Directories mismatch'), end

%-Check path validity, prune directories without a readable spm.m
%-----------------------------------------------------------------------
tmp = zeros(size(SPMverPaths,1),1);
for SPMver = 1:length(tmp)
	tmp(SPMver) = exist([deblank(SPMverPaths(SPMver,:)),'/spm.m'])==2;
end
SPMverNames = SPMverNames(tmp,:);
SPMverCodes = SPMverCodes(tmp,:);
SPMverPaths = SPMverPaths(tmp,:);
if isempty(SPMverNames)
	fprintf(['%s\nNo valid paths to choose an SPM version from!\n',...
		'- set version names and paths correctly in spm.m\n',...
		'  or set your MATLABPATH appropriately\n\n'],7)
	return
end

%-Version specified as parameter by SPMverCode
%-----------------------------------------------------------------------
if nargin>0
	n    = size(SPMverCodes,1);
	tmp  = str2mat(SPMverCodes,SPMverCode);
	SPMver = ...
		find(all(tmp(1:n,:)'==...
		setstr(ones(size(SPMverCodes,1),1)*tmp(n+1,:))'));
	if SPMver
		SPMverNames = deblank(SPMverNames(SPMver,:));
		SPMverPaths = deblank(SPMverPaths(SPMver,:));
	end
end

%-If only one version specified then launch it and return
%-----------------------------------------------------------------------
if size(SPMverNames,1)==1
	if exist([SPMverPaths,'/spm.m'])~=2
		error(['Invalid path for "',SPMverNames,'"']), end
	path(SPMverPaths,path)
	eval('spm')
	return
end

%-Open chooser window
%=======================================================================
if size(SPMverNames,1) > 4, error('Not enough space for >4 buttons!'), end
S = get(0,'ScreenSize');
close(findobj(get(0,'Children'),'Tag','SPMchooser'))
F = figure('Color',[1 1 1]*.8,...
	'Name','',...
	'NumberTitle','off',...
	'Position',[S(3)/2-200,S(4)/2-140,300,280],...
	'Resize','off',...
	'Tag','SPMchooser',...
	'Pointer','Watch',...
	'Visible','off');

%-Frames and text
%-----------------------------------------------------------------------
axes('Position',[0 0 80/300 280/280],'Visible','Off')
text(0.5,0.5,'SPM',...
	'FontName','Times','FontSize',96,...
	'Rotation',90,...
	'VerticalAlignment','middle','HorizontalAlignment','center',...
	'Color',[1 1 1]*.6);

text(0.3,0.95,'Statistical Parametric Mapping',...
	'FontName','Times','FontSize',16,'FontAngle','Italic',...
	'FontWeight','Bold',...
	'Color',[1 1 1]*.6);

uicontrol(F,'Style','Frame','Position',[110 020 180 230]);

%-Buttons to launch SPM versions
%-----------------------------------------------------------------------
uicontrol(F,'Style','Text',...
	'String','Select version...',...
	'HorizontalAlignment','Center',...
	'Position',[120 210 160 030],...
	'ForegroundColor','k');
for i=1:size(SPMverNames,1)
	uicontrol(F,'String',deblank(SPMverNames(i,:)),...
		'Position',[120 175-(i-1)*35 160 025],...
		'CallBack',[...
			'path(get(gco,''UserData''),path),',...
			'close(gcf),',...
			'clear functions,',...
			'eval(''spm'')'],...
		'UserData',deblank(SPMverPaths(i,:)),...
		'Interruptible','yes',...
		'ForegroundColor',[0 1 1]);
end

uicontrol(F,'String','Quit',...
	'Position',[120 030 160 030],...
	'CallBack','close(gcf)',...
	'ForegroundColor','r');

set(F,'Pointer','Arrow','Visible','on')
