
% Display and analysis of regional effects
% FORMAT spm_results_ui
%_______________________________________________________________________
%
% spm_results_ui prompts for details of a SPM{Z} or SPM{F} that is then
% displayed in the results window.  By moving (dragging) the pointer on
% the maximum intensity projection one can select a particular
% voxel (or planes passing though that voxel) for further display and
% analysis:
%
% 1) plot the activities at the point selected
% 2) render the SPM{Z/F} on transverse slices at the level selected
% 3) render the SPM{Z/F} through 3 orthogonal sections selected
% 4) list and characterize all the maxima in a selected region
% 5) write out the filtered image as an Analyze volume
%
% The SPM{Z/F} is thresholded at a specified threshold and displayed
% on top of a background image.  This background image can be arbitrary
% but should have ORIGIN correctly specified in its header.  ORIGIN
% is the [i j k] voxel corresponding to [0 0 0] mm in XYZ [typically
% the anterior commissure in the space defined by the atlas of Talairach
% and Tournoux (1988)].
%
% The specification of the contrasts to use and the height and size
% thresholds are the same as that described in spm_projections_ui.
%_______________________________________________________________________
% %W% Karl Friston, Andrew Holmes %E%

%-Programers notes
%-----------------------------------------------------------------------
% spm_results_ui is not a function so that data and results can be
% held in the base workspace. Thus the spm_*.m routines invoked are
% also M-files. This scenario enables the user easy access to the data
% structures.
%
% The interactive MIP is handled by spm_mip_ui
%_______________________________________________________________________

% Initialise 
%-----------------------------------------------------------------------
clear
WS     = spm('GetWinScale');
Fgraph = spm_figure('FindWin','Graphics');
Finter = spm_figure('FindWin','Interactive');
spm_clf(Fgraph)
spm_clf(Finter)
set(Finter,'Name','SPM results')
global CWD

% Which SPM
%-----------------------------------------------------------------------
SPMZ     = spm_input('which SPM',1,'b','SPM{Z}|SPM{F}',[1 0]);
SPMF     = ~SPMZ;

% Get thresholded data, thresholds and parameters
%-----------------------------------------------------------------------
if SPMZ

	[t,XYZ,QQ,U,k,s,w] = spm_projections_ui('Results');

elseif SPMF

	[t,XYZ,QQ,U,k,s,w] = spm_projectionsF_ui('Results');

end

% proceed if there are voxels
%-----------------------------------------------------------------------
if ~length(t); return; end
set(Finter,'Name','SPM results','Pointer','Watch')

%-Load SPM.mat file from appropriate directory
% (CWD, set by spm_projections*_ui)
%-----------------------------------------------------------------------
K        = [];
FLIP     = 0;
load([CWD,'/SPM'])
S        = s;
W        = w;

if SPMF, df = Fdf; end

% Make description strings
%-----------------------------------------------------------------------
if SPMZ
	tmp     = 1 - spm_Ncdf(U);
elseif SPMF
	tmp     = 1 - spm_Fcdf(U,df);
end

Descrip = spm_str_manip(CWD,'a30');
Descrip = str2mat(Descrip,...
	sprintf('Height threshold {u} = %0.2f, p = %0.6f',U,tmp));
Descrip = str2mat(Descrip,...
	sprintf('Extent threshold {k} = %i voxels',k));

% Display SPM{Z/F} & setup buttons
%=======================================================================
spm_mip_ui(t,XYZ,V,Descrip)

% Objects with calls to spm_*.m routines
%-----------------------------------------------------------------------
uicontrol(Fgraph,'Style','Frame','Position',[20 470 100 185].*WS);

tmp = 'L = spm_mip_ui(''GetCoords'')'';';

uicontrol(Fgraph,'Style','Text',...
	'String','Results',...
	'Position',[40 630 60 20].*WS,...
	'HorizontalAlignment','Center',...
	'ForegroundColor','w')
uicontrol(Fgraph,'Style','Pushbutton','String','plot',...
	'Callback',[tmp 'spm_graph;'],...
	'Position',[40 610 60 20].*WS,...
	'Interruptible','yes');
uicontrol(Fgraph,'Style','Pushbutton','String','maxima',...
	'Callback',[tmp 'spm_maxima;'],...
	'Position',[40 550 60 20].*WS,...
	'Interruptible','yes');
uicontrol(Fgraph,'Style','Pushbutton','String','slices',...
	'Callback',[tmp 'spm_transverse;'],...
	'Position',[40 580 60 20].*WS,...
	'Interruptible', 'yes');
uicontrol(Fgraph,'Style','Pushbutton','String','write',...
	'Callback','spm_write_filtered(t,XYZ,U,k,V,SPMZ);',...
	'Position',[40 520 60 20].*WS,...
	'Interruptible', 'yes');
if V(3) > 1
	uicontrol(Fgraph,'Style','Pushbutton','String','sections',... 
		'Callback',[tmp 'spm_sections;'],...
		'Position',[40 490 60 20].*WS,...
		'Interruptible','yes');
end


% Finished
%-----------------------------------------------------------------------
set(Finter,'Pointer','Arrow')

















