
% Display and analysis of regional effects
% FORMAT spm_results
%_______________________________________________________________________
%
% The SPM results section is for the interactive exploration and
% characterisation of the results of a statistical analysis.
% 
% The user is prompted to select a SPM{Z} or SPM{F}, which is
% thresholded at user specified levels. The specification of the
% contrasts to use and the height and size thresholds are the same as
% that described in spm_projections_ui.m (for SPM{Z}) &
% spm_projectionsF_ui.m (fpr SPM{F}). The resulting SPM is then
% displayed in the graphics window as a maximum intensity projection,
% alongside the design matrix.
%
% The cursors in the MIP can be moved (dragged) to select a particular
% voxel. The three mouse buttons give different drag and drop behavoir:
% Button 1 - point & drop; Button 2 - "dynamic" drag & drop with
% co-ordinate & SPM value updating; Button 3 - "magnetic" drag & drop,
% where the cursor jumps to the nearest suprathreshold voxel in the
% MIP, and shows the value there. See spm_mip_ui.m, the MIP GUI handling
% function.
%
% The current voxel specifies the voxel, suprathreshold cluster, or
% orthogonal planes (planes passing through that voxel) for subsequent
% localised utilities. Only voxels saved by the statistical analysis
% which survives the specified thresholds are available for analysis.
%
% A control panel in the interactive window enables interactive
% exploration of the results, divided into volume, cluster, and voxel
% level utilities:
%
% Volume level:
%	  (i) Multiplanar - Launches a MultiPlanar lightbox viewing window
%                                                       - see spm_????.m
%	 (ii)   
%	(iii) Write - write out thresholded statistic image
%	 (iv) 
%	  (v) 
% Cluster level:
%         (i) Maxima - tables all local maxima in the suprathreshold region
%             containing the cursor, going over onto multiple pages if necessary.
%             The table lists all local maxima, their locations, statistic
%             values and p-values. (Note that the cursor is repositioned to the
%             nearest local maxima by this option.)
%                                                       - see spm_maxima.m
%	 (ii)   
%	(iii) 
%	 (iv) 
%	  (v) 
% Voxel level:
%         (i) Plot - Graphs of adjusted and fitted activity against various
%             ordinates. (Note that the cursor is repositioned to the nearest
%             voxel with data by this option.) Additionally, writes out adjusted
%             data to the MatLab command window.
%                                                       - see spm_graph.m
%
%	 (ii) Slices - slices of the thresholded statistic image overlaid
%             on a secondary image chosen by the user. Three transverse
%             slices are shown, being those at the level of the cursor
%             in the z-axis and the two adjacent to it.
%                                                       - see spm_transverse.m
%
%	(iii) Sections - orthogonal sections of the thresholded statistic image
%             overlaid on a secondary image chosen by the user. The sections
%             are through the cursor position.
%                                                       - see spm_sections.m
%	 (iv) 
%	  (v) 
%
% Graphics appear in the bottom half of the graphics window, additional
% controls and questions in the interactive window.
%
%                           ----------------
%
% The MIP uses a template outline in Talairach space. Consequently for
% the results section to display properly the input images to the
% statistics section should either be in Talairach space (with the
% ORIGIN correctly specified), or the ORIGIN header fields should be
% set to the voxel coordinates of the anterior commissure in the input
% images. See spm_format.man ("Data Format" in the help facility) for
% further details of the Analyze image format used by SPM.
%
% Similarly, secondary images should be aligned with the input images
% used for the statistical analysis. In particular the ORIGIN must
% correspond to (0,0,0) in XYZ, the vector of locations.
%
%-----------------------------------------------------------------------
%
% spm_results is an M-file, a script rather than a function, so that
% all the results variables (location, p-values, cluster definitions
% etc.) are available to the (advanced) user in the MatLab workspace.
% Many of the spm_*.m routines invoked are also M-files.
%
%_______________________________________________________________________
% %W% Karl Friston, Andrew Holmes %E%

%-Programers notes
%-----------------------------------------------------------------------
% The GUI is handled by spm_results_ui, which uses an XYZ registry handled
% by spm_XYZreg to synchronise the various GUI location controls.
%
% The interactive MIP is handled by spm_mip_ui.
%_______________________________________________________________________

% Initialise 
%-----------------------------------------------------------------------
clear
WS     = spm('GetWinScale');
Fgraph = spm_figure('FindWin','Graphics');
if isempty(Fgraph), Fgraph=spm_figure('Create','Graphics','Graphics'); end
Finter = spm_figure('FindWin','Interactive');
if isempty(Finter), Finter=spm('CreateIntWin'); end
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

%-Proceed only if there are voxels
%-----------------------------------------------------------------------
if ~length(t); return; end
set(Finter,'Name','SPM results','Pointer','Watch')

%-Load SPM.mat file from appropriate directory
% (CWD, set by spm_projections*_ui)
%-----------------------------------------------------------------------
load([CWD,'/SPM'])
S        = s;
W        = w;
M        = [ [diag(V(4:6)), -(V(7:9).*V(4:6))]; [zeros(1,3) ,1]];
DIM      = V(1:3);

if SPMF, df = Fdf; end

%-Load ER.mat (event-related) file if it exists
%-----------------------------------------------------------------------
str      = [CWD,'/ER.mat'];
if exist(str,'file'); load(str); end



%-Setup Resuts User Interface; Display MIP, design matrix & parameters
%=======================================================================
%-Setup results GUI
hReg = spm_results_ui('SetupGUI',M,DIM);
figure(Finter)

%-Setup MIP & register
%-----------------------------------------------------------------------
figure(Fgraph)
hMIPax = axes('Position',[0.05 0.55 0.55 0.4],'Visible','off');
hMIPax = spm_mip_ui(t,XYZ,V,hMIPax);
spm_XYZreg('XReg',hReg,hMIPax,'spm_mip_ui');

%-Design matrix
%-----------------------------------------------------------------------
hDax = axes('Position',[0.625 0.55 0.35 0.2]);
image((spm_DesMtx('sca',[H C B G],Dnames)+1)*32)
xlabel 'design matrix'
%if SPMZ
%	dy = 0.1/size(CONTRAST,1);
%	for i = 1:size(CONTRAST,1)
%		axes('Position',[0.625 (0.75 +dy*(i-1)) 0.35 dy])
%		h = bar(CONTRAST(i,:),1);
%		set(h,'FaceColor',[1 1 1]*.8)
%		set(gca,'XLim',[0.5,size(CONTRAST,2)+0.5],'Visible','off')
%		text(size(CONTRAST,2)+.55,0,sprintf('%d',i),...
%			'FontSize',10,'Color','r')
%	end
%	axes('Position',[0.6 0.75 0.025 0.1],'Visible','off')
%	text(0.5,0.5,'contrast',...
%		'FontSize',10,'Rotation',90,...
%		'HorizontalAlignment','Center','VerticalAlignment','middle')
%end

%-Text
%-----------------------------------------------------------------------
hTax = axes('Position',[0.6 0.9 0.4 0.05],'Visible','off');
if SPMZ
	SPMdist = 'Z';
	tmp = 1 - spm_Ncdf(U);
else
	SPMdist = 'F';
	tmp = 1 - spm_Fcdf(U,df);
end
text(0.5,1,['SPM{',SPMdist,'}'],...
	'FontSize',16,'Fontweight','Bold','HorizontalAlignment','Center')
text(0,0.6,spm_str_manip(CWD,'a30'),...
	'FontSize',10,'FontWeight','Bold')
text(0,0.3,sprintf('Height threshold {u} = %0.2f, p = %0.3f',U,tmp),...
	'FontSize',10)
text(0,0,sprintf('Extent threshold {k} = %0.0f voxels',k),...
	'FontSize',10)


% Finished
%=======================================================================
set(Finter,'Pointer','Arrow')

















