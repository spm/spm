
% Display and analysis of regional effects
% FORMAT spm_results
%_______________________________________________________________________
%
% The SPM results section is for the interactive exploration and
% characterisation of the results of a statistical analysis.
% 
% The user is prompted to select a SPM{Z}, SPM{T} or SPM{F}, that is
% thresholded at user specified levels. The specification of the
% contrasts to use and the height and size thresholds are the same as
% that described in spm_getSPM.m. The resulting SPM is then
% displayed in the graphics window as a maximum intensity projection,
% alongside the design matrix and contrasts employed.
%
% The cursors in the MIP can be moved (dragged) to select a particular
% voxel. The three mouse buttons give different drag and drop behaviour:
% Button 1 - point & drop; Button 2 - "dynamic" drag & drop with
% co-ordinate & SPM value updating; Button 3 - "magnetic" drag & drop,
% where the cursor jumps to the nearest suprathreshold voxel in the
% MIP, and shows the value there. See spm_mip_ui.m, the MIP GUI handling
% function.
%
% The current voxel specifies the voxel, suprathreshold cluster, or
% orthogonal planes (planes passing through that voxel) for subsequent
% localised utilities. Only voxels saved by the statistical analysis
% that survive the specified thresholds are available for analysis.
%
% A control panel in the interactive window enables interactive
% exploration of the results and is divided into volume, cluster, and 
% voxel level utilities:
%
% Volume level:
%	  (i) List - tabulates p values and statistics
%						- see spm_list.m
%
%	 (ii) Multiplanar - Launches a MultiPlanar lightbox viewing window
%                                               - see spm_????.m
%
%	(iii) Write - write out thresholded SPM
%                                    		- see spm_write_filtered.m
%	 (iv) 
%	  (v)
%
% Cluster level:
%         (i) Maxima - tables all local maxima in the suprathreshold region
%             containing the cursor, going onto multiple pages if necessary.
%             The table lists all local maxima, their locations, statistic
%             values and p-values. (Note that the cursor is repositioned to the
%             nearest local maxima by this option.)
%                                  		- see spm_maxima.m
%
%	 (ii) Small Volume Correction (S.V.C) computes a p value corrected
%             for a small specified volume of interest for, and centred on,
%             the current voxel.
%                                 		- see spm_VOI.m
%	(iii) 
%	 (iv) 
%	  (v)
%
% Voxel level:
%         (i) Plot - Graphs of adjusted and fitted activity against various
%             ordinates. (Note that the cursor is repositioned to the nearest
%             voxel with data by this option.) Additionally, writes out 
%             adjusted data to the MatLab command window.
%                                 		- see spm_graph.m
%
%	 (ii) Slices - slices of the thresholded statistic image overlaid
%             on a secondary image chosen by the user. Three transverse
%             slices are shown, being those at the level of the cursor
%             in the z-axis and the two adjacent to it.
%                                 		- see spm_transverse.m
%
%	(iii) Sections - orthogonal sections of the thresholded statistic image
%             overlaid on a secondary image chosen by the user. The sections
%             are through the cursor position.
%                                    		- see spm_sections.m
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
%---------------------------------------------------------------------------
%
% spm_results is an M-file, a script rather than a function, so that
% the results structures (SPM, VOL and DES) are available to the (advanced)
% user in the MatLab workspace.
%
%___________________________________________________________________________
% %W% Karl Friston, Andrew Holmes %E%

%-Programers notes
%---------------------------------------------------------------------------
% The GUI is handled by spm_results_ui, which uses an XYZ registry handled
% by spm_XYZreg to synchronise the various GUI location controls.
%
% The interactive MIP is handled by spm_mip_ui.
%___________________________________________________________________________

% Initialise 
%---------------------------------------------------------------------------
Fgraph = spm_figure('FindWin','Graphics');
Finter = spm_figure('FindWin','Interactive');
spm_clf(Fgraph)
spm_clf(Finter)
set(Finter,'Name','SPM results')
global CWD

%-Get thresholded data and parameters of design
%===========================================================================

%-Get SPM
%---------------------------------------------------------------------------
[SPM,VOL,DES] = spm_getSPM;


%-Setup Results User Interface; Display MIP, design matrix & parameters
%===========================================================================

%-Setup results GUI
%---------------------------------------------------------------------------
hReg   = spm_results_ui('SetupGUI',VOL.M,VOL.DIM);


%-Setup Maximium intensity projection (MIP) & register
%---------------------------------------------------------------------------
figure(Fgraph)
hMIPax = axes('Position',[0.05 0.55 0.55 0.4],'Visible','off');
hMIPax = spm_mip_ui(SPM.Z,VOL.XYZ,[VOL.DIM;VOL.VOX;VOL.ORG],hMIPax);
spm_XYZreg('XReg',hReg,hMIPax,'spm_mip_ui');
text(260,260,['SPM{' SPM.STAT '}'],'FontSize',16,'Fontweight','Bold')


%-Design matrix & contrast
%---------------------------------------------------------------------------
axes('Position',[0.65 0.55 0.25 0.25],'Tag','Empty')
image(spm_DesMtx('sca',DES.X)*64)
xlabel 'Design matrix'
dy     = 0.1/size(DES.C,1);
for  i = 1:size(DES.C,1)
	axes('Position',[0.65 (0.8 + dy*(i - 1)) 0.25 dy])
	bar(DES.C(i,:));
	axis tight, axis off
end
str    = str2mat(spm_str_manip(CWD,'a30'),...
	sprintf('Height threshold {u} = %0.2f',SPM.u),...
	sprintf('Extent threshold {k} = %0.0f resels',SPM.k));
title(str,'FontSize',9)


% Finished
%===========================================================================
set(Finter,'Pointer','Arrow')

















