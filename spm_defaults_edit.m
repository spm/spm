function spm_defaults_edit(arg1, arg2)
% Modify defaults
% FORMAT spm_defaults_edit
%_______________________________________________________________________
%
% spm_defaults_edit allows the current defaults to be edited.
%
% These changes do not persist across sessions. SPMs startup defaults
% are specified in the first spm_defaults on the MATLABPATH.
%
% The defaults which can be modified are:
% 
% Printing Options
%     Allows a number of different printing defaults to be specified.
% 
% Miscellaneous Defaults
%     This includes:
%     * Specification of a file for logging dialogue between
%       the user and SPM.
%     * Command line input option. Rather than clicking
%       buttons on the interface, input can be typed to
%       the Matlab window.
%     * The intensity of any grid which superimposed on any
%       displayed images.
%     * Specification of paging option for tabular output of cluster
%       statistics, enabling more exhaustive classifications
% 
% Header Defaults (for the currnet Modality - PET or fMRI)
%     The values to be taken as default when there are no Analyze
%     image headers. There are two different sets which depend on
%     the modality in which SPM is running.
%     * image size in x,y and z {voxels}
%     * voxel size in x,y and z {mm}
%     * scaling co-efficient applied to *.img data on entry
%       into SPM. 
%     * data type.  (see spm_type.m for supported types
%       and specifiers)
%     * offest of the image data in file {bytes}
%     * the voxel corresponding the [0 0 0] in the location
%       vector XYZ
%     * a string describing the nature of the image data.
% 
% Realignment & Coregistration Defaults
%     An assortment of defaults.
%
% Spatial Normalisation Defaults
%     An assortment of defaults.
%
% The 'reset' option re-loads the startup defaults from spm_defaults.m
%
%_______________________________________________________________________
% @(#)spm_defaults_edit.m	1.9 John Ashburner 96/09/27

global MODALITY
global CWD PRINTSTR LOGFILE CMDLINE GRID proj_MultiPage
global UFp DIM VOX TYPE SCALE OFFSET ORIGIN DESCRIP
global PET_UFp PET_DIM PET_VOX PET_TYPE PET_SCALE PET_OFFSET PET_ORIGIN PET_DESCRIP
global fMRI_UFp fMRI_DIM fMRI_VOX fMRI_TYPE fMRI_SCALE fMRI_OFFSET fMRI_ORIGIN fMRI_DESCRIP

if nargin == 0
	spm_figure('Clear','Interactive');
	set(spm_figure('FindWin','Interactive'),...
		'Name','Defaults Edit');
	spm_help('!ContextHelp','spm_defaults_edit.m')
	pos = 1;

	callbacks = str2mat(...
		'spm_defaults_edit(''Printing'');',...
		'spm_defaults_edit(''Misc'');',...
		'spm_defaults_edit(''Hdr'');',...
		'spm_realign(''Defaults'');',...
		'spm_sn3d(''Defaults'');',...
		'spm_defaults_edit(''Statistics'');',...
		'spm_defaults_edit(''Reset'');'...
		);

	a1 = spm_input('Defaults Area?',pos,'m',...
		['Printing Options|'...
		 'Miscellaneous Defaults|'...
		 'Header Defaults - ',MODALITY,'|'...
		 'Realignment & Coregistration|'...
		 'Spatial Normalisation|'...
		 'Statistics - ',MODALITY,'|'...
		 'Reset All']...
		);

	eval(deblank(callbacks(a1,:)));
	spm_figure('Clear','Interactive');

elseif strcmp(arg1, 'Directory')

	% Default directory for results files etc..
	%---------------------------------------------------------------
	CWD = deblank(spm_get(-1,'*','Directory'));
	chdir(CWD);


elseif strcmp(arg1, 'Misc')

	% Miscellaneous
	%---------------------------------------------------------------
	if ~isempty(LOGFILE), tmp='yes'; def=1; else, tmp='no'; def=2; end
	if spm_input(['Log to file? (' tmp ')'],2,'y/n',[1,0],def)
		LOGFILE = ...
			deblank(spm_input('Logfile Name:',2,'s', LOGFILE));
	else
		LOGFILE = '';
	end

	if CMDLINE ~= 0, tmp='yes'; def=1; else, tmp='no'; def=2; end
	CMDLINE = ...
	    spm_input(['Command Line Input (' tmp ')?'],3,'y/n',[1,0],def);
	GRID = spm_input('Grid value (0-1):', 4, 'e', GRID);

	if proj_MultiPage, tmp='yes'; def=1; else, tmp='no'; def=2; end
	proj_MultiPage = ...
		spm_input(['Paging of stats (' tmp ')?'],5,'y/n',[1,0],def);

elseif strcmp(arg1, 'Printing')

	% Printing Defaults
	%---------------------------------------------------------------
	a0 = spm_input('Printing Mode?', 2, 'm', [...
			'Postscript to File|'...
			'Postscript to Printer|'...
			'Other Format to File|'...
			'Custom'...
			]);
	if (a0 == 1)
		fname = date; fname(find(fname=='-')) = []; fname = ['spmfig_' fname];
		fname = spm_str_manip(spm_input('Postscript filename:',3,'s', fname),'rtd');
		a1    = spm_input('Postscript Type?', 4, 'm', [...
			'PostScript for black and white printers|'...
			'PostScript for colour printers|'...
			'Level 2 PostScript for black and white printers|'...
			'Level 2 PostScript for colour printers|'...
			'Encapsulated PostScript (EPSF)|'...
			'Encapsulated Colour PostScript (EPSF)|'...
			'Encapsulated Level 2 PostScript (EPSF)|'...
			'Encapsulated Level 2 Color PostScript (EPSF)|'...
			'Encapsulated        with 1-bit preview (EPSI)|'...
			'Encapsulated Colour with 1-bit preview (EPSI)|'...
			'Encapsulated Level 2   w 1-bit preview (EPSI)|'...
			'Encapsulated Level 2 Colour  w preview (EPSI)|'...
			]);
		prstr1 = str2mat(...
			['print(''-dps'' ,''-append'',''' fname '.ps'');'],...
			['print(''-dpsc'',''-append'',''' fname '.ps'');'],...
			['print(''-dps2'',''-append'',''' fname '.ps'');'],...
			['print(''-dpsc2'',''-append'',''' fname '.ps'');']);
		prstr1 = str2mat(prstr1,...
			['global PAGENUM;if isempty(PAGENUM),PAGENUM = 1;end;'...
			 'print(''-deps'',[''' fname '_'' num2str(PAGENUM) ''.ps'']); PAGENUM = PAGENUM + 1;'],...
			['global PAGENUM;if isempty(PAGENUM),PAGENUM = 1;end;'...
			 'print(''-depsc'',[''' fname '_'' num2str(PAGENUM) ''.ps'']); PAGENUM = PAGENUM + 1;'],...
			['global PAGENUM;if isempty(PAGENUM),PAGENUM = 1;end;'...
			 'print(''-deps2'',[''' fname '_'' num2str(PAGENUM) ''.ps'']); PAGENUM = PAGENUM + 1;'],...
			['global PAGENUM;if isempty(PAGENUM),PAGENUM = 1;end;'...
			 'print(''-depsc2'',[''' fname '_'' num2str(PAGENUM) ''.ps'']); PAGENUM = PAGENUM + 1;'],...
			['global PAGENUM;if isempty(PAGENUM),PAGENUM = 1;end;'...
			 'print(''-deps'',''-epsi'',[''' fname '_'' num2str(PAGENUM) ''.ps'']); PAGENUM = PAGENUM + 1;'],...
			['global PAGENUM;if isempty(PAGENUM),PAGENUM = 1;end;'...
			 'print(''-depsc'',''-epsi'',[''' fname '_'' num2str(PAGENUM) ''.ps'']); PAGENUM = PAGENUM + 1;'],...
			['global PAGENUM;if isempty(PAGENUM),PAGENUM = 1;end;'...
			 'print(''-deps2'',''-epsi'',[''' fname '_'' num2str(PAGENUM) ''.ps'']); PAGENUM = PAGENUM + 1;'],...
			['global PAGENUM;if isempty(PAGENUM),PAGENUM = 1;end;'...
			 'print(''-depsc2'',''-epsi'',[''' fname '_'' num2str(PAGENUM) ''.ps'']); PAGENUM = PAGENUM + 1;']);
		PRINTSTR = deblank(prstr1(a1,:));
	elseif (a0 == 2)
		printer = '';
		if (spm_input('Default Printer?', 3, 'y/n') == 'n')
			printer = spm_input('Printer Name:',3,'s');
			printer = [' -P' printer];
		end
		a1 = spm_input('Postscript Type:',4,'b','B & W|Colour', str2mat('-dps', '-dpsc'));
		PRINTSTR = ['print ' a1 printer];
	elseif (a0 == 3)
		fname = date; fname(find(fname=='-')) = []; fname = ['spmfig_' fname];
		fname = spm_str_manip(spm_input('Graphics filename:',3,'s', fname),'rtd');
		a1    = spm_input('Graphics Type?', 4, 'm', [...
			'HPGL compatible with Hewlett-Packard 7475A plotter|'...
			'Adobe Illustrator 88 compatible illustration file|'...
			'M-file (and Mat-file, if necessary)|'...
			'Baseline JPEG image|'...
			'TIFF with packbits compression|'...
			'Color image format|'...
			]);
		prstr1 = str2mat(...
			['global PAGENUM;if isempty(PAGENUM),PAGENUM = 1;end;'...
			 'print(''-dhpgl'',[''' fname '_'' num2str(PAGENUM) ''.hpgl'']); PAGENUM = PAGENUM + 1;'],...
			['global PAGENUM;if isempty(PAGENUM),PAGENUM = 1;end;'...
			 'print(''-dill'',[''' fname '_'' num2str(PAGENUM) ''.ill'']); PAGENUM = PAGENUM + 1;'],...
			['global PAGENUM;if isempty(PAGENUM),PAGENUM = 1;end;'...
			 'print(''-dmfile'',[''' fname '_'' num2str(PAGENUM) ''.m'']); PAGENUM = PAGENUM + 1;'],...
			['global PAGENUM;if isempty(PAGENUM),PAGENUM = 1;end;'...
			 'print(''-djpeg'',[''' fname '_'' num2str(PAGENUM) ''.jpeg'']); PAGENUM = PAGENUM + 1;'],...
			['global PAGENUM;if isempty(PAGENUM),PAGENUM = 1;end;'...
			 'print(''-dtiff'',[''' fname '_'' num2str(PAGENUM) ''.tiff'']); PAGENUM = PAGENUM + 1;'],...
			['global PAGENUM;if isempty(PAGENUM),PAGENUM = 1;end;'...
			 'print(''-dtiffnocompression'',[''' fname '_'' num2str(PAGENUM) ''.tiff'']); PAGENUM = PAGENUM + 1;']);
		PRINTSTR = deblank(prstr1(a1,:));
	else
		PRINTSTR = spm_input('Print String',3,'s');
	end

elseif strcmp(arg1, 'Hdr')

	% Header Defaults
	%---------------------------------------------------------------

	n = 0;
	while n ~= 3
		tmp      = spm_input('Image size {voxels}',2,'s',...
			[num2str(DIM(1)) ' ' num2str(DIM(2)) ' ' num2str(DIM(3))]);
		[dim, n] = sscanf(tmp,'%d %d %d');
	end
	DIM = dim;

	n = 0;
	while n ~= 3
		tmp      = spm_input('Voxel Size {mm}',3,'s',...
			[num2str(VOX(1)) ' ' num2str(VOX(2)) ' ' num2str(VOX(3))]);
		[vox, n] = sscanf(tmp,'%g %g %g');
	end
	VOX = vox;

	SCALE = spm_input('Scaling Coefficient',4,'e',[SCALE]);

	type_val = [2 4 8 16 64];
	type_str = str2mat('Unsigned Char','Signed Short','Signed Integer','Floating Point','Double Precision');
	TYPE = spm_input(['Data Type (' deblank(type_str(find(type_val==TYPE),:)) ')'],5,'m',...
		'Unsigned Char	(8  bit)|Signed Short	(16 bit)|Signed Integer	(32 bit)|Floating Point|Double Precision',...
		[2 4 8 16 64]);
	OFFSET = spm_input('Offset  {bytes}',6,'e',[OFFSET]);
	n = 0;
	while n ~= 3
		tmp      = spm_input('Origin {voxels}',7,'s',...
			[num2str(ORIGIN(1)) ' ' num2str(ORIGIN(2)) ' ' num2str(ORIGIN(3))]);
		[origin, n] = sscanf(tmp,'%d %d %d');
	end
	ORIGIN = origin;
	DESCRIP = spm_input('Description',8,'s', DESCRIP);

	if strcmp(MODALITY,'PET')
		PET_DIM       = DIM;
		PET_VOX       = VOX;
		PET_TYPE      = TYPE;
		PET_SCALE     = SCALE;
		PET_OFFSET    = OFFSET;
		PET_ORIGIN    = ORIGIN;
		PET_DESCRIP   = DESCRIP;
	elseif strcmp(MODALITY,'FMRI')
		fMRI_DIM      = DIM;
		fMRI_VOX      = VOX;
		fMRI_TYPE     = TYPE;
		fMRI_SCALE    = SCALE;
		fMRI_OFFSET   = OFFSET;
		fMRI_ORIGIN   = ORIGIN;
		fMRI_DESCRIP  = DESCRIP;
	end

elseif strcmp(arg1, 'Statistics')
	UFp = spm_input('Upper tail F prob. threshold',2,'e',UFp);
	if strcmp(MODALITY,'PET')
		PET_UFp       = UFp;
	elseif strcmp(MODALITY,'FMRI')
		fMRI_UFp      = UFp;
	end


elseif strcmp(arg1, 'Reset')
	if exist('spm_defaults')==2
		spm_defaults;
	end
	spm('chmod',MODALITY);
end
