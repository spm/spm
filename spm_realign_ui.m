function spm_realign_ui(arg1)
% Within Mode Image Realignment
%___________________________________________________________________________
%
% This routine realigns a time-series of images acquired from the same subject 
% using a least squares approach and a 6 parameter (rigid body) spatial 
% transformation.  The first image in the list specified by the user is used
% as a reference to which all subsequent scans are realigned.  The reference
% scan does not have to the the first chronologically and it may be wise to
% chose a 'representative scan' in this role.
%
% For fMRI data an additional adjustment is made to the data that removes
% a tiny amount of the movement-related confounds of these effects.
% However, it may be preferable to include the functions of the estimated
% movement parameters as confounds in the statistics part.
%
%
% Uses
% Primarily to remove movement artefact in fMRI and PET time-series (or more 
% generally longitudinal studies)
%
%
% Inputs
% A series of *.img conforming to SPM data format (see 'Data Format').  The 
% relative displacement of the images should be small with respect to their 
% resolution.  This is usually easy to ensure for functional images (e.g. 
% fMRI, PET SPECT).
%
%
% Outputs
% The parameter estimation part writes out ".mat" files for each of the
% input images.  The part of the routine that writes the resliced images
% uses information in these ".mat" files and writes the realigned *.img
% files to the same subdirectory prefixed with an 'r' (i.e. r*.img).  The
% details of the transformation are displayed in the results window as
% plots of translation and rotation.
% A set of realignment parameters are saved for each session, named:
% realignment_params_*.txt.
%
%__________________________________________________________________________
% Refs:
%
% Friston KJ, Ashburner J, Frith CD, Poline J-B, Heather JD & Frackowiak
% RSJ (1995) Spatial registration and normalization of images Hum. Brain
% Map. 2:165-189
%
% Friston KJ, Williams SR, Howard R Frackowiak RSJ and Turner R (1995)
% Movement-related effect in fMRI time-series.  Mag. Res. Med. 35:346-355
%
% W. F. Eddy, M. Fitzgerald and D. C. Noll (1996) Improved Image
% Registration by Using Fourier Interpolation. Mag. Res. Med. 36(6):923-931
%
% R. W. Cox and A. Jesmanowicz (1999)  Real-Time 3D Image Registration
% for Functional MRI.  Submitted to MRM (April 1999) and avaliable from:
% http://varda.biophysics.mcw.edu/~cox/index.html.
%
%__________________________________________________________________________
%
% --- The Prompts Explained ---
%
% 'number of subjects'
% Enter the number of subjects you wish to realign.
%
% For fMRI, it will ask you the number of sessions for each subject.
% In the coregistration step, the sessions are first realigned to
% each other, by aligning the first scan from each session to the
% first scan of the first session.  Then the images within each session
% are aligned to the first image of the session.
% The parameter estimation is performed this way because it is assumed
% (rightly or not) that there may be systematic differences
% in the images between sessions.
% The adjustment step (correcting for resampling artifacts) is also
% performed completely independantly between each of the fMRI sessions.
%
% 'select scans for subject ..'
% Select the scans you wish to realign. All operations are relative
% to the first image selected.
%
% ......... Note that not all of the following prompts may be used: .........
%
% 'Which option?'
% 	'Coregister only'
% 	Only determine the parameters required to transform each of the
% 	images 2..n to the same space as image 1.
% 	The determined parameters for image XXXX.img are saved in the
%	file XXXX.mat. This is a Matlab file, containing the matrix 'M'.
% 	The location of an image voxel (in mm) can be determined by
% 	computing M(1:3,:)*[xcoord ycoord zcoord 1].
%	Note that if the coregistration is performed more than once on
%	the unresliced data, the starting estimates are obtained from
% 	parameters stored in the '.mat' files.
%	Note that for PET, the coregistration is a two step process.
%	First of all, the images are all realigned to the first in
%	the series.  A mean of these realigned images is created, and
%	a second pass realignment is performed to realign all the
%	images to the mean. Finally, the parameters are corrected
%	for any differences estimated by registering the first image in
%	the series to the mean image.
%
% 	'Reslice Only'
% 	Reslice the specified images according to the contents of the
% 	previously determined parameters. The images are resliced to be
% 	in the same space as the first one selected.  For fMRI, this is
%	the first image of the first session.
%
% 	'Coregister & Reslice'
% 	Combine the above two steps together.
%
%
% Options for reslicing:
%
% 'Create what?'
% 	'All Images (1..n)'
% 	This reslices all the images - including the first image selected
% 	- which will remain in it's original position.
%
%	'Images 2..n'
% 	Reslices images 2..n only. Useful for if you wish to reslice
% 	(for example) a PET image to fit a structural MRI, without
% 	creating a second identical MRI volume.
%
%	'All Images + Mean Image'
% 	In addition to reslicing the images, it also creates a mean of the
% 	resliced image.
%
%	'Mean Image Only'
% 	Creates the mean image only.
%
% 'Mask the images?'
% To avoid artifactual movement-related variance.
% Because of subject motion, different images are likely to have different
% patterns of zeros from where it was not possible to sample data.
% With masking enabled, the program searches through the whole time series
% looking for voxels which need to be sampled from outside the original
% images. Where this occurs, that voxel is set to zero for the whole set
% of images (unless the image format can represent NaN, in which case
% NaNs are used where possible).
%
% 'Adjust sampling errors?' (fMRI only)
% Adjust the data (fMRI) to remove interpolation errors arising from the
% reslicing of the data.  The adjustment for each fMRI session is performed
% independantly of any other session.  Bayesian statistics are used to
% attempt to regularize the adjustment in order to prevent an excessive
% amount of signal from being removed.  A priori variances for coefficients
% are assumed to be stationary and are estimated by translating the first
% image by a number of different distances using both Fourier and sinc
% interpolation.  This gives a ball park figure on how much error is
% likely to arise because of the approximations in the sinc interpolation.
% The certainty of the solution is obtained from the residuals after
% fitting the optimum linear combination of the basis functions through
% the data.  Estimates of certainty based on the residuals are
% unfortunately just an approximation.   
% We still don't fully understand the nature of the movement artifacts
% that arise using fMRI.  The current model is simply attempting to remove
% interpolation errors.  There are many other sources of error that the
% model does not attempt to remove.
% It is possible that adjusting the data without taking into account
% the design matrix for the statistics may be problematic when there are
% stimulous correlated movements, since adjusting seperately requires the
% assumption that the movements are independant from the paradigm.  It
% MAY BE BE BETTER TO INCLUDE THE ESTIMATED MOTION PARAMETERS AS CONFOUNDS
% WHEN THE STATISTICS ARE RUN.  The motion parameters are saved for each
% session, so this should be easily possible.
%
% 'Reslice Interpolation Method?'
% 	'Trilinear Interpolation'
% 	Use trilinear interpolation (first-order hold) to sample the images
%       during the writing of realigned images.
%
% 	'Sinc Interpolation'
% 	Use a sinc interpolation to sample the images during the writing
%	of realigned images.
% 	This is slower than bilinear interpolation, but produces better
% 	results. It is especially recommended for fMRI time series.
%	An 9x9x9 kernel is used to resample the images. 
%
%	'Fourier space Interpolation' (fMRI only)
%	Rigid body rotations are executed as a series of shears, which
%	are performed in Fourier space (Eddy et. al. 1996).  This routine
%	only supports cubic voxels (since zooms can not be done by
%	convolution in Fourier space).
%	No adjustment is available for this.
%
%__________________________________________________________________________
%
% The `.mat' files.
%
% This simply contains a 4x4 affine transformation matrix in a variable `M'.
% These files are normally generated by the `realignment' and
% `coregistration' modules.  What these matrixes contain is a mapping from
% the voxel coordinates (x0,y0,z0) (where the first voxel is at coordinate
% (1,1,1)), to coordinates in millimeters (x1,y1,z1).  By default, the
% the new coordinate system is derived from the `origin' and `vox' fields
% of the image header.
%  
% x1 = M(1,1)*x0 + M(1,2)*y0 + M(1,3)*z0 + M(1,4)
% y1 = M(2,1)*x0 + M(2,2)*y0 + M(2,3)*z0 + M(2,4)
% z1 = M(3,1)*x0 + M(3,2)*y0 + M(3,3)*z0 + M(3,4)
%
% Assuming that image1 has a transformation matrix M1, and image2 has a
% transformation matrix M2, the mapping from image1 to image2 is: M2\M1
% (ie. from the coordinate system of image1 into millimeters, followed
% by a mapping from millimeters into the space of image2).
%
% These `.mat' files allow several realignment or coregistration steps to be
% combined into a single operation (without the necessity of resampling the
% images several times).  The `.mat' files are also used by the spatial
% normalisation module.
%__________________________________________________________________________
% %W% John Ashburner - with input from Oliver Josephs %E%

% Programmers Guide
% Batch system implemented on this routine. See spm_bch.man
% If inputs are modified in this routine, try to modify spm_bch.man
% and spm_bch_bchmat (if necessary) accordingly. 
% Calls to spm_input in this routine use the BCH gobal variable.  
%    BCH.bch_mat 
%    BCH.index0  = {'realign',index_of_Analysis};
% or 
%    BCH.index0  = {'RealignCoreg',index_of_Analysis}; (when
%                   spm_realign is launched for edit_defaults 
%
%_______________________________________________________________________


global MODALITY sptl_WhchPtn sptl_DjstFMRI sptl_CrtWht sptl_MskOptn SWD
global BCH sptl_RlgnQlty sptl_WghtRg;

if (nargin == 0)
	% User interface.
	%_______________________________________________________________________
	SPMid = spm('FnBanner',mfilename,'2.21');
	[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Realign');
	spm_help('!ContextHelp','spm_realign_ui.m');

	pos = 1;

	n     = spm_input('number of subjects', pos, 'e', 1,...
                          'batch',{},'subject_nb');
	if (n < 1)
		spm_figure('Clear','Interactive');
		return;
	end

	P = cell(n,1);
	pos = pos + 1;
	for i = 1:n,
		if strcmp(MODALITY,'FMRI'),
			ns = spm_input(['num sessions for subject ' num2str(i)], pos,...
                                       'e', 1,'batch',{},'num_sessions');
			pp = cell(1,ns);
			for s=1:ns,
				p = '';
				while size(p,1)<1,
					if isempty(BCH),
						p = spm_get(Inf,'.img',...
						['scans for subj ' num2str(i) ', sess' num2str(s)]);
					else,
						p = spm_input('batch',{'sessions',i},'images',s);
					end;
				end;
				pp{s} = p;
			end;
			P{i} = pp;
		else, %- no batch mode for 'PET'
			p  = cell(1,1);
			p{1} = '';
			while size(p{1},1)<1,
			      p{1} = spm_get(Inf,'.img',...
				  ['select scans for subject ' num2str(i)]);
			end;
			P{i} = p;
		end;
	end;


	if strcmp(MODALITY,'PET'),
		FlagsC = struct('quality',sptl_RlgnQlty,'fwhm',8,'rtm',[]);
	else,
		FlagsC = struct('quality',sptl_RlgnQlty,'fwhm',6);
	end;

	if sptl_WhchPtn == 1,
		WhchPtn = 3;
	else,
		WhchPtn = spm_input('Which option?', pos, 'm',...
			'Coregister only|Reslice Only|Coregister & Reslice',...
			[1 2 3],3,'batch',{},'option');
		pos = pos + 1;
	end;

	PW = '';
	if (WhchPtn == 1 | WhchPtn == 3) & sptl_WghtRg,
		if spm_input(...
			['Weight the reference image(s)?'],...
			2, 'm',...
			['Dont weight registration|'...
			 'Weight registration'], [0 1], 1,...
			 'batch',{},'weight_reg'),

			if isempty(BCH),
				PW = spm_get(n,'.img',...
					'Weight images for each subj');
			else,
				PW = spm_input('batch',{'sessions',i},'weights',s);
			end;
		end;
	end;

	% Reslicing options
	%-----------------------------------------------------------------------
	if WhchPtn == 2 | WhchPtn == 3,
		FlagsR = struct('hold',1,'mask',0,'which',2,'mean',1);
		FlagsR.hold = spm_input('Reslice interpolation method?',pos,'m',...
			     'Trilinear Interpolation|Sinc Interpolation|Fourier space Interpolation',...
			     [1 -9 Inf],2,'batch',{},'reslice_method');
		pos = pos + 1;

		if sptl_CrtWht == 1,
			p = 3;
		else
			p = spm_input('Create what?',pos,'m',...
				[' All Images (1..n)| Images 2..n|'...
				 ' All Images + Mean Image| Mean Image Only'],...
				[1 2 3 4],3,'batch',{},'create');
			pos = pos + 1;
		end
		if (p == 1) FlagsR.which = 2; FlagsR.mean = 0; end
		if (p == 2) FlagsR.which = 1; FlagsR.mean = 0; end
		if (p == 3) FlagsR.which = 2; FlagsR.mean = 1; end
		if (p == 4) FlagsR.which = 0; FlagsR.mean = 1; end
		if FlagsR.which > 0,
			if sptl_MskOptn == 1,
				FlagsR.mask = 1;
			else,
				if spm_input('Mask the resliced images?',pos,'y/n',...
                                              'batch',{},'mask') == 'y',
					FlagsR.mask = 1;
				end;
				pos = pos + 1;
			end;
			if strcmp(MODALITY, 'FMRI'),
				if finite(FlagsR.hold),
					if sptl_DjstFMRI == 1,
						FlagsR.fudge = 1;
					elseif sptl_DjstFMRI ~= 0,
						if spm_input(...
							'Adjust sampling errors?',pos,'y/n','batch',...
                                                        {},'adjust_sampling_errors') == 'y',
						FlagsR.fudge = 1;
						end;
						pos = pos + 1;
					end;
				end
			end
		end
	end

	spm('Pointer','Watch');
	for i = 1:n
		spm('FigName',['Realign: working on subject ' num2str(i)],Finter,CmdLine);
		fprintf('\rRealigning Subject %d: ', i);
		if WhchPtn==1 | WhchPtn==3,
			flagsC = FlagsC;
			if ~isempty(PW), flagsC.PW = deblank(PW(i,:)); end;
			spm_realign(P{i},flagsC);
		end
		if WhchPtn==2 | WhchPtn==3,
			spm_reslice(P{i},FlagsR)
		end;
	end
	fprintf('\r%60s%s', ' ',sprintf('\b')*ones(1,60));
	spm('FigName','Realign: done',Finter,CmdLine);
	spm('Pointer');
	return;

elseif nargin == 1 & strcmp(arg1,'Defaults'),
	edit_defaults;
	return;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function edit_defaults
global MODALITY sptl_WhchPtn sptl_DjstFMRI sptl_CrtWht sptl_MskOptn SWD
global sptl_RlgnQlty sptl_WghtRg

%- in batch mode, the top level variable here is 'RealignCoreg',iA

tmp = 1;
if sptl_WhchPtn == 1, tmp = 2; end;
sptl_WhchPtn = spm_input(...
	['Coregistration and reslicing?'],...
	2, 'm',...
	['Allow separate coregistration and reslicing|'...
	 'Combine coregistration and reslicing'], [-1 1], tmp,...
         'batch',{},'separate_combine');

tmp = 2;
if sptl_CrtWht == 1,
	tmp = 1;
end
sptl_CrtWht   = spm_input(['Images to create?'], 3, 'm',...
	       'All Images + Mean Image|Full options', [1 -1], tmp,...
               'batch',{},'create');

tmp = 3;
if sptl_DjstFMRI == 1,
	tmp = 1;
elseif sptl_DjstFMRI == 0
	tmp = 2;
end
sptl_DjstFMRI = spm_input(['fMRI adjustment for interpolation?'],4,'m',...
	          '   Always adjust |    Never adjust|Optional adjust',...
	          [1 0 -1], tmp,'batch',{},'adjust');

tmp = 2;
if sptl_MskOptn == 1,
	tmp = 1;
end
sptl_MskOptn  = spm_input(['Option to mask images?'], 5, 'm',...
		'  Always mask|Optional mask', [1 -1], tmp,...
	        'batch',{},'mask');

tmp2 = [1.00 0.90 0.75 0.50 0.25 0.10 0.05 0.01 0.005 0.001];
tmp = find(sptl_RlgnQlty == tmp2);
if isempty(tmp) tmp = length(0.5); end
sptl_RlgnQlty = spm_input('Registration Quality?','+1','m',...
	['Quality 1.00  (slowest/most accurate) |Quality 0.90|' ...
	 'Quality 0.75|Quality 0.50|Quality 0.25|Quality 0.10|' ...
	 'Quality 0.05|Quality 0.01|' ...
	 'Quality 0.005|Quality 0.001 (fastest/poorest)'],tmp2, tmp,...
                'batch',{},'reg_quality');


tmp = 0; if sptl_WghtRg == 1, tmp = 1; end;
sptl_WghtRg = spm_input(...
	['Allow weighting of reference image?'],...
	'+1', 'm',...
	['Allow weighting|'...
	 'Dont allow weighting'], [1 0], tmp,...
         'batch',{},'weight_reg');

return;
%_______________________________________________________________________
