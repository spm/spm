function spm_sn3d(P,matname,bb,Vox,params,spms,brainmask,objmask)
% Spatial (stereotactic) normalization
% ___________________________________________________________________________
% 
% This module spatially (stereotactically) normalizes MRI, PET or SPECT
% images into a standard space defined by some ideal model or template
% image[s].  The template images supplied with SPM conform to the space
% defined by the ICBM, NIH P-20 project, and approximate that of the
% the space described in the atlas of Talairach and Tournoux (1988).
% The transformation can also be applied to any other image that has
% been coregistered with these scans.
% 
% 
% Mechanism
% Generally, the algorithms work by minimising the sum of squares
% difference between the image which is to be normalised, and a linear
% combination of one or more template images.  For the least squares
% registration to produce an unbiased estimate of the spatial
% transformation, the image contrast in the templates (or linear
% combination of templates) should be similar to that of the image from
% which the spatial normalization is derived.  The registration simply
% searches for an optimum solution.  If the starting estimates are not
% good, then the optimum it finds may not find the global optimum.
% 
% The first step of the normalization is to determine the optimum
% 12-parameter affine transformation.  Initially, the registration is
% performed by matching the whole of the head (including the scalp) to
% the template.  Following this, the registration proceeded by only
% matching the brains together, by appropriate weighting of the template
% voxels.  This is a completely automated procedure (that does not
% require ``scalp editing'') that discounts the confounding effects of
% skull and scalp differences.   A Bayesian framework is used, such that
% the registration searches for the solution that maximizes the a
% posteriori probability of it being correct.  i.e., it maximizes the
% product of the likelihood function (derived from the residual squared
% difference) and the prior function (which is based on the probability
% of obtaining a particular set of zooms and shears).
% 
% The affine registration is followed by estimating nonlinear deformations,
% whereby the deformations are defined by a linear combination of three
% dimensional discrete cosine transform (DCT) basis functions.  The default
% options result in each of the deformation fields being described by 1176
% parameters, where these represent the coefficients of the deformations in
% three orthogonal directions.  The matching involved simultaneously
% minimizing the membrane energies of the deformation fields and the
% residual squared difference between the images and template(s).
% 
% An option is provided for allowing weighting images (consisting of pixel
% values between the range of zero to one) to be used for registering
% abnormal or lesioned brains.  These images should match the dimensions
% of the image from which the parameters are estimated, and should contain
% zeros corresponding to regions of abnormal tissue.
% 
% 
% Uses
% Primarily for stereotactic normalization to facilitate inter-subject
% averaging and precise characterization of functional anatomy.  It is
% not necessary to spatially normalise the data (this is only a
% pre-requisite  for  intersubject averaging or reporting in the
% Talairach space).  If you wish to circumnavigate this step  (e.g. if
% you have single slice data or do not have an appropriate high
% resolution MRI scan) simply specify where you think the  anterior
% commissure  is  with  the  ORIGIN in the header of the first scan
% (using the 'Display' facility) and proceed directly  to  'Smoothing'
% or 'Statistics'.
% 
% 
% Inputs
% The first input is the image which is to be normalised. This image
% should be of the same modality (and MRI sequence etc) as the template
% which is specified. The same spatial transformation can then be
% applied to any other images of the same subject.  These files should
% conform to the SPM data format (See 'Data Format'). Many subjects can
% be entered at once, and there is no restriction on image dimensions
% or voxel size.
% 
% Providing that the images have a correct ".mat" file associated with
% them, which describes the spatial relationship between them, it is
% possible to spatially normalise the images without having first
% resliced them all into the same space. The ".mat" files are generated
% by "spm_realign" or "spm_coregister".
% 
% Default values of parameters pertaining to the extent and sampling of
% the standard space can be changed, including the model or template
% image[s].
% 
% 
% Outputs
% All normalized *.img scans are written to the same subdirectory as
% the original *.img, prefixed with a 'n' (i.e. n*.img).  The details
% of the transformations are displayed in the results window, and the
% parameters are saved in the "*_sn3d.mat" file.
% 
% 
%____________________________________________________________________________
% Refs:
% K.J. Friston, J. Ashburner, C.D. Frith, J.-B. Poline,
% J.D. Heather, and R.S.J. Frackowiak
% Spatial Registration and Normalization of Images.
% Human Brain Mapping 2:165-189(1995)
% 
% J. Ashburner, P. Neelin, D.L. Collins, A.C. Evans and K. J. Friston
% Incorporating Prior Knowledge into Image Registration.
% NeuroImage 6:344-352 (1997)
%
% J. Ashburner and K. J. Friston
% Nonlinear Spatial Normalization using Basis Functions.
% Human Brain Mapping 7(4):in press (1999)
%
%_______________________________________________________________________
%
% --- The Prompts Explained ---
%
% 'Which option?'
% 	'Determine Parameters Only'
% 	Computes the parameters which best fit the image to the
% 	template(s),
% 	and saves them to the file imagename'_sn3d.mat'.
%
% 	'Write Normalised Only'
% 	Select the appropriate '_sn3d.mat' file, and the images you
% 	wish to have normalised. The normalised images will be written
% 	prefixed with 'r'.
%
% 	'Determine Parameters & Write Normalised'
% 	Combines the above two steps into one.
%
% Options for Determine Parameters:
%
% 'select Template(s) '
% Select one or more templates. The image will be fitted (in the least 
% squares sense) to the optimum linear combination of these templates.
%
% 'Normalisation Type?'
% 	'Default Normalisation'
% 	This is a 12 parameter affine normalisation, followed by
% 	5 iterations of the nonlinear normalisation using 4x5x4 basis
% 	functions. If this doesn't work, or if you wish to push the
% 	normalisation a bit harder, try a custom normalisation.
% 	The default arguments for the custom normalisation have
% 	(D) next to them.
%
% 	'Affine Only'
% 	Only perform an affine ('brain in a box') normalisation.
% 	The affine normalisation corrects for gross differences in
% 	brain size, and position.
%
% 	'Custom Affine & Nonlinear'
% 	Perform an affine and also a 3D nonlinear normalisation. The
% 	nonlinear part of the normalisation corrects for more subtle
% 	differences in brain shape.
%
% 'Affine Starting Estimates?'
% The normalised images should use neurological convention for their
% orientation. If any flipping needs doing, it should be done at this
% stage.
%
%	'Neurological Convention (R is R)'
% 	The starting estimates do not include a left right flip.
% 	[0 0 0 0 0 0  1 1 1 0 0 0].
%
%	'Radiological Convention (L is R)'
% 	The starting estimates flip left right (since the template uses
% 	Neurological Convention).
% 	[0 0 0 0 0 0 -1 1 1 0 0 0].
%
% 	'Custom Starting Estimates'
% 	Enter starting estimates in the following order:
% 		x translation (mm)
% 		y translation (mm)
% 		z translation (mm)
% 		x rotation about - {pitch} (radians)
% 		y rotation about - {roll}  (radians)
% 		z rotation about - {yaw}   (radians)
% 		x scaling
% 		y scaling
% 		z scaling
% 		x affine
% 		y affine
% 		z affine
%	For images aquired in different orientations, it is possible
% 	to use starting estimates which reflect these orientations.
% 	eg. a pitch of +/- pi/2 radians for coronal images.
% 	    a roll  of +/- pi/2 radians for saggital images.
% 	For volumes which are flipped, then the appropriate scaling
% 	can be set to -1 in the starting estimates.
%
% '# Nonlinear Iterations?'
% Try about 12 iterations.
%
% '# Basis Functions (x y z)'
% The deformation field is computed from a linear combination of (3D)
% basis images. The basis images used are those which make up the
% lowest frequency components of a 3D discrete cosine transform.
% What is entered here is the dimensions of this transform.
%
% 'Nonlinear Regularization'
% Pick a medium value.  However, if your normalized images appear
% distorted, then it may be an idea to increase the amount of
% regularization - or even just use an affine normalization.
% The regularization influences the smoothness of the deformation
% fields.
%
% 'Mask brain when registering?'
% Applies a weighting mask to the template(s) during the parameter
% estimation.  Weights in and around the brain have values of one
% whereas those clearly outside the brain have values of zero.
% This is an attempt to base the normalization purely upon the
% shape of the brain, rather than the shape of the head (since
% low frequency basis functions can not really cope with variations
% in skull thickness).
%
% 'Mask object brain when registering?'
% Applies a weighting mask to the object image(s) during the parameter
% estimation.  Weights are as for the template mask (0-1).  Used
% (usually) to prevent unusual or abnormal areas of brain (e.g.
% stroke, tumour) influencing normalisation to normal brain
%
% Options for Write Normalised:
% 'select Normalisation Parameter Set'
% Select the '_sn3d.mat' file.
%
% 'Interpolation Method?'
% The method by which the images are sampled when being written in a
% different space.
% 	'Nearest Neighbour'
% 		- Fastest, but not normally recommended.
% 	'Bilinear Interpolation'
% 		- OK for PET, or realigned fMRI.
%
% 	'Sinc Interpolation'
%		- With sinc interpolation, a sinc function with a
%		  Hanning filter envelope is used to resample the data
%		  (11x11x11 kernel).
%
% 'Bounding Box?'
% The bounding box (in mm) of the volume which is to be written
% (relative to the anterior commissure).
%
% 'Voxel Sizes '
% The voxel sizes (x, y & z, in mm) of the written normalised images.
%
%_______________________________________________________________________
% %W% John Ashburner MRCCU/FIL %E%
% With suggested modifications by Matthew Brett of the MRCCU

% Programmers notes
%-----------------------------------------------------------------------
%
% FORMAT spm_sn3d('Defaults')
% acts as a user interface for setting normalisation defaults.
%
%-----------------------------------------------------------------------
%
% FORMAT spm_sn3d(P,matname,bb,Vox,params,spms,brainmask, objmask)
% P         - image(s) to normalize
% matname   - name of file to store deformation definitions
% bb        - bounding box for normalized image
% Vox       - voxel sizes for normalized image
% params(1) - number of basis functions in X
% params(2) - "      "  "     "         "  Y
% params(3) - "      "  "     "         "  Z
% params(4) - number of iterations for elastic deformation
%             Setting any of these parameters to 0 will force
%             the program to perform only the affine
%             normalization.
% params(5) - smoothing for image (mm).
% params(6) - amount of regularization.  Higher values produce smoother deformation
%             fields.
% spms      - template image(s).
% brainmask - Weighting image for template(s)
% objmask   - Weighting image for object images

global SWD sptl_Vx sptl_BB sptl_NBss sptl_Ornt sptl_CO sptl_NItr sptl_Rglrztn;
global sptl_MskBrn sptl_MskObj;

bboxes  = [   -78 78 -112 76  -50 85         
	      -64 64 -104 68  -28 72         
	      -90 91 -126 91  -72 109];
bbprompt =  [' -78:78 -112:76  -50:85  (Default)|'...
	    ' -64:64 -104:68  -28:72  (SPM95)   |'...
	    ' -90:91 -126:91  -72:109 (Template)'];
voxdims    = [ 1   1   1 ; 1.5 1.5 1.5 ; 2   2   2 ; 3   3   3 ; 4   4   4 ; 1   1   2 ; 2   2   4];
voxprompts = ' 1   1   1 | 1.5 1.5 1.5 | 2   2   2 | 3   3   3 | 4   4   4 | 1   1   2 | 2   2   4';
bases =     [0 0 0;2 2 2;2 3 2;3 3 3;3 4 3;4 4 4;4 5 4;5 5 5;5 6 5;6 6 6;6 7 6;7 7 7;7 8 7;8 8 8];
basprompt = 'none|2x2x2|2x3x2|3x3x3|3x4x3|4x4x4|4x5x4|5x5x5|5x6x5|6x6x6|6x7x6|7x7x7|7x8x7|8x8x8';


if (nargin == 0)
	% With no arguments, act as spm_sn3d_ui
	SPMid = spm('FnBanner',mfilename,'%I%');
	[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Normalize');
	spm_help('!ContextHelp','spm_sn3d.m');

	%-----------------------------------------------------------------------
	a1 = spm_input('Which option?',1,'m',...
		'Determine Parameters Only|Write Normalised Only|Determine Parameters & Write Normalised',...
		[1 2 3],3);
	nsubjects = spm_input('# Subjects','+1','w',1, 1);
	if nsubjects < 1,
		spm_figure('Clear','Interactive');
		return;
	end;

	object_masking = sptl_MskObj;
	if (a1 == 1 | a1 == 3) & sptl_CO ~= 1,
		object_masking = spm_input('Mask object brain when registering?', '+1', 'm',...
			'Mask object|Dont Mask object',[1 0],find([1 0] == object_masking));
	end;

	% Select images..
	%-----------------------------------------------------------------------
	for i=1:nsubjects,
		subj(i) = struct('P','','PP','','objmask','','matname','');
		if a1 == 1 | a1 == 3,
			subj(i).P = spm_get(1,'.img',['subj ' num2str(i) ...
				' - Image to determine parameters from']);

			% object weight
 	 		if object_masking==1,
				subj(i).objmask = spm_get(1+sqrt(-1),'.img',...
					['Select object masking image (or Done for none)']);
			end;

			subj(i).matname = [spm_str_manip(subj(i).P,'sd') '_sn3d.mat'];
		else,
			subj(i).matname = spm_get(1,'_sn3d.mat',['subj ' num2str(i) ...
				' - Normalisation parameter set:']);
		end;
	
		if a1 == 2 | a1 == 3,
			subj(i).PP = spm_get(Inf,'.img',['subj ' num2str(i) ...
				' - Images to write normalised']);
		end;
	end;

	if a1 == 1 | a1 == 3,

		% Get template(s)
		ok = 0;
		while ~ok,
			Template = spm_get(Inf,'.img',['Template image(s)'],...
				fullfile(SWD,'templates'));
			vv=spm_vol(Template);
			if prod(size(vv))==1,
				ok = 1;
			elseif prod(size(vv)) ~= 0,
				tmp1 = cat(1,vv.dim);
				tmp2 = reshape(cat(3,vv.mat),4*4,prod(size(vv)));
				if ~any(any(diff(tmp1(:,1:3)))) & ~any(any(diff(tmp2,1,2))),
					ok=1;
				end;
			end;
		end;

		nbasis     = sptl_NBss;
		iterations = sptl_NItr;
		rglrztn    = sptl_Rglrztn;

		brainmask = '';
		if sptl_MskBrn==1, brainmask = fullfile(SWD,'apriori','brainmask.img'); end;

		if sptl_CO ~= 1,
			% Customise the normalisation
			%-----------------------------------------------------------------------
			a2  = spm_input('Normalisation Type?', '+1', 'm',...
				'Default Normalisation|Custom Normalisation', [0 1],1);

			if a2 == 1,

				% Nonlinear options - # basis functions
				%-----------------------------------------------------------------------
				if prod(size(nbasis)) == 3,
					tmp = find(bases(:,1) == nbasis(1) & bases(:,2) == nbasis(2) & ...
						bases(:,3) == nbasis(3));
					if isempty(tmp), tmp = size(bases,1)+1; end
				else,
					tmp = size(bases,1)+2;
				end;
				nb = spm_input('# Nonlinear Basis Functions?','+1','m',...
					[basprompt '|Custom'],[1:size(bases,1) 0], tmp);
				if (nb>0), nbasis = bases(nb,:);
				elseif nb == 0,
					if (prod(size(sptl_NBss)) ~= 3) sptl_NBss = [5 6 5]; end;
					tmp = sprintf('%d %d %d', sptl_NBss(1), sptl_NBss(2), sptl_NBss(3));
					NBss = spm_input('# Basis Functions (x y z)','+0', 'w', tmp, 3);
					NBss = NBss(:)';
					nbasis = NBss(:)';
				else, nbasis = 0; end;

				if prod(nbasis==0),
					iterations = 0;
					rglrztn    = 0;
				else,
					% Nonlinear options - # iterations
					%-----------------------------------------------------------------------
					tmp2 = [1 3 5 8 12 16];
					tmp = find(iterations == tmp2);
					if isempty(tmp) tmp = length(tmp2); end
					iterations = spm_input('# Nonlinear Iterations?','+1','m',...
						['1  nonlinear iteration |3  nonlinear iterations|' ...
						 '5  nonlinear iterations|8  nonlinear iterations|' ...
						 '12 nonlinear iterations|16 nonlinear iterations'],tmp2, tmp);

					% Get the amount of regularization
					%-----------------------------------------------------------------------
					tmp2 = [1 0.1 0.01 0.001 0.0001];
					tmp = find(tmp2 == rglrztn);
					if isempty(tmp) tmp = length(tmp2); end;
					rglrztn = spm_input('Nonlinear Regularization','+1','m',...
						['Extremely Heavy regularization|Heavy regularization|'...
						 'Medium regularization|Light regularization|'...
						 'Very Light regularization'], tmp2, tmp);
				end;

				brainmask = '';
				tmp = spm_input('Mask brain when registering?', '+1', 'm',...
					'Mask Brain|Dont Mask Brain',[1 0],find([1 0] == sptl_MskBrn));
				if tmp == 1, brainmask = fullfile(SWD,'apriori','brainmask.img'); end;
			end;
		end;

	end;

	if a1 == 2 | a1 == 3,

		% Get interpolation method (for writing images)
		%-----------------------------------------------------------------------
		Hold = spm_input('Interpolation Method?','+1','m',...
			['Nearest Neighbour|Bilinear Interpolation|'...
			'Sinc Interpolation (11x11x11)'],[0 1 -11], 2);

		% Get bounding box.
		%-----------------------------------------------------------------------
		if prod(size(sptl_BB)) == 6, bb = sptl_BB;
		else,
			ans = spm_input('Bounding Box?','+1','m',...
				[ bbprompt '|Customise'], [1:size(bboxes,1) 0], 1);
			if ans>0, bb=reshape(bboxes(ans,:),2,3);
			else,
				directions = 'XYZ';
				bb = zeros(2,1);
				for d=1:3,
					str = sprintf('%d %d', bboxes(1,d*2-1), bboxes(1,d*2));
					bb(:,d) = spm_input(['Bounding Box ' directions(d) ],...
						'+1', 'e', str,2);
				end;
			end;
		end


		% Get output voxel sizes.
		%-----------------------------------------------------------------------
		if prod(size(sptl_Vx)) == 3, Vox = sptl_Vx;
		else,
			ans = spm_input('Voxel Sizes?','+1','m',...
				[ voxprompts '|Customise'], [1:size(voxdims,1) 0],3);
			if ans>0, Vox = voxdims(ans,:);
			else, Vox = spm_input('Voxel Sizes ','+0', 'e', '2 2 2',3); end;
		end;
	else,
		bb     = sptl_BB;
		Vox    = sptl_Vx;
	end;

	% Go and do the work
	%-----------------------------------------------------------------------
	spm('Pointer','Watch')
	if a1 == 1 | a1 == 3,
		for i=1:length(subj),
			spm('FigName',['Normalize: working on subj ' num2str(i)],Finter,CmdLine);
			fprintf('\rSubject %d: ', i);
			spm_sn3d(subj(i).P,subj(i).matname,bb,Vox,...
				[nbasis(:)' iterations 8 rglrztn],Template,...
				brainmask, subj(i).objmask);
		end;
	end;
	
	if a1 == 2 | a1 == 3,
		for i=1:length(subj),
			spm('FigName',['Normalize: writing subj ' num2str(i)],Finter,CmdLine);
			fprintf('\rSubject %d: Writing Normalized..', i);
			spm_write_sn(subj(i).PP,subj(i).matname,bb,Vox,Hold);
		end;
	end;
	fprintf('\r%60s%s', ' ',sprintf('\b')*ones(1,60));
	spm('FigName','Normalize: done',Finter,CmdLine);
	spm('Pointer');
	return;

elseif strcmp(P,'Defaults')

	% Edit defaults
	%_______________________________________________________________________


	% Starting estimates for image position
	%-----------------------------------------------------------------------
	if sptl_Ornt == 0, tmp1 = 4;
	elseif sum(sptl_Ornt==[0 0 0 0 0 0  1 1 1 0 0 0])==12, tmp1 = 1;
	elseif sum(sptl_Ornt==[0 0 0 0 0 0 -1 1 1 0 0 0])==12, tmp1 = 2;
	else, tmp1 = 3; end;
	tmp = spm_input(['Affine Starting Estimates?'],1,'m',...
		['Neurological Convention (R is R)|'...
		 'Radiological Convention (L is R)|'...
		 'Custom Affine Starting Estimates'],...
		[0 1 2], tmp1);

	if tmp == 0,     sptl_Ornt = [0 0 0 0 0 0  1 1 1 0 0 0];
	elseif tmp == 1, sptl_Ornt = [0 0 0 0 0 0 -1 1 1 0 0 0];
	elseif tmp == 2,
		se = 0;
		if prod(size(sptl_Ornt)) > 12 | prod(size(sptl_Ornt)) < 8,
			str = '[0 0 0 0 0 0  1 1 1 0 0 0]';
		else,
			str = [num2str(sptl_Ornt(1))];
			for i=2:prod(size(sptl_Ornt))
				str = [str ' ' num2str(sptl_Ornt(i))];
			end
		end;
		sptl_Ornt = spm_input('Enter Affine Starting Estimates:','+1', 'e', str, 12)';
	end;


	% Give option to customise the normalisation
	%-----------------------------------------------------------------------
	tmp = 1;
	if sptl_CO == 1, tmp = 2; end;
	sptl_CO  = spm_input(['Allow customised normalisation?'], '+1', 'm',...
		'   Allow customised|Disallow Customised',[-1 1], tmp);

	% Get number of nonlinear basis functions
	%-----------------------------------------------------------------------
	if prod(size(sptl_NBss)) == 3,
		tmp = find(bases(:,1) == sptl_NBss(1) & bases(:,2) == sptl_NBss(2) & bases(:,3) == sptl_NBss(3));
		if isempty(tmp),tmp = size(bases,1)+1; end;
	else, tmp = size(bases,1)+2; end;

	nb = spm_input(['# Nonlinear Basis Functions?'],'+1','m',[basprompt '|Custom'],[1:size(bases,1) 0], tmp);
	if nb>0, sptl_NBss = bases(nb,:);
	elseif nb == 0,
		if (prod(size(sptl_NBss)) ~= 3) sptl_NBss = [5 6 5]; end;
		tmp = sprintf('%d %d %d', sptl_NBss(1), sptl_NBss(2), sptl_NBss(3));
		NBss = spm_input('# Basis Functions (x y z)','+0', 'w', tmp, 3);
		NBss = NBss(:)';
		sptl_NBss = NBss(:)';
	else, sptl_NBss = 0; end;


	% Get number of nonlinear iterations
	%-----------------------------------------------------------------------
	if prod(sptl_NItr) > 0,
		tmp2 = [1 3 5 8 12 16];
		tmp = find(tmp2 == sptl_NItr);
		if isempty(tmp) tmp = length(tmp2); end;
		sptl_NItr = spm_input(['# Nonlinear Iterations?'],'+1','m',...
			['1  nonlinear iteration |3  nonlinear iterations'...
			'|5  nonlinear iterations|8  nonlinear iterations'...
			'|12 nonlinear iterations|16 nonlinear iterations'],tmp2, tmp);
	else, sptl_NItr = 0; end;

	% Get the amount of regularization
	%-----------------------------------------------------------------------
	tmp2 = [1 0.1 0.01 0.001 0.0001];
	tmp = find(tmp2 == sptl_Rglrztn);
	if isempty(tmp) tmp = length(tmp2); end;
	sptl_Rglrztn = spm_input('Nonlinear Regularization','+1','m',...
		['Extremely Heavy regularization|Heavy regularization|'...
		 'Medium regularization|Light regularization|'...
		 'Very Light regularization'], tmp2, tmp);

	tmp = [1 0];
	sptl_MskBrn = spm_input('Mask brain when registering?', '+1', 'm',...
		'Mask Brain|Dont Mask Brain',[1 0],find(tmp == sptl_MskBrn));
   
	% ask for object image weighting
	sptl_MskObj = spm_input('Mask object brain when registering?', '+1', 'm',...
		'Mask object|Dont Mask object',[1 0],find([1 0] == sptl_MskObj));
   
	% Get default bounding box
	%-----------------------------------------------------------------------
	if prod(size(sptl_BB)) == 6,
		tmp = find(	sptl_BB(1) == bboxes(:,1) & sptl_BB(2) == bboxes(:,2) & ...
				sptl_BB(3) == bboxes(:,3) & sptl_BB(4) == bboxes(:,4) & ...
				sptl_BB(5) == bboxes(:,5) & sptl_BB(6) == bboxes(:,6));
		if isempty(tmp), tmp = size(bboxes,1)+1; end;
	else,
		tmp = size(bboxes,1)+2;
		sptl_BB = reshape(bboxes(1,:),2,3);
	end;

	ans = spm_input('Bounding Box?','+1','m',[ bbprompt '|Customise|Runtime option'], [1:size(bboxes,1) 0 -1], tmp);
	if ans>0, sptl_BB=reshape(bboxes(ans,:),2,3);
	elseif ans == 0,
		if prod(size(sptl_BB)) ~= 6, sptl_BB = reshape(bboxes(1,:),2,3); end;
		directions = 'XYZ';
		bb = zeros(2,1);
		for d=1:3,
			str = sprintf('%d %d', sptl_BB(1,d), sptl_BB(2,d));
			sptl_BB(:,d) = spm_input(['Bounding Box ' directions(d) ], '+1', 'e', str, 2);
		end;
	else, sptl_BB = 0; end;


	% Get default voxel sizes
	%-----------------------------------------------------------------------
	if prod(size(sptl_Vx)) == 3,
		tmp = find(voxdims(:,1) == sptl_Vx(1) & voxdims(:,2) == sptl_Vx(2) & voxdims(:,3) == sptl_Vx(3));
		if isempty(tmp), tmp = size(voxdims,1)+1; end;
	else, tmp = size(voxdims,1)+2; end;
	ans = spm_input(...
		['Voxel Sizes?'], '+1','m', [ voxprompts '|Customise|Runtime option'], [1:size(voxdims,1) 0 -1], tmp);

	if ans>0, sptl_Vx = voxdims(ans,:);
	elseif ans == 0,
		Vox = [];
		if (prod(size(sptl_Vx)) ~= 3) sptl_Vx = [2 2 2]; end
		sptl_Vx = spm_input('Voxel Sizes ','+0', 'e', sprintf('%d %d %d', sptl_Vx(1), sptl_Vx(2), sptl_Vx(3)), 3);
	else, sptl_Vx = 0; end;
	return;
end;



% Perform the spatial normalisation
%_______________________________________________________________________

linfun = inline('fprintf(''  %-60s%s'', x,sprintf(''\b'')*ones(1,62))');

VG = spm_vol(spms);

linfun('Smoothing..');
spm_smooth(P(1,:),fullfile('.','spm_sn3d_tmp.img'),params(5));
VF = spm_vol(fullfile('.','spm_sn3d_tmp.img'));


% Affine Normalisation
%-----------------------------------------------------------------------
spm_chi2_plot('Init','Affine Registration','Convergence');
% use object mask, if present, for rough coregistration
linfun('Coarse Affine Registration..');
p1 = spm_affsub3('affine3',spms,fullfile('.','spm_sn3d_tmp.img'),1,8,[],'', objmask);
% call to affsub3 allows empty brainmask, objmask
linfun('Fine Affine Registration..');
p1 = spm_affsub3('affine3',spms,fullfile('.','spm_sn3d_tmp.img'),1,6,p1,brainmask,objmask);

spm_chi2_plot('Clear');
prms   = p1(1:12);
scales = p1(13:length(p1));
Affine = inv(inv(VG(1).mat)*spm_matrix(prms')*VF(1).mat);

fov = VF(1).dim(1:3).*sqrt(sum(VF(1).mat(1:3,1:3).^2));
if any(fov<60),
	warning('The field of view is too small to attempt nonlinear registration\n');
	params(1:4)=0;
end;

if ~any(params(1:4)==0) & params(6)~=Inf,
	len = fprintf('  %s', '3D Cosine Transform Normalization: ');
	% check for masks and initialise as nec
	VW  = [];
	VW2 = [];
	if nargin >= 8,
		if ~isempty(brainmask), VW=spm_vol(brainmask); end;
		if nargin >= 9,
			if ~isempty(objmask), VW2=spm_vol(objmask); end;
		end;
	end;
   
	[Transform,Dims,scales] = spm_snbasis(VG,VF,Affine,params,VW,VW2);
	fprintf('%s', sprintf('\b')*ones(1,len));
else,
	Transform = [];
	Dims = [VG(1).dim(1:3) ; 0 0 0];
end

% Save parameters for future use.
%-----------------------------------------------------------------------
linfun('Saving Parameters..');
MF     = VF(1).mat;
MG     = VG(1).mat;
vox    = sqrt(sum(MG(1:3,1:3).^2));
origin = MG\[0 0 0 1]';
origin = round(origin(1:3)');
Dims   = [Dims ; vox ; origin ; VF(1).dim(1:3) ; sqrt(sum(MF(1:3,1:3).^2))];
mgc    = 960209;
eval(['save ' matname ' mgc Affine Dims Transform MF MG -v4']);

delete_image(fullfile('.','spm_sn3d_tmp.img'));

% Do the display stuff
%-----------------------------------------------------------------------
fg = spm_figure('FindWin','Graphics');
if ~isempty(fg),
	spm_figure('Clear','Graphics');
	ax=axes('Position',[0.1 0.51 0.8 0.45],'Visible','off','Parent',fg);
	text(0,0.90, 'Spatial Normalisation','FontSize',16,'FontWeight','Bold',...
		'Interpreter','none','Parent',ax);
	text(0,0.75, [ 'Image		: ' P(1,:)],'FontSize',12,'FontWeight','Bold',...
		'Interpreter','none','Parent',ax);
	text(0,0.7, [ 'Parameters	: ' spm_str_manip(matname,'sd')],'FontSize',12,...
		'Interpreter','none','Parent',ax);

	str = 'no flipping'; 
	if det(Affine(1:3,1:3))<0, str = 'image flipped'; end;
	Q = spm_matrix(prms');
	text(0,0.6, ['Linear {affine} component - ' str],'FontWeight','Bold',...
		'Interpreter','none','Parent',ax);
	text(0,0.55, sprintf('X1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(1,:)),...
		'Interpreter','none','Parent',ax);
	text(0,0.50, sprintf('Y1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(2,:)),...
		'Interpreter','none','Parent',ax);
	text(0,0.45, sprintf('Z1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(3,:)),...
		'Interpreter','none','Parent',ax);

	if (~any(params(1:4)==0) & params(6)~=Inf),
		text(0,0.35, sprintf('%d nonlinear iterations',params(4)),...
			'Interpreter','none','Parent',ax);
		text(0,0.30, sprintf('%d x %d x %d basis functions',params(1:3)),...
			'Interpreter','none','Parent',ax);
	else,
		text(0,0.35, 'No nonlinear components',...
			'Interpreter','none','Parent',ax);
	end;

	h1=spm_orthviews('Image',deblank(spms(1,:)),[0. 0.1 .5 .5]);
	linfun('Writing Image for Display..');
	spm_write_sn(P(1,:),matname,bb,Vox,1);
	p  = spm_str_manip(P(1,:), 'd');
	q  = max([find(p == spm_platform('sepchar')) 0]);
	q  = [p(1:q) 'n' p((q + 1):length(p))];
	h2=spm_orthviews('Image',q,[.5 0.1 .5 .5]);
	spm_orthviews('Space',h2);
	spm_print;
	drawnow;
end;
linfun(' ');
return;
%_______________________________________________________________________
%_______________________________________________________________________
function delete_image(iname)
iname = spm_str_manip(iname,'sd');
spm_unlink([iname '.img'], [iname '.hdr'], [iname '.mat'], [iname '.mnc']);
return;

