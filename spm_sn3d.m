function spm_sn3d(P,matname,bb,Vox,params,spms,aff_parms,brainmask,objmask)
% Perform 3D spatial normalization
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
%
% The routine attempts to match the input image to an
% optimum linear combination of template images. This
% provides additional flexibility in the type of input images
% which can be used. Typical template images consist of:
%     gray matter
%     white matter
%     scalp
%     striatum
%     ventricles
%     etc...
%
% First of all, an affine normalization stage is used to estimate
% the overall size, orientation etc of the image.
% Following this, an elastic deformation is computed which
% will match the image to the template.
% The deformation is in 3-dimensions, and is constrained
% to consist of a linear combination of basis functions.
% The basis functions chosen were the lowest frequency
% components of a 3 dimensional discrete cosine transform.
%
% The parameters for the affine transformation, and 3D
% basis function transformation are saved.
%
% for a complete description of this approach see Friston et al (1994)
% ref: Friston et al (1994) The spatial registration and normalization
% of images HBM 0:00-00
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
% FORMAT spm_sn3d(P,matname,bb,Vox,params,spms,aff_parms,brainmask, objmask)
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
% aff_parms - starting parameters for affine normalisation.
% brainmask - Weighting image for template(s)
% objmask   - Weighting image for object images

global SWD CWD sptl_Vx sptl_BB sptl_NBss sptl_Ornt sptl_CO sptl_NItr sptl_Rglrztn;
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

	spm_figure('Clear','Interactive');
	set(spm_figure('FindWin','Interactive'),'Name','Spatial Normalisation')
	spm_help('!ContextHelp','spm_sn3d.m');

	pos = 1;

	%-----------------------------------------------------------------------
	a1 = spm_input('Which option?',pos,'m',...
		'Determine Parameters Only|Write Normalised Only|Determine Parameters & Write Normalised',...
		[1 2 3],3);
	pos = pos + 1;
	nsubjects = spm_input('# Subjects',pos,'e',1);
	if (nsubjects < 1)
		spm_figure('Clear','Interactive');
		return;
	end

	pos = pos + 1;

	% Select images..
	%-----------------------------------------------------------------------
	for i=1:nsubjects
		if (a1 == 1 | a1 == 3)
			P = spm_get(1,'.img',['subj ' num2str(i) ' - Image to determine parameters from']);
			eval(['P' num2str(i) '=P;']);
			matname = [spm_str_manip(P,'sd') '_sn3d.mat'];
		else
			matname = spm_get(1,'_sn3d.mat',['subj ' num2str(i) ' - Normalisation parameter set:']);
		end
		eval(['matname' num2str(i) '=matname;']);
	
		if (a1 == 2 | a1 == 3)
			P = spm_get(Inf,'.img',['subj ' num2str(i) ' - Images to write normalised']);
			eval(['PP' num2str(i) '=P;']);
		end
	end

	aff_parms = [0 0 0 0 0 0 1 1 1 0 0 0];

	if (a1 == 1 | a1 == 3)

		if sptl_MskBrn==1,
			brainmask = [SWD '/apriori/brainmask.img'];
		else,
			brainmask = '';
		end;
		% object weight
  		if sptl_MskObj==1,
			objmask = spm_get(1,'.img',['Object masking image']);
		else,
			objmask = '';
		end;

		nbasis     = sptl_NBss;
		iterations = sptl_NItr;
		rglrztn = sptl_Rglrztn;
		aff_parms = sptl_Ornt;


		% Get template(s)
		ok = 0;
		while ~ok,
			Template = spm_get(Inf,'.img',['Template image(s)'], [SWD '/templates']);
			vv=spm_vol(Template);
			if prod(size(vv))==1,
				ok = 1;
			else,
				tmp1 = cat(1,vv.dim);
				tmp2 = reshape(cat(3,vv.mat),4*4,prod(size(vv)));
				if ~any(any(diff(tmp1(:,1:3)))) & ~any(any(diff(tmp2,1,2))),
					ok=1;
				end;
			end;
		end;

		if sptl_CO ~= 1
			% Customise the normalisation
			%-----------------------------------------------------------------------
			a2  = spm_input('Normalisation Type?', pos, 'm',...
				'Default Normalisation|Affine Only|Custom Affine & Nonlinear', [0 1 2],1);
			pos = pos + 1;

			if (a2 == 2)
				tmp2 = [1 3 5 8 12 16];
				tmp = find(iterations == tmp2);
				if isempty(tmp) tmp = length(tmp2); end

				% Nonlinear options - # iterations
				%-----------------------------------------------------------------------
				tmp2 = [1 3 5 8 12 16];
				tmp = find(tmp2 == sptl_NItr);
				if isempty(tmp) tmp = length(tmp2); end;
				iterations = spm_input('# Nonlinear Iterations?',pos,'m','1  nonlinear iteration |3  nonlinear iterations|5  nonlinear iterations|8  nonlinear iterations|12 nonlinear iterations|16 nonlinear iterations',tmp2, tmp);
				pos = pos + 1;

				% Nonlinear options - # basis functions
				%-----------------------------------------------------------------------
				nbasis = [];
				if (prod(size(sptl_NBss)) == 3)
					tmp = find(bases(:,1) == sptl_NBss(1) & bases(:,2) == sptl_NBss(2) & bases(:,3) == sptl_NBss(3));
					if isempty(tmp)
						tmp = size(bases,1)+1;
					end
				else
					tmp = size(bases,1)+2;
				end

				nb = spm_input('# Nonlinear Basis Functions?',pos,'m',[basprompt '|Custom'],[1:size(bases,1) 0], tmp);
				if (nb>0), nbasis = bases(nb,:); end
				while prod(size(nbasis)) ~= 3 | any(nbasis < 1) | prod(nbasis) > 1000
					tmp = sprintf('%d %d %d', sptl_NBss(1), sptl_NBss(2), sptl_NBss(3));
					nbasis = spm_input('# Basis Functions (x y z)',pos, 'e', tmp);
					nbasis = nbasis(:)';
				end
				pos = pos+1;

				% Get the amount of regularization
				%-----------------------------------------------------------------------
				tmp2 = [1 0.1 0.01 0.001 0.0001];
				tmp = find(tmp2 == sptl_Rglrztn);
				if isempty(tmp) tmp = length(tmp2); end;
					rglrztn = spm_input('Nonlinear Regularization',pos,'m',...
					['Extremely Heavy regularization|Heavy regularization|'...
					 'Medium regularization|Light regularization|'...
					 'Very Light regularization'], tmp2, tmp);
				pos = pos+1;

			elseif (a2 == 1)
				nbasis     = [0 0 0];
				iterations = 0;
				rglrztn    = 0;
			end

			if a2==1 | a2==2,
				tmp = [1 0];
				tmp = spm_input('Mask brain when registering?', pos, 'm',...
					'Mask Brain|Dont Mask Brain',[1 0],find(tmp == sptl_MskBrn));
				pos = pos + 1;
				if tmp == 1,
					brainmask = [SWD '/apriori/brainmask.img'];
				else,
					brainmask = '';
				end;
            
				% get object weighting 
				tmp = [1 0];
				tmp = spm_input('Mask object brain when registering?', pos, 'm',...
					'Mask object|Dont Mask object',[1 0],find(tmp == sptl_MskObj));
				pos = pos + 1;
				if tmp == 1,
					objmask = spm_get(1,'.img',['Object masking image']);
				else,
					objmask = '';
				end;
			end;
		end;


		% Affine starting estimate
		%-----------------------------------------------------------------------
		if prod(size(sptl_Ornt))<12
			tmp = spm_input('Affine Starting Estimates?',pos,'m',...
				['Neurological Convention (R is R)|'...
				 'Radiological Convention (L is R)|'...
				 'Custom Affine Starting Estimates'],...
				[0 1 2]);
			pos = pos + 1;

			if (tmp == 0)
				aff_parms = [0 0 0 0 0 0  1 1 1 0 0 0];
			elseif (tmp == 1)
				aff_parms = [0 0 0 0 0 0 -1 1 1 0 0 0];
			elseif (tmp == 2)
				aff_parms = [0 0 0 0 0 0  1 1 1 0 0 0];
				se = 0;
				while prod(size(se)) > 12 | prod(size(se)) < 8
					se = spm_input('Enter Affine Starting Estimates:',pos);
					se = se(:)';
				end
				aff_parms(1:prod(size(se))) = se(:);
			end
		end
	end

	if (a1 == 2 | a1 == 3)

		% Get interpolation method (for writing images)
		%-----------------------------------------------------------------------
		Hold = spm_input('Interpolation Method?',pos,'m',['Nearest Neighbour|Bilinear Interpolation|'...
			'Sinc Interpolation (11x11x11)'],[0 1 -11], 3);
		pos = pos + 1;


		% Get bounding box.
		%-----------------------------------------------------------------------
		if prod(size(sptl_BB)) ~= 6
			ans = spm_input('Bounding Box?',pos,'m',[ bbprompt '|Customise'], [1:size(bboxes,1) 0], 1);
			if (ans>0)
				pos = pos + 1;
				bb=reshape(bboxes(ans,:),2,3);
			else
				directions = 'XYZ';
				bb = zeros(2,1);
				for d=1:3
					bbx = [];
					while size(bbx,1) ~= 2
						str = sprintf('%d %d', bboxes(1,d*2-1), bboxes(1,d*2));
						bbx = spm_input(['Bounding Box ' directions(d) ], pos, 'e', str);
						bbx = bbx(:);
					end
					bb(:,d) = bbx;
					pos = pos + 1;
				end
			end
		else
			bb = sptl_BB;
		end


		% Get output voxel sizes.
		%-----------------------------------------------------------------------
		if prod(size(sptl_Vx)) ~= 3
			ans = spm_input('Voxel Sizes?',pos,'m',[ voxprompts '|Customise'], [1:size(voxdims,1) 0],3);
			if (ans>0)
				Vox = voxdims(ans,:);
			else
				Vox = [];
				while size(Vox,2) ~= 3
					Vox = spm_input('Voxel Sizes ',pos, 'e', '2 2 2');
					Vox = Vox(:)';
				end
			end
			pos = pos + 1;
		else
			Vox = sptl_Vx;
		end

	else
		bb     = sptl_BB;
		Vox    = sptl_Vx;
	end


	% Go and do the work
	%-----------------------------------------------------------------------
	set(spm_figure('FindWin','Interactive'),'Name','Normalising','Pointer','Watch'); drawnow;
	if (a1 == 1 | a1 == 3)
		for i=1:nsubjects
			eval(['matname=matname' num2str(i) ';']);
			eval(['P=P' num2str(i) ';']);
			spm_sn3d(P,matname,bb,Vox,[nbasis iterations 8 rglrztn],Template,aff_parms,...
			brainmask, objmask);
		end
	end
	set(spm_figure('FindWin','Interactive'),'Name','Writing     Normalised','Pointer','Watch');
	drawnow;
	
	if (a1 == 2 | a1 == 3)
		for i=1:nsubjects
			eval(['matname=matname' num2str(i) ';']);
			eval(['P=PP' num2str(i) ';']);
			spm_write_sn(P,matname,bb,Vox,Hold);
		end
	end
	
	spm_figure('Clear','Interactive');
	return;

elseif strcmp(P,'Defaults')

	% Edit defaults
	%_______________________________________________________________________

	pos = 2;

	% Starting estimates for image position
	%-----------------------------------------------------------------------
	if sptl_Ornt == 0
		tmp1 = 4;
	elseif sum(sptl_Ornt==[0 0 0 0 0 0  1 1 1 0 0 0])==12
		tmp1 = 1;
	elseif sum(sptl_Ornt==[0 0 0 0 0 0 -1 1 1 0 0 0])==12
		tmp1 = 2;
	else
		tmp1 = 3;
	end
	tmp = spm_input(['Affine Starting Estimates?'],pos,'m',...
		['Neurological Convention (R is R)|'...
		 'Radiological Convention (L is R)|'...
		 'Custom Affine Starting Estimates|'...
		 'Runtime option'],...
		[0 1 2 -1], tmp1);
	pos = pos + 1;

	if (tmp == 0)
		sptl_Ornt = [0 0 0 0 0 0  1 1 1 0 0 0];
	elseif (tmp == 1)
		sptl_Ornt = [0 0 0 0 0 0 -1 1 1 0 0 0];
	elseif (tmp == 2)
		se = 0;

		if prod(size(sptl_Ornt)) > 12 | prod(size(sptl_Ornt)) < 8
			str = '[0 0 0 0 0 0  1 1 1 0 0 0]';
		else
			str = [num2str(sptl_Ornt(1))];
			for i=2:prod(size(sptl_Ornt))
				str = [str ' ' num2str(sptl_Ornt(i))];
			end
		end
		while prod(size(se)) > 12 | prod(size(se)) < 8
			se = spm_input('Enter Affine Starting Estimates:',pos, 'e', str);
			se = se(:)';
		end
		sptl_Ornt = [0 0 0 0 0 0  1 1 1 0 0 0];
		sptl_Ornt(1:prod(size(se))) = se(:);
	elseif (tmp == -1)
		sptl_Ornt = 0;
	end


	% Give option to customise the normalisation
	%-----------------------------------------------------------------------
	tmp = 1;
	if sptl_CO == 1, tmp = 2; end;
	sptl_CO  = spm_input(['Allow customised normalisation?'], pos, 'm',...
		'   Allow customised|Disallow Customised',[-1 1], tmp);
	pos = pos + 1;

	% Get number of nonlinear basis functions
	%-----------------------------------------------------------------------
	if (prod(size(sptl_NBss)) == 3)
		tmp = find(bases(:,1) == sptl_NBss(1) & bases(:,2) == sptl_NBss(2) & bases(:,3) == sptl_NBss(3));
		if isempty(tmp)
			tmp = size(bases,1)+1;
		end
	else
		tmp = size(bases,1)+2;
	end

	nb = spm_input(['# Nonlinear Basis Functions?'],pos,'m',[basprompt '|Custom'],[1:size(bases,1) 0], tmp);
	if (nb>0)
		sptl_NBss = bases(nb,:);
	elseif nb == 0
		if (prod(size(sptl_NBss)) ~= 3) sptl_NBss = [5 6 5]; end;
		NBss = [];
		while prod(size(NBss)) ~= 3 | any(NBss < 0) | prod(NBss) > 1000
			tmp = sprintf('%d %d %d', sptl_NBss(1), sptl_NBss(2), sptl_NBss(3));
			NBss = spm_input('# Basis Functions (x y z)',pos, 'e', tmp);
			NBss = NBss(:)';
		end
		sptl_NBss = NBss(:)';
	else
		sptl_NBss = 0;
	end
	pos = pos+1;


	% Get number of nonlinear iterations
	%-----------------------------------------------------------------------
	if prod(sptl_NItr) > 0
		tmp2 = [1 3 5 8 12 16];
		tmp = find(tmp2 == sptl_NItr);
		if isempty(tmp) tmp = length(tmp2); end;
		sptl_NItr = spm_input(['# Nonlinear Iterations?'],pos,'m',...
			['1  nonlinear iteration |3  nonlinear iterations'...
			'|5  nonlinear iterations|8  nonlinear iterations'...
			'|12 nonlinear iterations|16 nonlinear iterations'],tmp2, tmp);
		pos = pos + 1;
	else
		sptl_NItr = 0;
	end

	% Get the amount of regularization
	%-----------------------------------------------------------------------
	tmp2 = [1 0.1 0.01 0.001 0.0001];
	tmp = find(tmp2 == sptl_Rglrztn);
	if isempty(tmp) tmp = length(tmp2); end;
	sptl_Rglrztn = spm_input('Nonlinear Regularization',pos,'m',...
		['Extremely Heavy regularization|Heavy regularization|'...
		 'Medium regularization|Light regularization|'...
		 'Very Light regularization'], tmp2, tmp);
	pos = pos+1;

	tmp = [1 0];
	sptl_MskBrn = spm_input('Mask brain when registering?', pos, 'm',...
		'Mask Brain|Dont Mask Brain',[1 0],find(tmp == sptl_MskBrn));
	pos = pos + 1;
   
   % ask for object image weighting
   tmp = [1 0];
	sptl_MskObj = spm_input('Mask object brain when registering?', pos, 'm',...
		'Mask object|Dont Mask object',[1 0],find(tmp == sptl_MskObj));
   pos = pos + 1;
   
   % Get default bounding box
	%-----------------------------------------------------------------------
	if prod(size(sptl_BB)) == 6
		tmp = find(	sptl_BB(1) == bboxes(:,1) & sptl_BB(2) == bboxes(:,2) & ...
				sptl_BB(3) == bboxes(:,3) & sptl_BB(4) == bboxes(:,4) & ...
				sptl_BB(5) == bboxes(:,5) & sptl_BB(6) == bboxes(:,6));
		if isempty(tmp) tmp = size(bboxes,1)+1; end;
	else
		tmp = size(bboxes,1)+2;
		sptl_BB = reshape(bboxes(1,:),2,3);
	end

	ans = spm_input('Bounding Box?',pos,'m',[ bbprompt '|Customise|Runtime option'], [1:size(bboxes,1) 0 -1], tmp);
	if (ans>0)
		pos = pos + 1;
		sptl_BB=reshape(bboxes(ans,:),2,3);
	elseif (ans == 0)
		if prod(size(sptl_BB)) ~= 6, sptl_BB = reshape(bboxes(1,:),2,3); end;
		directions = 'XYZ';
		bb = zeros(2,1);
		for d=1:3
			bbx = [];
			while size(bbx,1) ~= 2
				str = sprintf('%d %d', sptl_BB(1,d), sptl_BB(2,d));
				bbx = spm_input(['Bounding Box ' directions(d) ], pos, 'e', str);
				bbx = bbx(:);
			end
			sptl_BB(:,d) = bbx;
			pos = pos + 1;
		end
	else
		sptl_BB = 0;
		pos = pos + 1;
	end


	% Get default voxel sizes
	%-----------------------------------------------------------------------
	if (prod(size(sptl_Vx)) == 3)
		tmp = find(voxdims(:,1) == sptl_Vx(1) & voxdims(:,2) == sptl_Vx(2) & voxdims(:,3) == sptl_Vx(3));
		if isempty(tmp)
			tmp = size(voxdims,1)+1;
		end
	else
		tmp = size(voxdims,1)+2;
	end
	ans = spm_input(...
		['Voxel Sizes?'], pos,'m', [ voxprompts '|Customise|Runtime option'], [1:size(voxdims,1) 0 -1], tmp);

	if (ans>0)
		sptl_Vx = voxdims(ans,:);
	elseif (ans == 0)
		Vox = [];
		if (prod(size(sptl_Vx)) ~= 3) sptl_Vx = [2 2 2]; end
		while size(Vox,2) ~= 3
			Vox = spm_input('Voxel Sizes ',pos, 'e', sprintf('%d %d %d', sptl_Vx(1), sptl_Vx(2), sptl_Vx(3)));
		end
		sptl_Vx = Vox(:)';
	else
		sptl_Vx = 0;
	end
	pos = pos + 1;

	return;
end



% Perform the spatial normalisation
%_______________________________________________________________________

VG = spm_vol(spms);

fprintf('smoothing..');
spm_smooth(P(1,:),'./spm_sn3d_tmp.img',params(5));
VF = spm_vol('./spm_sn3d_tmp.img');
fprintf(' ..done\n');


% Affine Normalisation
%-----------------------------------------------------------------------
spm_chi2_plot('Init','Affine Registration','Convergence');
% use object mask, if present, for rough coregistration
p1 = spm_affsub3('affine3',spms,'./spm_sn3d_tmp.img',1,8,[],'', objmask);
% call to affsub3 allows empty brainmask, objmask
p1 = spm_affsub3('affine3',spms,'./spm_sn3d_tmp.img',1,6,p1,brainmask,objmask);

spm_chi2_plot('Clear');
prms   = p1(1:12);
scales = p1(13:length(p1));
Affine = inv(inv(VG(1).mat)*spm_matrix(prms')*VF(1).mat);

fov = VF(1).dim(1:3).*sqrt(sum(VF(1).mat(1:3,1:3).^2));
if any(fov<60),
	fprintf('The field of view is too small to attempt nonlinear registration\n');
	params(1:4)=0;
end;

if ~any(params(1:4)==0) & params(6)~=Inf,
	fprintf('3D Cosine Transform Normalization\n');
	% check for masks and initialise as nec
	VW  = [];
	VW2 = [];
	if nargin >= 8,
		if ~isempty(brainmask),
			VW=spm_vol(brainmask);
		end;
		if nargin >= 9,
			if ~isempty(objmask),
				VW2=spm_vol(objmask);
			end;
		end;
	end;
   
	[Transform,Dims,scales] = spm_snbasis(VG,VF,Affine,params,VW,VW2);
else,
	Transform = [];
	Dims = [VG(1).dim(1:3) ; 0 0 0];
end

% Save parameters for future use.
%-----------------------------------------------------------------------
MF     = VF(1).mat;
MG     = VG(1).mat;
vox    = sqrt(sum(MG(1:3,1:3).^2));
origin = MG\[0 0 0 1]';
origin = origin(1:3)';
Dims   = [Dims ; vox ; origin ; VF(1).dim(1:3) ; sqrt(sum(MF(1:3,1:3).^2))];
mgc    = 960209;
eval(['save ' matname ' mgc Affine Dims Transform MF MG -v4']);

delete_image('./spm_sn3d_tmp.img');

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

	Q = spm_matrix(prms');
	text(0,0.6, 'Linear {affine} component','FontWeight','Bold',...
		'Interpreter','none','Parent',ax);
	text(0,0.55, sprintf('X1 = %0.2f*X + %0.2f*Y + %0.2f*Z + %0.2f',Q(1,:)),...
		'Interpreter','none','Parent',ax);
	text(0,0.50, sprintf('Y1 = %0.2f*X + %0.2f*Y + %0.2f*Z + %0.2f',Q(2,:)),...
		'Interpreter','none','Parent',ax);
	text(0,0.45, sprintf('Z1 = %0.2f*X + %0.2f*Y + %0.2f*Z + %0.2f',Q(3,:)),...
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
	spm_write_sn(P(1,:),matname,bb,Vox,1);
	p  = spm_str_manip(P(1,:), 'd');
	q  = max([find(p == '/') 0]);
	q  = [p(1:q) 'n' p((q + 1):length(p))];
	h2=spm_orthviews('Image',q,[.5 0.1 .5 .5]);
	spm_orthviews('Space',h2);
	spm_print;
	drawnow;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function delete_image(iname)
iname = spm_str_manip(iname,'sd');
spm_unlink([iname '.img'], [iname '.hdr'], [iname '.mat'], [iname '.mnc']);
return;

