function spm_sn3d(P,matname,bb,Vox,params,spms,aff_parms,free)
% Perform 3D stereotactic Normalization
%
% FORMAT spm_sn3d
% With no arguments, spm_sn3d acts as its own user interface.
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
% 	'Custom Affine Only'
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
% '# Affine Params?'
% Most good data sets should be stable for 9 or 12 parameters.
%
%
% '# Nonlinear Iterations?'
% Try 5 iterations.
%
% '# Basis Functions (x y z)'
% The deformation field is computed from a linear combination of (3D)
% basis images. The basis images used are those which make up the
% lowest frequency components of a 3D discrete cosine transform.
% What is entered here is the dimensions of this transform.
% Try 4 5 4.
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
%	With sinc interpolation, a sinc function with a Hanning filter
% 	envelope is used to resample the data.
% 	'Rough   Sinc Interpolation'            -  7x 7x 7 kernel.
% 	'Medium  Sinc Interpolation (slow)'     - 11x11x11 kernel.
% 	'Better Sinc Interpolation (very slow)' - 17x17x17 kernel.
%
% 'Bounding Box?'
% The bounding box (in mm) of the volume which is to be written
% (relative to the anterior commissure).
%
% 'Voxel Sizes '
% The voxel sizes (x, y & z, in mm) of the written normalised images.
%
%-----------------------------------------------------------------------
%
% FORMAT spm_sn3d('Defaults')
% acts as a user interface for setting normalisation defaults.
%
%-----------------------------------------------------------------------
%
% FORMAT spm_sn3d(P,matname,bb,Vox,params,spms,aff_parms,aff_free)
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
% params(6) - 'smoothness' of deformation field.
% spms      - template image(s).
% aff_parms - starting parameters for affine normalisation.
% aff_free  - a vector of ones and zeros representing which parameters
%             for the affine component should remain fixed at the
%             starting estimates.
%             - ones for parameters to be fitted.
%             - zeros for fixed parameters.
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

% %W% John Ashburner MRCCU/FIL %E%


global SWD CWD sptl_Vx sptl_BB sptl_NBss sptl_NAP sptl_Ornt sptl_CO sptl_NItr;


bboxes  = [   -78 78 -112 76  -50 85         
	      -64 64 -104 68  -28 72         
	      -90 91 -126 91  -72 109];
bbprompt =  [' -78:78 -112:76  -50:85  (MNI)     |'...
	    ' -64:64 -104:68  -28:72  (SPM95)   |'...
	    ' -90:91 -126:91  -72:109 (Template)'];
voxdims    = [ 1   1   1 ; 1.5 1.5 1.5 ; 2   2   2 ; 3   3   3 ; 4   4   4 ; 1   1   2 ; 2   2   4];
voxprompts = ' 1   1   1 | 1.5 1.5 1.5 | 2   2   2 | 3   3   3 | 4   4   4 | 1   1   2 | 2   2   4';
bases =     [2 2 2;2 3 2;3 3 3;3 4 3;4 4 4;4 5 4;5 5 5;5 6 5;6 6 6;6 7 6;7 7 7;7 8 7;8 8 8];
basprompt = '2x2x2|2x3x2|3x3x3|3x4x3|4x4x4|4x5x4|5x5x5|5x6x5|6x6x6|6x7x6|7x7x7|7x8x7|8x8x8';


if (nargin == 0)
	% With no arguments, act as spm_sn3d_ui

	spm_figure('Clear','Interactive');
	pos = 1;

	%-----------------------------------------------------------------------
	a1 = spm_input('Which option?',pos,'m',...
		'Determine Parameters Only|Write Normalised Only|Determine Parameters & Write Normalised',...
		[1 2 3]);
	pos = pos + 1;
	nsubjects = spm_input('# Subjects',pos);
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
		% Get template(s)
		ok = 0;
		while (~ok)
			Template = spm_get(Inf,'.img',['Template image(s)'],'', [SWD '/templates']);
			if (size(Template,1)>0)
				dims = zeros(size(Template,1),9);
				for i=1:size(Template,1)
					[dim vox dummy dummy dummy origin dummy] = spm_hread(deblank(Template(i,:)));
					dims(i,:) = [dim vox origin];
				end
				if size(dims,1) == 1 | ~any(diff(dims))
					ok = 1;
				end
			end
		end

		if sptl_CO ~= 1
			% Customise the normalisation
			%-----------------------------------------------------------------------
			a2  = spm_input('Normalisation Type?', pos, 'm',...
				'Default Normalisation|Custom Affine Only|Custom Affine & Nonlinear', [0 1 2]);
			pos = pos + 1;
			if (a2 == 0)
				nap        = sptl_NAP;
				iterations = sptl_NItr;
				nbasis     = sptl_NBss;
			end


			% Number of affine parameters
			%-----------------------------------------------------------------------
			if (a2 == 1 | a2 == 2)
				nap = spm_input('# Affine Params?',pos,'m',[
					' 2 params (X & Y Translations)|'...
				' 3 params (Translations)|'...
				' 5 params (Translations + Pitch & Roll)|'...
				' 6 params (Rigid Body)|'...
				' 8 params (Rigid Body + X & Y Zoom)|'...
				' 9 params (Rigid Body + Zooms)|'...
				'11 params (Rigid Body,  Zooms + X & Y Affine)|'...
				'12 params (Rigid Body,  Zooms & Affine)'
				],[2 3 5 6 8 9 11 12]);
				pos = pos + 1;
			end

			if (a2 == 2)
				% Nonlinear options - # iterations
				%-----------------------------------------------------------------------
				iterations = spm_input('# Nonlinear Iterations?',pos,'m','1  nonlinear iteration |3  nonlinear iterations|5  nonlinear iterations|8  nonlinear iterations|12 nonlinear iterations|16 nonlinear iterations',[1 3 5 8 12 16]);
				pos = pos + 1;

				% Nonlinear options - # basis functions
				%-----------------------------------------------------------------------
				nbasis = [];
				nb = spm_input('# Nonlinear Basis Functions?',pos,'m',[basprompt '|Custom'],[1:size(bases,1) 0]);
				if (nb>0), nbasis = bases(nb,:); end
				while prod(size(nbasis)) ~= 3 | any(nbasis < 1) | prod(nbasis) > 1000
					nbasis = spm_input('# Basis Functions (x y z)',pos);
					nbasis = nbasis(:)';
				end
				pos = pos+1;
			elseif (a2 == 1)
				nbasis     = [0 0 0];
				iterations = 0;
				smoothness = 0;
			end
		else
			nbasis     = sptl_NBss;
			iterations = sptl_NItr;
			nap        = sptl_NAP;
		end
	

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
				while prod(size(se)) > 12 & prod(size(se)) < 8
					se = spm_input('Enter Affine Starting Estimates:',pos);
					se = se(:)';
				end
				aff_parms(1:prod(size(se))) = se(:);
			end
		else
			aff_parms = sptl_Ornt;
		end

	end

	if (a1 == 2 | a1 == 3)

		% Get interpolation method (for writing images)
		%-----------------------------------------------------------------------
		Hold = spm_input('Interpolation Method?',pos,'m',['Nearest Neighbour|Bilinear Interpolation|'...
			'Rough   Sinc Interpolation|Medium  Sinc Interpolation (slow)|'...
			'Better Sinc Interpolation (very slow)'],[0 1 3 5 8]);
		pos = pos + 1;


		% Get bounding box.
		%-----------------------------------------------------------------------
		if prod(size(sptl_BB)) ~= 6
			ans = spm_input('Bounding Box?',pos,'m',[ bbprompt '|Customise'], [1:size(bboxes,1) 0]);
			if (ans>0)
				pos = pos + 1;
				bb=reshape(bboxes(ans,:),2,3);
			else
				directions = 'XYZ';
				bb = zeros(2,1);
				for d=1:3
					bbx = [];
					while size(bbx,1) ~= 2
						bbx = spm_input(['Bounding Box ' directions(d) ], pos);
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
			ans = spm_input('Voxel Sizes?',pos,'m',[ voxprompts '|Customise'], [1:size(voxdims,1) 0]);
			if (ans>0)
				Vox = voxdims(ans,:);
			else
				Vox = [];
				while size(Vox,2) ~= 3
					Vox = spm_input('Voxel Sizes ',pos);
					Vox = Vox(:)';
				end
			end
			pos = pos + 1;
		else
			Vox = sptl_Vx;
		end

	else
		bb     = reshape(bboxes(1,:),2,3);
		Vox    = [2 2 4];
	end


	% Go and do the work
	%-----------------------------------------------------------------------
	set(spm_figure('FindWin','Interactive'),'Name','Normalising','Pointer','Watch'); drawnow;
	if (a1 == 1 | a1 == 3)
		for i=1:nsubjects
			eval(['matname=matname' num2str(i) ';']);
			eval(['P=P' num2str(i) ';']);
			spm_sn3d(P,matname,bb,Vox,[nbasis iterations 8 0.01],Template,aff_parms,ones(nap,1));
		end
	end
	set(spm_figure('FindWin','Interactive'),'Name','Writing     Normalised','Pointer','Watch'); drawnow;
	
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
		tmp1 = 'runtime';
	elseif sum(sptl_Ornt==[0 0 0 0 0 0  1 1 1 0 0 0])==12
		tmp1 = 'R is R';
	elseif sum(sptl_Ornt==[0 0 0 0 0 0 -1 1 1 0 0 0])==12
		tmp1 = 'L is R';
	else
		tmp1 = 'custom';
	end
	tmp = spm_input(['Affine Starting Estimates (' tmp1 ')?'],pos,'m',...
		['Neurological Convention (R is R)|'...
		 'Radiological Convention (L is R)|'...
		 'Custom Affine Starting Estimates|'...
		 'Runtime option'],...
		[0 1 2 -1]);
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
			str = ['[' num2str(sptl_Ornt(1))];
			for i=2:prod(size(sptl_Ornt))
				str = [str ' ' num2str(sptl_Ornt(i))]
			end
			str = [str ']']
		end
		while prod(size(se)) > 12 | prod(size(se)) < 8
			se = spm_input('Enter Affine Starting Estimates:',pos, 'e', str);
			se = se(:)';
		end
		sptl_Ornt = [0 0 0 0 0 0  1 1 1 0 0 0];
		aff_parms(1:prod(size(se))) = se(:);
	elseif (tmp == -1)
		sptl_Ornt = 0;
	end
	pos = pos + 1;


	% Give option to customise the normalisation
	%-----------------------------------------------------------------------
	tmp = 'allow';
	if sptl_CO == 1, tmp = 'disallow'; end;
	sptl_CO  = spm_input(['Allow customised normalisation (' tmp ')?'], pos, 'm',...
		'   Allow customised|Disallow Customised',[-1 1]);
	pos = pos + 1;


	% Get number of affine parameters
	%-----------------------------------------------------------------------
	sptl_NAP = spm_input(['# Affine Params (' num2str(sptl_NAP) ')?'],pos,'m',[
		' 2 params (X & Y Translations)|'...
		' 3 params (Translations)|'...
		' 5 params (Translations + Pitch & Roll)|'...
		' 6 params (Rigid Body)|'...
		' 8 params (Rigid Body + X & Y Zoom)|'...
		' 9 params (Rigid Body + Zooms)|'...
		'11 params (Rigid Body,  Zooms + X & Y Affine)|'...
		'12 params (Rigid Body,  Zooms & Affine)'
		],[2 3 5 6 8 9 11 12]);
	pos = pos + 1;


	% Get number of nonlinear basis functions
	%-----------------------------------------------------------------------
	nb = spm_input(['# Nonlinear Basis Functions (' num2str(sptl_NBss(1)) 'x'  num2str(sptl_NBss(2)) 'x'  num2str(sptl_NBss(3)) ')?'],pos,'m',[basprompt '|Custom'],[1:size(bases,1) 0]);
	if (nb>0)
		sptl_NBss = bases(nb,:);
	elseif nb == 0
		sptl_NBss = [];
		while prod(size(sptl_NBss)) ~= 3 | any(sptl_NBss < 0) | prod(sptl_NBss) > 1000
			sptl_NBss = spm_input('# Basis Functions (x y z)',pos);
			sptl_NBss = sptl_NBss(:)';
		end
	else
		sptl_NBss = 0;
	end
	pos = pos+1;


	% Get number of nonlinear iterations
	%-----------------------------------------------------------------------
	if prod(sptl_NBss) > 0
		sptl_NItr = spm_input(['# Nonlinear Iterations (' num2str(sptl_NItr) ')?'],pos,'m',...
			['1  nonlinear iteration |3  nonlinear iterations'...
			'|5  nonlinear iterations|8  nonlinear iterations'...
			'|12 nonlinear iterations|16 nonlinear iterations'],[1 3 5 8 12 16]);
		pos = pos + 1;
	else
		sptl_NItr = 0;
	end


	% Get default bounding box
	%-----------------------------------------------------------------------
	ans = spm_input('Bounding Box?',pos,'m',[ bbprompt '|Customise|Runtime option'], [1:size(bboxes,1) 0 -1]);
	if (ans>0)
		pos = pos + 1;
		sptl_BB=reshape(bboxes(ans,:),2,3);
	elseif (ans == 0)
		directions = 'XYZ';
		bb = zeros(2,1);
		for d=1:3
			bbx = [];
			while size(bbx,1) ~= 2
				bbx = spm_input(['Bounding Box ' directions(d) ], pos);
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
	if prod(size(sptl_Vx))~=3
		str = 'runtime';
	else
		str = [num2str(sptl_Vx(1)) ' '  num2str(sptl_Vx(2)) ' '  num2str(sptl_Vx(3))];
	end
	ans = spm_input(...
		['Voxel Sizes (' str ')?'], pos,'m', [ voxprompts '|Customise|Runtime option'], [1:size(voxdims,1) 0 -1]);

	if (ans>0)
		sptl_Vx = voxdims(ans,:);
	elseif (ans == 0)
		Vox = [];
		while size(Vox,2) ~= 3
			Vox = spm_input('Voxel Sizes ',pos);
			sptl_Vx = Vox(:)';
		end
	else
		sptl_Vx = 0;
	end
	pos = pos + 1;

	return;
end



% Perform the spatial normalisation
%_______________________________________________________________________



% Map the template(s).
%-----------------------------------------------------------------------
VG = [];
for i=1:size(spms,1)
	fname = spms(i,:);
	fname = fname(fname ~= ' ');
	VG = [VG spm_map(fname)];
	if (i==1)
		[dim,vox,scale,type,offset,origin,descrip] = spm_hread(fname);
	end
end

% Check for consistancy
%-----------------------------------------------------------------------
if (size(VG,2)> 1 & any(diff(VG(1:6,:)')))
	error('Templates must have identical dimensions');
end

% Map the image to normalize.
%-----------------------------------------------------------------------
fprintf('smoothing..');
spm_smooth(P,'./spm_sn3d_tmp.img',params(5));
VF = spm_map('./spm_sn3d_tmp.img');
fprintf(' ..done\n');

MF = spm_get_space(P(1,:));
MG = spm_get_space(spms(1,:));

[prms,scales] = spm_affine(VG,VF,MG,MF,aff_parms,free,1,0.25);
Affine = inv(inv(MG)*spm_matrix(prms')*MF);

if (~any(params==0))
	fprintf('3D Cosine Transform Normalization\n');
	[Transform,Dims,scales] = spm_snbasis_map(VG,VF,Affine,params);
else
	Transform = [];
	Dims = [VG(1:3,1)' ; 0 0 0];
end

% Save parameters for future use.
%-----------------------------------------------------------------------
Dims = [Dims ; vox ; origin ; VF(1:3,1)' ; VF(4:6,1)'];
mgc = 960209;
eval(['save ' matname ' mgc Affine Dims Transform MF MG']);

spm_unmap_vol(VF);

spm_write_sn('./spm_sn3d_tmp.img',matname,bb,Vox,1);
spm_unlink ./spm_sn3d_tmp.img ./spm_sn3d_tmp.hdr ./spm_sn3d_tmp.mat
[tmp1,tmp2,tmp3,tmp4,tmp5,origin2,tmp6] = spm_hread('./nspm_sn3d_tmp.img');

% Do the display stuff
%-----------------------------------------------------------------------
centre = [10 0 0];
VGT = VG;
VGT(7,:) = VGT(7,:).*scales(:)';
figure(spm_figure('FindWin','Graphics'));
spm_figure('Clear','Graphics');
axes('Position',[0.1 0.51 0.8 0.45],'Visible','off');
text(0,0.90, 'Spatial Normalisation','FontSize',16,'FontWeight','Bold');
text(0,0.75, [ 'Image		: ' P(1,:)],'FontSize',12,'FontWeight','Bold');
text(0,0.7, [ 'Parameters	: ' [spm_str_manip(P,'sd') '_sn3d.mat']],'FontSize',12);

Q = spm_matrix(prms');
text(0,0.6, 'Linear {affine} component','FontWeight','Bold');
text(0,0.55, sprintf('X1 = %0.2f*X + %0.2f*Y + %0.2f*Z + %0.2f',Q(1,:)));
text(0,0.50, sprintf('Y1 = %0.2f*X + %0.2f*Y + %0.2f*Z + %0.2f',Q(2,:)));
text(0,0.45, sprintf('Z1 = %0.2f*X + %0.2f*Y + %0.2f*Z + %0.2f',Q(3,:)));

if (~any(params==0))
	text(0,0.35, sprintf('%d nonlinear iterations',params(4)));
	text(0,0.30, sprintf('%d x %d x %d basis functions',params(1:3)));
else
	text(0,0.35, 'No nonlinear components');
end

VN = spm_map('./nspm_sn3d_tmp.img');
spm_orthviews(VGT,centre,origin ,bb,1,[0. 0.1 .5 .5],'Template');
spm_orthviews(VN ,centre,origin2,bb,1,[.5 0.1 .5 .5],'Normalized');
drawnow;
spm_print

spm_unlink ./nspm_sn3d_tmp.img ./nspm_sn3d_tmp.hdr ./nspm_sn3d_tmp.mat

for v=[VG VN]
	spm_unmap_vol(v);
end

