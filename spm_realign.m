function spm_realign(arg1,arg2,arg3,arg4,arg5,arg6)
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
% %W% Karl Friston & John Ashburner - with input from Oliver Josephs %E%


global MODALITY sptl_WhchPtn sptl_DjstFMRI sptl_CrtWht sptl_MskOptn SWD

if (nargin == 0)
	% User interface.
	%_______________________________________________________________________
	SPMid = spm('FnBanner',mfilename,'%I%');
	[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Realign');
	spm_help('!ContextHelp','spm_realign.m');

	pos = 1;

	n     = spm_input('number of subjects', pos, 'e', 1);
	if (n < 1)
		spm_figure('Clear','Interactive');
		return;
	end

	P = cell(n,1);
	sessions = cell(n,1);

	pos = pos + 1;
	for i = 1:n
		if strcmp(MODALITY,'FMRI'),
			ns = spm_input(['num sessions for subject 'num2str(i)], pos, 'e', 1);
			sess = zeros(1,ns);
			pp = [];
			for s=1:ns
				p = '';
				while size(p,1)<1,
					p = spm_get(Inf,'.img',...
						['scans for subj ' num2str(i) ', sess' num2str(s)]);
				end;
				if s==1,
					pp = p;
				else
					pp=str2mat(pp,p);
				end
				sess(s) =  size(p,1);
			end
			sessions{i} = cumsum(sess);
			P{i} = pp;
		else
			ns = 1;
			p  = '';
			while size(p,1)<1,
				p = spm_get(Inf,'.img',...
					['select scans for subject ' num2str(i)]);
			end;
			P{i} = p;
			sessions{i} = size(P{i},1);
		end
	end

	FlagsC = ' ';
	FlagsR = ' ';

	if sptl_WhchPtn == 1
		WhchPtn = 3;
	else
		WhchPtn = spm_input('Which option?', pos, 'm',...
			'Coregister only|Reslice Only|Coregister & Reslice',...
			[1 2 3],3);
		pos = pos + 1;
	end

	% Co-registration options
	%-----------------------------------------------------------------------
	if WhchPtn == 1 | WhchPtn == 3,
		% 'Mask brain when registering?'
		% The object is to match brains together using a rigid body transformation.
		% It is possible that nonrigid movement of structures outside the brain
		% (eyeballs, toung, cheeks etc) may interfere with this registration.
		% To get around this problem, non-brain areas of the image can be masked
		% out during the parameter estimation step.  In order to acheive this, a
		% template (of the appropriate modality) is first registered with the images.
		% This allows a brain/non-brain mask to be overlayed on the images during
		% the registration process.
		%
		Q = '';
		tmp = 0;
		%tmp = spm_input('Mask brain when registering?', pos, 'm',...
		%	'Mask Brain|Dont Mask Brain',...
		%	[1 0],2);
		if tmp==1
			if (strcmp(MODALITY, 'FMRI')),
				def = 1;
			else
				def = 2;
			end
			DIR1 = fullfile(SWD,'coreg');
			templates = str2mat(	fullfile(DIR1,'EPI.img'),...
						fullfile(DIR1,'PET.img'),...
						fullfile(DIR1,'T1.img'),...
						fullfile(DIR1,'T2.img'),...
						fullfile(DIR1,'Transm.img'));
			tmp = spm_input('Modality of images?',3,'m',...
				'EPI MR images|PET images|T1 MR images|T2 MR images|Transm images',...
				[1 2 3 4 5],def);
			Q = str2mat(deblank(templates(tmp,:)),...
				fullfile(SWD,'apriori','brainmask.img'));
		end
		pos = pos + 1;
	end

	% Reslicing options
	%-----------------------------------------------------------------------
	if (WhchPtn == 2 | WhchPtn == 3),
		p = spm_input('Reslice interpolation method?',pos,'m',...
			'Trilinear Interpolation|Sinc Interpolation|Fourier space Interpolation',...
			[1 2 3],2);

		pos = pos + 1;
		if (p == 2) FlagsR = [FlagsR 'S']; end
		if (p == 3) FlagsR = [FlagsR 'F']; end

		if sptl_CrtWht == 1
			p = 3;
		else
			p = spm_input('Create what?',pos,'m',...
				[' All Images (1..n)| Images 2..n|'...
				 ' All Images + Mean Image| Mean Image Only'],...
				[1 2 3 4],3);
			pos = pos + 1;
		end
		if (p == 2) FlagsR = [FlagsR 'n']; end
		if (p == 3) FlagsR = [FlagsR 'i']; end
		if (p == 4) FlagsR = [FlagsR 'Ni']; end
		if (~any(FlagsR == 'N'))
			if sptl_MskOptn == 1
				FlagsR = [FlagsR 'k'];
			else
				if (spm_input('Mask the resliced images?',pos,'y/n') == 'y')
					FlagsR = [FlagsR 'k'];
				end
				pos = pos + 1;
			end
			if (strcmp(MODALITY, 'FMRI'))
				FlagsR = [FlagsR 'a'];

				if ~any(FlagsR == 'F'),
					if (sptl_DjstFMRI == 1)
						FlagsR = [FlagsR 'c'];
					elseif sptl_DjstFMRI ~= 0
						if (spm_input(...
							'Adjust sampling errors?'...
							,pos,'y/n') == 'y')
							FlagsR = [FlagsR 'c'];
						end
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
			realign_images(P{i},Q,sessions{i});
		end
		if WhchPtn==2 | WhchPtn==3,
			reslice_images(P{i},FlagsR,sessions{i})
		end;
	end
	fprintf('\r%60s%s', ' ',sprintf('\b')*ones(1,60));
	spm('FigName','Realign: done',Finter,CmdLine);
	spm('Pointer');
	return;

elseif nargin == 1 & strcmp(arg1,'Defaults'),
	edit_defaults;
	return;
elseif nargin >= 1 & strcmp(arg1,'Reslice'),
	P = arg2;
	Flags = 'n';
	Sessions = size(P,1);
	if nargin>=3,
		Flags = arg3;
		if nargin>=4,
			Sessions = arg4;
		end;
	end;
	reslice_images(P,Flags,Sessions);
elseif nargin >= 1 & strcmp(arg1,'Realign'),
	P = arg2;
	Q = '';
	if nargin<3, sessions = length(P);
	else, sessions = arg3; end;
	realign_images(P,Q,sessions);
end;
%_______________________________________________________________________

%_______________________________________________________________________
function reslice_images(P,Flags,sessions)
% FORMAT reslice_images(P,Flags,sessions)
%
% P     - matrix of filenames {one string per row}
%         All operations are performed relative to the first image.
%         ie. Coregistration is to the first image, and resampling
%         of images is into the space of the first image.
%
% sessions - the last scan in each of the sessions.  For example,
%            the images in the second session would be
%            P((sessions(2-1)+1):sessions(2),:).
%
% Flags - options flags
%
%         c - adjust the data (fMRI) to remove movement-related components.
%
%         k - mask output images
%             To avoid artifactual movement-related variance the realigned
%             set of images can be internally masked, within the set (i.e.
%             if any image has a zero value at a voxel than all images have
%             zero values at that voxel).  Zero values occur when regions
%             'outside' the image are moved 'inside' the image during
%             realignment.
%
%         i - write mean image
%             The average of all the realigned scans is written to
%             mean*.img.
%
%         S - use sinc interpolation for reslicing (9x9x9).
%
%         n - don't reslice the first image
%             The first image is not actually moved, so it may not be
%             necessary to resample it.
%
%         N - don't reslice any of the images - except possibly create a
%             mean image.
%
%         a - write absolute values in images - more appropriate for fMRI.
%
%             The spatially realigned and adjusted images are written to
%             the orginal subdirectory with the same filename but prefixed
%             with a 'r'.
%             They are all aligned with the first.

if any(Flags == 'F') | ~any(Flags == 'c'),
	if any(Flags == 'c'), warning('No adjustment will be done'); end;
	reslice_images_volbyvol(P,Flags,sessions);
	return;
end;

take_abs = 0; if any(Flags == 'a'), take_abs = 1; end;

P=spm_vol(P);

linfun = inline('fprintf(''  %-60s%s'', x,sprintf(''\b'')*ones(1,60))');
linfun('Reslicing images..');

Hold = 1;
if (any(Flags == 'S')) Hold = -9; end

start_vol = 1;
if (any(Flags == 'n') & ~any(Flags == 'i')) start_vol = 2; end

% Write headers and matrixes
%------------------------------------------------------------------
if ~any(Flags == 'N'),
	PO = P;
	for i = start_vol:prod(size(P)),
		PO(i).fname   = prepend(P(i).fname,'r');
		PO(i).dim     = [P(1).dim(1:3) P(i).dim(4)];
		PO(i).mat     = P(1).mat;
		PO(i).descrip = 'spm - realigned';
		spm_create_image(PO(i));
	end
end

spm_progress_bar('Init',P(1).dim(3),'Reslicing','planes completed');

if any(Flags == 'i')
	Integral = zeros(P(1).dim(1)*P(1).dim(2),P(1).dim(3));
end

tiny = 5e-2; % From spm_vol_utils.c

if any(Flags == 'c'),
	linfun('Estimating a priori covariance matrix..');
	Y1 = zeros(prod(P(1).dim(1:2)),prod(size(P)));
	Y2 = zeros(prod(P(1).dim(1:2)),prod(size(P)));
	Y3 = zeros(prod(P(1).dim(1:2)),prod(size(P)));

	% Estimate a priori covariance matrix describing the distribution of
	% the basis functions for the adjustment.
	gl = spm_global(P(1));
	varerr = 0;
	npix   = 0;
	for x3 = 1:P(1).dim(3),
		tmp = spm_slice_vol(P(1),spm_matrix([0 0 x3]),P(1).dim(1:2),0);
		msk = find(tmp > gl*0.8);
		tmp = reshape(interp_errors(tmp,Hold),prod(P(1).dim(1:2)),4);
		varerr = varerr + sum(sum(tmp(msk,:).^2));
		npix = npix + prod(size(msk));
	end;
	varerr = varerr/npix;
	IC0  = eye(6)/varerr;
end;

X  = zeros(prod(P(1).dim(1:2)),prod(size(P)));


x1=repmat((1:P(1).dim(1))',1,P(1).dim(2));x1=x1(:);
x2=repmat( 1:P(1).dim(2)  ,P(1).dim(1),1);x2=x2(:);

for x3 = 1:P(1).dim(3)
	linfun(['Reslicing plane ' num2str(x3) '..']);
	Count = zeros(prod(P(1).dim(1:2)),1);
	for i = 1:prod(size(P))
		M = inv(P(1).mat\P(i).mat);
		y1=M(1,1)*x1+M(1,2)*x2+(M(1,3)*x3+M(1,4));
		y2=M(2,1)*x1+M(2,2)*x2+(M(2,3)*x3+M(2,4));
		y3=M(3,1)*x1+M(3,2)*x2+(M(3,3)*x3+M(3,4));

		Mask =        (y1 >= (1-tiny) & y1 <= (P(i).dim(1)+tiny));
		Mask = Mask & (y2 >= (1-tiny) & y2 <= (P(i).dim(2)+tiny));
		Mask = Mask & (y3 >= (1-tiny) & y3 <= (P(i).dim(3)+tiny));
		Count = Count + Mask;

		d  = spm_sample_vol(P(i),y1,y2,y3,Hold);
		if take_abs==1, d = abs(d); end;
		X(:,i)= d;

		if any(Flags == 'c'),
			Y1(:,i)=y1;
			Y2(:,i)=y2;
			Y3(:,i)=y3;
		end;
	end;
	Mask = (Count == prod(size(P)));

	if any(Flags == 'c'),
		ss = 1;
		for s = 1:length(sessions),
			se = sessions(s);
			degrf = se-ss+1-7;

			for xx=find(Mask)',
				% Basis functions appropriate for sinc interpolation
				A = [	cos(2*pi*Y1(xx,ss:se)') sin(2*pi*Y1(xx,ss:se)') ...
			  		cos(2*pi*Y2(xx,ss:se)') sin(2*pi*Y2(xx,ss:se)') ...
			   		cos(2*pi*Y3(xx,ss:se)') sin(2*pi*Y3(xx,ss:se)')];

				Tac   = X(xx,ss:se)';
				% centre the signal and covariates
				Tac   = Tac - mean(Tac);
				A     = A - repmat(mean(A,1),size(A,1),1);
				alpha = A'*A;
				beta  = A'*Tac;

				% Estimate the distribution of the errors.
				sumsq = sum((Tac - A*(pinv(alpha)*beta)).^2)/degrf;

				% Subtract a MAP estimate of the interpolation errors
				% from the data.
				Tac   = X(xx,ss:se)' - A*((alpha + IC0*sumsq)\beta);
				X(xx,ss:se) = Tac';
			end
			ss = se+1;
		end;
	end;


	if any(Flags == 'i'),
		Integral(:,x3) = sum(X,2)./Count;
	end;

	if ~any(Flags == 'N'),
		if any(Flags == 'k'), notmsk = find(~Mask); else, notmsk=[]; end;
		start_vol = 1;
		if (any(Flags == 'n')) start_vol = 2; end

		for i = start_vol:prod(size(P)),
			tmp = reshape(X(:,i),PO(i).dim(1:2));
			tmp(notmsk) = NaN;
			spm_write_plane(PO(i),tmp,x3);
		end;
	end
	spm_progress_bar('Set',x3);
end;


if any(Flags == 'i'),
	% Write integral image (16 bit signed)
	%-----------------------------------------------------------

	PO         = P(1);
	PO.fname   = prepend(P(1).fname, 'mean');
	PO.pinfo   = [max(max(Integral))/32767 0 0]';
	PO.descrip = 'spm - mean image';
	PO.dim(4)  = 4;
	spm_write_vol(PO,reshape(Integral,PO.dim(1:3)));
end;

spm_figure('Clear','Interactive');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function realign_images(P,Q,sessions)
% FORMAT realign_images(P,Q,sessions)
%
% P     - matrix of filenames {one string per row}
%         All operations are performed relative to the first image.
%         ie. Coregistration is to the first image, and resampling
%         of images is into the space of the first image.
%
% sessions - the last scan in each of the sessions.  For example,
%            the images in the second session would be
%            P((sessions(2-1)+1):sessions(2),:).
%
% Q     - an optional matrix of filenames.  These are used for masking
%         out regions which are (roughly) considered to be outside the
%         brain.  The last filename is an image containing values
%         between zero and one, where each value is a weight.  The
%         other images are template images.  The way that this works is
%         that the first image in the series is matched (using an affine
%         transformation) to a linear combinantion of the template images.
%         This affine mapping can then be used to overlay the weight image
%         over the first image of the series.
%

linfun = inline('fprintf(''  %-60s%s'', x,sprintf(''\b'')*ones(1,60))');

global MODALITY
P=spm_vol(P);
Flags = ' ';
if strcmp(MODALITY, 'PET'), Flags = [Flags 'p']; end;
if nargin<2, Q = ''; end
if ~isempty(Q),
	tmp = clock;
	ref = [pwd '/spm_realign_tmp_' sprintf('%.2d%.2d%.2d%.2d%.4d',tmp(3),tmp(4),tmp(5),round(tmp(6)),round(rem(tmp(6),1)*10000)) '.img'];
	spm_smooth(P(1).fname,ref,8);
	linfun('Initial registration to templates..');
	params  =spm_affsub3('affine2',Q(1:(end-1),:),ref,1,8);
	params = spm_affsub3('affine2',Q(1:(end-1),:),ref,1,6,params);
	tmp = spm_str_manip(ref,'rd');
	spm_unlink(ref,[tmp '.hdr'],[tmp '.mat']);
	MW = inv(spm_matrix(params(1:12)));
	PW = spm_vol(Q(end,:));
	PW.mat = inv(spm_matrix(params(1:12)))*PW.mat;
else
	PW=[];
end

if (length(sessions)==1),
	linfun('Registering images..');
	P = realign_series(P,PW,Flags);
	save_parameters(P);
else
	linfun('Registering together the first image of each session..');
	tmp = [1 sessions(1:(end-1))+1];
	Ptmp = realign_series(P(tmp),PW,Flags);
	ss=1;
	for s=1:length(sessions)
		es = sessions(s);
		M = Ptmp(s).mat*inv(P(ss).mat);
		for i=ss:es
			P(i).mat = M*P(i).mat;
		end;
		ss = es+1;
	end;

	ss=1;
	for s=1:length(sessions)
		es = sessions(s);
		linfun(['Registering together images from session ' num2str(s) '..']);
		P(ss:es) = realign_series(P(ss:es),PW,Flags);
		save_parameters(P(ss:es));
		ss = es+1;
	end;
end

% Save Realignment Parameters
%---------------------------------------------------------------------------
linfun('Saving parameters..');
for i=1:prod(size(P)),
	spm_get_space(P(i).fname, P(i).mat);
end;

plot_parameters(P);

return;
%_______________________________________________________________________

%_______________________________________________________________________
function P = realign_series(P,Wt,Flags)
% Realign a time series of 3D images to the first of the series.
% FORMAT P = realign_series(P,Wt,Flags)
% P  - a vector of volumes (see spm_vol)
% Wt - an optional masking volume
%-----------------------------------------------------------------------
% P(i).mat is modified to reflect the modified position of the image i.
% The scaling (and offset) parameters are also set to contain the
% optimum scaling required to match the images.
%_______________________________________________________________________

if nargin < 2, Wt = []; end;
if nargin < 3, Flags = ' '; end;

Hold1 = -8;	% Interpolation method
fwhm  = 6;	% Smoothing kernel size
separation = 4.5;	% Spacing (in millimeters) between samples
register_to_mean = 0;

% because PET images are noisier:
if any(Flags == 'p'), register_to_mean = 1; fwhm = 8; end;

if any(Flags == 'd'), separation = 2.5; end;
if any(Flags == 's'), fwhm = 10; end;

lkp = [1 2 3 4 5 6];
if P(1).dim(3) < 3, lkp = [1 2 6]; end;

% Points to sample in reference image
%-----------------------------------------------------------------------
skip = sqrt(sum(P(1).mat(1:3,1:3).^2)).^(-1)*separation;
d    = P(1).dim(1:3);
[x1,x2,x3]=meshgrid(1:skip(1):d(1),1:skip(2):d(2),1:skip(3):d(3));
x1 = x1(:);
x2 = x2(:);
x3 = x3(:);


% Possibly mask an area of the sample volume.
%-----------------------------------------------------------------------
if ~isempty(Wt),
	[y1,y2,y3]=coords([0 0 0  0 0 0],P(1).mat,Wt.mat,x1,x2,x3);
	wt = spm_sample_vol(Wt,y1,y2,y3,1);
	msk = find(wt>0.01);
	x1=x1(msk);
	x2=x2(msk);
	x3=x3(msk);
	wt=wt(msk);
else
	wt = [];
end
n=prod(size(x1));


% Compute rate of change of chi2 w.r.t changes in parameters (matrix A)
% Also compute A'A for the whole volume, so that A'A for the masked
% volume can be computed more efficiently. 
%-----------------------------------------------------------------------
V=smooth_vol(P(1),fwhm);
[b,dG1,dG2,dG3]=spm_sample_vol(V,x1,x2,x3,Hold1);
clear V
A0 = make_A(P(1).mat,x1,x2,x3,dG1,dG2,dG3,wt,lkp);
if ~isempty(wt), b = b.*wt; end

if register_to_mean,
	count = ones(n,1);
	ave   = b;
	grad1 = dG1;
	grad2 = dG2;
	grad3 = dG3;
end;

spm_progress_bar('Init',length(P)-1,'Registering Images');
% Loop over images
%-----------------------------------------------------------------------
for i=2:length(P),
	V=smooth_vol(P(i),fwhm);
	ss = Inf;
	countdown = -1;
	Hold = 1;	% Begin with bi-linear interpolation.
	for iter=1:64,
		[y1,y2,y3] = coords([0 0 0  0 0 0],P(1).mat,P(i).mat,x1,x2,x3);
		msk        = find((y1>=1 & y1<=d(1) & y2>=1 & y2<=d(2) & y3>=1 & y3<=d(3)));
		if length(msk)<32, error_message(P(i)); end;

		F          = spm_sample_vol(V, y1(msk),y2(msk),y3(msk),Hold);
		if ~isempty(wt), F = F.*wt(msk); end;

		A          = [A0(msk,:) F];
		Alpha      = spm_atranspa(A);
		Beta       = A'*b(msk);
		soln       = Alpha\Beta;

		p          = [0 0 0  0 0 0  1 1 1  0 0 0];
		p(lkp)     = soln(1:(end-1));
		P(i).mat   = inv(spm_matrix(p))*P(i).mat;

		pss        = ss;
		ss         = sum((F*soln(end)-b(msk)).^2)/length(msk);
		if (pss-ss)/pss < 1e-8,	% Stopped converging.
			if Hold == 1
				% Switch to a better (slower) interpolation
				% to finish with
				Hold = Hold1;
			elseif countdown == -1
				% Do two final iterations to finish off with
				countdown = 2;
			end
		end;
		if countdown ~= -1,
			if countdown==0, break; end;
			countdown = countdown -1;
		end;
	end;

	if register_to_mean,
		% Generate mean and derivatives of mean
		tiny = 5e-2; % From spm_vol_utils.c
		msk        = find((y1>=(1-tiny) & y1<=(d(1)+tiny) &...
		                   y2>=(1-tiny) & y2<=(d(2)+tiny) &...
		                   y3>=(1-tiny) & y3<=(d(3)+tiny)));
		count(msk) = count(msk) + 1;
		[G,dG1,dG2,dG3] = spm_sample_vol(V,y1(msk),y2(msk),y3(msk),Hold1);
		ave(msk)   = ave(msk)   + G.*soln(end);
		grad1(msk) = grad1(msk) + dG1.*soln(end);
		grad2(msk) = grad2(msk) + dG2.*soln(end);
		grad3(msk) = grad3(msk) + dG3.*soln(end);
	end;
	spm_progress_bar('Set',i-1);
end;
spm_progress_bar('Clear');

if ~register_to_mean, return; end;
%_______________________________________________________________________
M=P(1).mat;
A0 = make_A(M,x1,x2,x3,grad1./count,grad2./count,grad3./count,wt,lkp);
if ~isempty(wt), b = (ave./count).*wt;
else, b = (ave./count); end

clear ave grad1 grad2 grad3

% Loop over images
%-----------------------------------------------------------------------
spm_progress_bar('Init',length(P),'Registering Images to Mean');
for i=1:length(P),
	V=smooth_vol(P(i),fwhm);
	ss = Inf;
	countdown = -1;
	for iter=1:64,
		[y1,y2,y3] = coords([0 0 0  0 0 0],M,P(i).mat,x1,x2,x3);
		msk        = find((y1>=1 & y1<=d(1) & y2>=1 & y2<=d(2) & y3>=1 & y3<=d(3)));
		if length(msk)<32, error_message(P(i)); end;

		F          = spm_sample_vol(V, y1(msk),y2(msk),y3(msk),Hold1);
		if ~isempty(wt), F = F.*wt(msk); end;

		A          = [A0(msk,:) F];
		Alpha      = spm_atranspa(A);
		Beta       = A'*b(msk);
		soln       = Alpha\Beta;
		p          = [0 0 0  0 0 0  1 1 1  0 0 0];
		p(lkp)     = soln(1:(end-1));
		P(i).mat   = inv(spm_matrix(p))*P(i).mat;

		pss        = ss;
		ss         = sum((F*soln(end)-b(msk)).^2)/length(msk);
		if (pss-ss)/pss < 1e-8 & countdown == -1 % Stopped converging.
			% Do three final iterations to finish off with
			countdown = 3;
		end;
		if countdown ~= -1
			if countdown==0, break; end;
			countdown = countdown -1;
		end;
	end;
	spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear');


% Since we are supposed to be aligning everything to the first
% image, then we had better do so
%-----------------------------------------------------------------------
M = M/P(1).mat;
for i=1:length(P)
	P(i).mat   = M*P(i).mat;
end

return;
%_______________________________________________________________________

%_______________________________________________________________________
function [y1,y2,y3]=coords(p,M1,M2,x1,x2,x3)
% Rigid body transformation of a set of coordinates.
M  = (inv(M2)*inv(spm_matrix(p(1:6)))*M1);
y1 = M(1,1)*x1 + M(1,2)*x2 + M(1,3)*x3 + M(1,4);
y2 = M(2,1)*x1 + M(2,2)*x2 + M(2,3)*x3 + M(2,4);
y3 = M(3,1)*x1 + M(3,2)*x2 + M(3,3)*x3 + M(3,4);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function V=smooth_vol(P,fwhm)
% Convolve the volume in memory.
s = sqrt(sum(P.mat(1:3,1:3).^2)).^(-1)*(fwhm/sqrt(8*log(2)));
x  = round(6*s(1)); x = [-x:x];
y  = round(6*s(2)); y = [-y:y];
z  = round(6*s(3)); z = [-z:z];
x  = exp(-(x).^2/(2*(s(1)).^2));
y  = exp(-(y).^2/(2*(s(2)).^2));
z  = exp(-(z).^2/(2*(s(3)).^2));
x  = x/sum(x);
y  = y/sum(y);
z  = z/sum(z);

i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;

V=zeros(P.dim(1:3));
spm_conv_vol(P,V,x,y,z,-[i j k]);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function A = make_A(M,x1,x2,x3,dG1,dG2,dG3,wt,lkp)
p0 = [0 0 0  0 0 0  1 1 1  0 0 0];
At=zeros(length(lkp),prod(size(x1)));
for i=1:length(lkp)
	pt         = p0;
	pt(lkp(i)) = pt(i)+1e-6;
	[y1,y2,y3] = coords(pt,M,M,x1,x2,x3);
	tmp        = sum([y1-x1 y2-x2 y3-x3].*[dG1 dG2 dG3],2)/(-1e-6);
	if ~isempty(wt), A(:,i) = tmp.*wt;
	else, A(:,i) = tmp; end
end
return;
%_______________________________________________________________________

%_______________________________________________________________________
function error_message(P)

str = {	'There is not enough overlap in the images',...
	'to obtain a solution.',...
	' ',...
	'Offending image:',...
	 P.fname,...
	' ',...
	'Please check that your header information is OK.'};
spm('alert*',str,mfilename,sqrt(-1));
error('insufficient image overlap')

return
%_______________________________________________________________________

%_______________________________________________________________________
function plot_parameters(P)
fg=spm_figure('FindWin','Graphics');
if ~isempty(fg)

	Params = zeros(prod(size(P)),12);
	for i=1:prod(size(P))
		Params(i,:) = spm_imatrix(P(i).mat/P(1).mat);
	end

	% display results
	% translation and rotation over time series
	%-------------------------------------------------------------------
	spm_figure('Clear','Graphics');
	ax=axes('Position',[0.1 0.65 0.8 0.2],'Parent',fg,'Visible','off');
	set(get(ax,'Title'),'String','Image realignment','FontSize',16,'FontWeight','Bold','Visible','on');
	x     =  0.1;
	y     =  0.9;
	for i = 1:min([prod(size(P)) 12])
		text(x,y,[sprintf('%-4.0f',i) P(i).fname],'FontSize',10,'Interpreter','none','Parent',ax);
		y = y - 0.08;
	end
	if prod(size(P)) > 12
		text(x,y,'................ etc','FontSize',10,'Parent',ax); end

	ax=axes('Position',[0.1 0.35 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on');
	plot(Params(:,1:3),'Parent',ax)
	s = ['x translation';'y translation';'z translation'];
	text([2 2 2], Params(2, 1:3), s, 'Fontsize',10,'Parent',ax)
	set(get(ax,'Title'),'String','translation','FontSize',16,'FontWeight','Bold');
	set(get(ax,'Xlabel'),'String','image');
	set(get(ax,'Ylabel'),'String','mm');


	ax=axes('Position',[0.1 0.05 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on');
	plot(Params(:,4:6)*180/pi,'Parent',ax)
	s = ['pitch';'roll ';'yaw  '];
	text([2 2 2], Params(2, 4:6)*180/pi, s, 'Fontsize',10,'Parent',ax)
	set(get(ax,'Title'),'String','rotation','FontSize',16,'FontWeight','Bold');
	set(get(ax,'Xlabel'),'String','image');
	set(get(ax,'Ylabel'),'String','degrees');

	% print realigment parameters
	spm_print
end
return;
%_______________________________________________________________________

%_______________________________________________________________________
function edit_defaults
global MODALITY sptl_WhchPtn sptl_DjstFMRI sptl_CrtWht sptl_MskOptn SWD

tmp = 1;
if sptl_WhchPtn == 1, tmp = 2; end;
sptl_WhchPtn = spm_input(...
	['Coregistration and reslicing?'],...
	2, 'm',...
	['Allow separate coregistration and reslicing|'...
	 'Combine coregistration and reslicing'], [-1 1], tmp);

tmp = 2;
if sptl_CrtWht == 1,
	tmp = 1;
end
sptl_CrtWht   = spm_input(['Images to create?'], 3, 'm',...
	'All Images + Mean Image|Full options', [1 -1], tmp);

tmp = 3;
if sptl_DjstFMRI == 1,
	tmp = 1;
elseif sptl_DjstFMRI == 0
	tmp = 2;
end
sptl_DjstFMRI = spm_input(['fMRI adjustment for interpolation?'],4,'m',...
		'   Always adjust |    Never adjust|Optional adjust',...
		[1 0 -1], tmp);

tmp = 2;
if sptl_MskOptn == 1,
	tmp = 1;
end
sptl_MskOptn  = spm_input(['Option to mask images?'], 5, 'm',...
		'  Always mask|Optional mask', [1 -1], tmp);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function reslice_images_volbyvol(P,Flags,sessions)
% Reslices images volume by volume
% FORMAT reslice_images_volbyvol(P,Flags,sessions)
%
% P     - matrix of filenames {one string per row}
%         All operations are performed relative to the first image.
%         ie. Coregistration is to the first image, and resampling
%         of images is into the space of the first image.
%
% sessions - the last scan in each of the sessions.  For example,
%            the images in the second session would be
%            P((sessions(2-1)+1):sessions(2),:).
%
% Flags - options flags
%
%         k - mask output images
%             To avoid artifactual movement-related variance the realigned
%             set of images can be internally masked, within the set (i.e.
%             if any image has a zero value at a voxel than all images have
%             zero values at that voxel).  Zero values occur when regions
%             'outside' the image are moved 'inside' the image during
%             realignment.
%
%         i - write mean image
%             The average of all the realigned scans is written to
%             mean*.img.
%
%         S - use sinc interpolation for reslicing (9x9x9).
%
%         n - don't reslice the first image
%             The first image is not actually moved, so it may not be
%             necessary to resample it.
%
%         N - don't reslice any of the images - except possibly create a
%             mean image.
%
%             The spatially realigned images are written to the orginal
%             subdirectory with the same filename but prefixed with an 'r'.
%             They are all aligned with the first.
%
% The only reason for reslicing a series plane by plane is to do the adjustment.


P=spm_vol(P);

if any(Flags == 'F'),
	% Check the matrixes
	for i=1:prod(size(P)),
		pp = spm_imatrix(P(1).mat\P(i).mat);
		if any(abs(pp(7:12)-[1 1 1 0 0 0]) > 1e-7),
			error('Can''t do non isotropic voxels!');
		end;
	end;
end;

linfun = inline('fprintf(''%-60s%s'', x,sprintf(''\b'')*ones(1,60))');

if any(Flags == 'k') | any(Flags == 'i'),
	linfun('Computing mask..');
	spm_progress_bar('Init',P(1).dim(3),'Computing available voxels','planes completed');
	x1    = repmat((1:P(1).dim(1))',1,P(1).dim(2));
	x2    = repmat( 1:P(1).dim(2)  ,P(1).dim(1),1);
	if any(Flags == 'i'), Count = zeros(P(1).dim(1:3)); end;
	if any(Flags == 'k'), msk   = cell(P(1).dim(3),1);  end;
	for x3 = 1:P(1).dim(3),
		tmp = zeros(P(1).dim(1:2));
		for i = 1:prod(size(P)),
			tmp = tmp + getmask(inv(P(1).mat\P(i).mat),x1,x2,x3,P(i).dim(1:3));
		end;
		if any(Flags == 'k'), msk{x3} = find(tmp ~= prod(size(P))); end;
		if any(Flags == 'i'), Count(:,:,x3) = tmp; end;
		spm_progress_bar('Set',x3);
	end;
end;

linfun('Reslicing images..');
Hold = 1;
if (any(Flags == 'S')) Hold = -9; end
spm_progress_bar('Init',prod(size(P)),'Reslicing','volumes completed');


start_vol = 1;
if any(Flags == 'n') & ~any(Flags == 'i'), start_vol = 2; end;
if any(Flags == 'i'), Integral = zeros(P(1).dim(1:3)); end;

tiny = 5e-2; % From spm_vol_utils.c

PO = P;
for i = start_vol:prod(size(P)),
	linfun(['Reslicing volume ' num2str(i) '..']);

	PO(i).fname   = prepend(P(i).fname,'r');
	PO(i).dim     = [P(1).dim(1:3) P(i).dim(4)];
	PO(i).mat     = P(1).mat;
	PO(i).descrip = 'spm - realigned';

	if any(Flags == 'F'),
		v = abs(kspace3d(loadvol(P(i)),P(1).mat\P(i).mat));
		for x3 = 1:P(1).dim(3),
			if any(Flags == 'i'),
				Integral(:,:,x3) = Integral(:,:,x3) + ...
					v(:,:,x3).*getmask(inv(P(1).mat\P(i).mat),x1,x2,x3,P(i).dim(1:3));
			end;
			if any(Flags == 'k'),
				tmp          = v(:,:,x3);
				tmp(msk{x3}) = NaN;
				v(:,:,x3)    = tmp;
			end;
		end;
		if ((i>1) | ~any(Flags == 'n')) & ~any(Flags == 'N'), spm_write_vol(PO(i),v); end;
	else,
		if ((i>1) | ~any(Flags == 'n')) & ~any(Flags == 'N'), spm_create_image(PO(i)); end;
		for x3 = 1:P(1).dim(3),
			if (i>1) | ~any(Flags == 'n') | any(Flags == 'i'),
				M = inv(spm_matrix([0 0 -x3 0 0 0 1 1 1])*inv(P(1).mat)*P(i).mat);
				v = spm_slice_vol(P(i),M,P(1).dim(1:2),Hold);
			end;
			if any(Flags == 'i'),
				Integral(:,:,x3) = Integral(:,:,x3) + v;
			end;

			if ((i>1) | ~any(Flags == 'n')) & ~any(Flags == 'N'),
				if any(Flags == 'k'), v(msk{x3}) = NaN; end;
				spm_write_plane(PO(i),v,x3);
			end;
		end;
	end;
	spm_progress_bar('Set',i);
end;

if any(Flags == 'i')
	% Write integral image (16 bit signed)
	%-----------------------------------------------------------
	Integral   = Integral./Count;
	PO         = P(1);
	PO.fname   = prepend(P(1).fname, 'mean');
	PO.pinfo   = [max(max(max(Integral)))/32767 0 0]';
	PO.descrip = 'spm - mean image';
	PO.dim(4)  = 4;
	spm_write_vol(PO,Integral);
end

spm_figure('Clear','Interactive');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function v = kspace3d(v,M)
% 3D rigid body transformation performed as shears in 1D Fourier space.
% FORMAT v1 = kspace3d(v,M)
% Inputs:
% v - the image stored as a 3D array.
% M - the rigid body transformation matrix.
% Output:
% v - the transformed image.
%
% The routine is based on the excellent papers:
% R. W. Cox and A. Jesmanowicz (1999)
% Real-Time 3D Image Registration for Functional MRI
% Submitted to MRM (April 1999) and avaliable from:
% http://varda.biophysics.mcw.edu/~cox/index.html.
% and:
% W. F. Eddy, M. Fitzgerald and D. C. Noll (1996)
% Improved Image Registration by Using Fourier Interpolation
% Magnetic Resonance in Medicine 36(6):923-931
%_______________________________________________________________________

[S0,S1,S2,S3] = shear_decomp(M);

d  = [size(v) 1 1 1];
g = 2.^ceil(log2(d));
if any(g~=d),
	tmp = v;
	v   = zeros(g);
	v(1:d(1),1:d(2),1:d(3)) = tmp;
	clear tmp;
end;

% XY-shear
tmp1 = -sqrt(-1)*2*pi*([0:((g(3)-1)/2) 0 (-g(3)/2+1):-1])/g(3);
for j=1:g(2),
	t        = reshape( exp((j*S3(3,2) + S3(3,1)*(1:g(1)) + S3(3,4)).'*tmp1) ,[g(1) 1 g(3)]);
	v(:,j,:) = real(ifft(fft(v(:,j,:),[],3).*t,[],3));
end;

% XZ-shear
tmp1 = -sqrt(-1)*2*pi*([0:((g(2)-1)/2) 0 (-g(2)/2+1):-1])/g(2);
for k=1:g(3),
	t        = exp( (k*S2(2,3) + S2(2,1)*(1:g(1)) + S2(2,4)).'*tmp1);
	v(:,:,k) = real(ifft(fft(v(:,:,k),[],2).*t,[],2));
end;

% YZ-shear
tmp1 = -sqrt(-1)*2*pi*([0:((g(1)-1)/2) 0 (-g(1)/2+1):-1])/g(1);
for k=1:g(3),
	t        = exp( tmp1.'*(k*S1(1,3) + S1(1,2)*(1:g(2)) + S1(1,4)));
	v(:,:,k) = real(ifft(fft(v(:,:,k),[],1).*t,[],1));
end;

% XY-shear
tmp1 = -sqrt(-1)*2*pi*([0:((g(3)-1)/2) 0 (-g(3)/2+1):-1])/g(3);
for j=1:g(2),
	t        = reshape( exp( (j*S0(3,2) + S0(3,1)*(1:g(1)) + S0(3,4)).'*tmp1) ,[g(1) 1 g(3)]);
	v(:,j,:) = real(ifft(fft(v(:,j,:),[],3).*t,[],3));
end;

if any(g~=d), v = v(1:d(1),1:d(2),1:d(3)); end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [S0,S1,S2,S3] = shear_decomp(A)
% Decompose rotation and translation matrix A into shears S0, S1, S2 and
% S3, such that A = S0*S1*S2*S3.  The original procedure is documented
% in:
% R. W. Cox and A. Jesmanowicz (1999)
% Real-Time 3D Image Registration for Functional MRI

A0 = A(1:3,1:3);
if any(abs(svd(A0)-1)>eps*1000), error('Can''t decompose matrix'); end;

t  = A0(2,3); if t==0, t=eps; end;
a0 = pinv(A0([1 2],[2 3])')*[(A0(3,2)-(A0(2,2)-1)/t) (A0(3,3)-1)]';
S0 = [1 0 0; 0 1 0; a0(1) a0(2) 1];
A1 = S0\A0;  a1 = pinv(A1([2 3],[2 3])')*A1(1,[2 3])';  S1 = [1 a1(1) a1(2); 0 1 0; 0 0 1];
A2 = S1\A1;  a2 = pinv(A2([1 3],[1 3])')*A2(2,[1 3])';  S2 = [1 0 0; a2(1) 1 a2(2); 0 0 1];
A3 = S2\A2;  a3 = pinv(A3([1 2],[1 2])')*A3(3,[1 2])';  S3 = [1 0 0; 0 1 0; a3(1) a3(2) 1];

s3 = A(3,4)-a0(1)*A(1,4)-a0(2)*A(2,4);
s1 = A(1,4)-a1(1)*A(2,4);
s2 = A(2,4);
S0 = [[S0 [0  0 s3]'];[0 0 0 1]];
S1 = [[S1 [s1 0  0]'];[0 0 0 1]];
S2 = [[S2 [0 s2  0]'];[0 0 0 1]];
S3 = [[S3 [0  0  0]'];[0 0 0 1]];
return;
%_______________________________________________________________________

%_______________________________________________________________________
function Mask = getmask(M,x1,x2,x3,dim)
tiny = 5e-2; % From spm_vol_utils.c
y1=M(1,1)*x1+M(1,2)*x2+(M(1,3)*x3+M(1,4));
y2=M(2,1)*x1+M(2,2)*x2+(M(2,3)*x3+M(2,4));
y3=M(3,1)*x1+M(3,2)*x2+(M(3,3)*x3+M(3,4));
Mask =        (y1 >= (1-tiny) & y1 <= (dim(1)+tiny));
Mask = Mask & (y2 >= (1-tiny) & y2 <= (dim(2)+tiny));
Mask = Mask & (y3 >= (1-tiny) & y3 <= (dim(3)+tiny));
return;
%_______________________________________________________________________

%_______________________________________________________________________
function v = loadvol(V)
v = zeros(V.dim(1:3));
for i=1:V.dim(3),
	v(:,:,i) = spm_slice_vol(V,spm_matrix([0 0 i]),V.dim(1:2),0);
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function A = interp_errors(img,Hold)
fimg = fft2(img);
t1 = 0.5;
t2 = 0.25;
A = zeros([size(img) 4]);
A(:,:,1) = interp_errors_sub1(fimg,[0 t1],Hold)-img;
A(:,:,2) = interp_errors_sub1(fimg,[0 t2],Hold)-img;
A(:,:,3) = interp_errors_sub1(fimg,[t1 0],Hold)-img;
A(:,:,4) = interp_errors_sub1(fimg,[t2 0],Hold)-img;
M = inv([1-cos(t1*2*pi) 1-cos(t2*2*pi); sin(t1*2*pi) sin(t2*2*pi)]);
M = [M zeros(2); zeros(2) M];
A = reshape(reshape(A,[prod(size(img)) 4])*M,[size(img) 4]);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function img2 = interp_errors_sub1(fimg,t,Hold)
d    = size(fimg);
tr   = t(2);
if tr ~= 0,
	c    = fftshift(exp(sqrt(-1)*2*pi*tr*((-d(1)/2):((d(1)-1)/2))/d(1)))';
	fimg = fimg.*repmat(c,1,d(2));
end;
tr   = -t(1);
if tr ~= 0,
	c    = fftshift(exp(sqrt(-1)*2*pi*tr*((-d(2)/2):((d(2)-1)/2))/d(2)));
	fimg = fimg.*repmat(c,d(1),1);
end;
img1 = real(ifft2(fimg));
img2 = spm_slice_vol(img1,spm_matrix([t(2) t(1) 1]),size(img1),Hold);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function save_parameters(V)
fname = [spm_str_manip(prepend(V(1).fname,'realignment_params_'),'s') '.txt'];
n = length(V);
Q = zeros(n,6);
for j=1:n,
	qq     = spm_imatrix(V(j).mat/V(1).mat);
	Q(j,:) = qq(1:6);
end;
save(fname,'Q','-ascii');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function PO = prepend(PI,pre)
[pth,nm,xt,vr] = fileparts(deblank(PI));
PO             = fullfile(pth,[pre nm xt vr]);
return;
%_______________________________________________________________________
