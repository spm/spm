function spm_coreg_ui(opt)
% Between modality coregistration using Mutual Information
%
% ____________________________________________________________________________
% 
%  The TARGET image is the image to which the OBJECT image is realigned.
%  If there are any OTHER images, then the same transformations are applied to
%  these images as are applied to the OBJECT image.
% 
%  eg 1) to realign a structural MR image to a sequence of PET images:
%   TARGET: meanPET1.img
%   OBJECT: MRI.img
%    OTHER: -
% 
%  eg 2) to realign a sequence of PET images to a structural MR image:
%   TARGET: MRI.img
%   OBJECT: meanPET1.img
%    OTHER: PET1.img PET2.img PET3.img etc...
% ____________________________________________________________________________
%
% The registration method used here is based on the work described in:
% A Collignon, F Maes, D Delaere, D Vandermeulen, P Suetens & G Marchal
% (1995) "Automated Multi-modality Image Registration Based On
% Information Theory". In the proceedings of Information Processing in
% Medical Imaging (1995).  Y. Bizais et al. (eds.).  Kluwer Academic
% Publishers.
%
% The original interpolation method described in this paper has been
% changed in order to give a smoother cost function.  The images are
% also smoothed slightly, as is the histogram.  This is all in order to
% make the cost function as smooth as possible, to give faster
% convergence and less chance of local minima.
%
% References
% ==========
% Mutual Information
% ------------------
% Collignon, Maes, Delaere, Vandermeulen, Suetens & Marchal (1995).
% "Automated multi-modality image registration based on information theory".
% In Bizais, Barillot & Di Paola, editors, Proc. Information Processing
% in Medical Imaging, pages 263--274, Dordrecht, The Netherlands, 1995.
% Kluwer Academic Publishers.
%
% Wells III, Viola, Atsumi, Nakajima & Kikinis (1996).
% "Multi-modal volume registration by maximisation of mutual information".
% Medical Image Analysis, 1(1):35-51, 1996.
%
% Entropy Correlation Coefficient
% -------------------------------
% F Maes, A Collignon, D Vandermeulen, G Marchal & P Suetens (1997).
% "Multimodality image registration by maximisation of mutual
% information". IEEE Transactions on Medical Imaging 16(2):187-198
%
% Normalised Mutual Information
% -----------------------------
% Studholme,  Hill & Hawkes (1998).
% "A normalized entropy measure of 3-D medical image alignment".
% in Proc. Medical Imaging 1998, vol. 3338, San Diego, CA, pp. 132-143.
%
% Optimisation
% ------------
% Press, Teukolsky, Vetterling & Flannery (1992).
% "Numerical Recipes in C (Second Edition)".
% Published by Cambridge.
%
%_______________________________________________________________________
%
% At the end of coregistration, the voxel-to-voxel affine transformation
% matrix is displayed, along with the histograms for the images in the
% original orientations, and the final orientations.  The registered
% images are displayed at the bottom.
% 
% Registration parameters are stored in the ".mat" files of the "source"
% and the "other" images.
%
%
%_______________________________________________________________________
%
%                        The Prompts Explained
%_______________________________________________________________________
%
% 'number of subjects'
% The number of pairs of images that are matched.
%
% 'Which option?'
% 	'Coregister only'
% 	 Only estimate the rigid body transform and write/modify the
% 	 .mat files.
% 	'Reslice Only'
% 	 Resample the images according to their .mat files.
% 	'Coregister & Reslice'
% 	 Estimate the transforms and resample the images.
%
% When doing the Coregister step:
%
% 'Target image, subj #'
% This is the image that is assumed to remain stationary (sometimes
% known as the reference or template image), while the source image
% is moved to match it.
%
% 'Source image, subj #'
% This is the image that is jiggled about to best match the target.
%
% 'Other images, subj #'
% These are any images that need to remain in alignment with the
% source image.
%
% At the end, the voxel-to-voxel affine transformation matrix is
% displayed, along with the histograms for the images in the original
% orientations, and the final orientations.  The registered images are
% displayed at the bottom.
%
% Registration parameters are stored in the ".mat" files of the "source"
% and the "other" images.                                                                                        
%
%
% When doing the Reslice step:
%
% 'Space defining image, subject #'
% This is analagous to the target image.  Images are resliced to match
% this image (providing they have been coregistered first).
%
% 'Images to reslice, subj #'
% These images are resliced to the same dimensions, voxel sizes,
% orientation etc as the space defining image.
%
%_______________________________________________________________________
%
%                           Defaults Options
%_______________________________________________________________________
%[   things in square brackets indicate corresponding defaults field   ]
% 
% Defaults for Coregistration:
%
% 'Cost Function?'
% For inter-modal registration, use:
% 	'mi'  - Mutual Information
% 	'nmi' - Normalised Mutual Information
% 	'ecc' - Entropy Correlation Coefficient
% For within modality, you could use:
% 	'ncc' - Normalised Cross Correlation
% [defaults.coreg.estimate.cost_fun]
%
%
% Defaults for Reslicing:
%
% 'Interpolation Method?'
% The method by which the images are sampled when being written in a
% different space.
%       'Nearest Neighbour'
%               - Fastest, but not normally recommended.
%       'Bilinear Interpolation'
%               - OK for PET, or realigned fMRI.
%
%       'B-spline Interpolation'
%               - Better quality (but slower) interpolation, especially
%                 with higher degree splines.  Don't use B-splines when
%                 there is any region of NaN or Inf in the images.
% [defaults.coreg.write.interp]
%
% 'Way to wrap images?'
% These are typically:
%       'No wrapping' - for PET or images that have already
%                       been spatially transformed.
%       'Wrap in  Y'  - for (un-resliced) MRI where phase encoding
%                       is in the Y direction (voxel space).
% [defaults.coreg.write.wrap]
%
% 'Mask images?'
% Because of subject motion, different images are likely to have different
% patterns of zeros from where it was not possible to sample data.
% With masking enabled, the program searches through the whole time series
% looking for voxels which need to be sampled from outside the original
% images. Where this occurs, that voxel is set to zero for the whole set
% of images (unless the image format can represent NaN, in which case
% NaNs are used where possible).
% [defaults.coreg.write.mask]
%
% Other Defaults are useful for making spm_coregister.m more flexible:
% [defaults.coreg.estimate.sep] - average distance between sampled points
% (in mm).  Can be a vector to allow a coarse registration followed by a
% more accurate one.
% [defaults.coreg.estimate.tol] - accuracy for each parameter.  Iterations
% stop when differences between sucessive estimates are less than the
% required tolerence.
% [defaults.coreg.estimate.fwhm] - Gaussian smoothing to apply to the
% 256x256 joint histogram. Other information theoretic coregistration
% methods use fewer bins, but Gaussian smoothing seems to be more elegant.
%
%__________________________________________________________________________
%
% The `.mat' files.
%
% This simply contains a 4x4 affine transformation matrix in a variable
% `mat'. What these matrixes contain is a mapping from
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
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id$


global defaults
 
if nargin==0 | strcmp(lower(opt),'ui'),
        run_ui(defaults.coreg);
elseif nargin>0 & strcmp(lower(opt),'defaults'),
        defaults.coreg = get_defs(defaults.coreg);
end;
return;

function run_ui(flags)
SPMid = spm('FnBanner',mfilename,'$Rev$');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Coregister');
spm_help('!ContextHelp',mfilename);

% get number of subjects
nsubjects = spm_input('number of subjects',1, 'e', 1);
if nsubjects < 1,
	spm_figure('Clear','Interactive');
	return;
end;

p = spm_input('Which option?',2,'m',...
	'Coregister only|Reslice Only|Coregister & Reslice', [1 2 3],3);

if p == 1 | p == 3,
	for i = 1:nsubjects,
		mireg(i)    = struct('VG',[],'VF',[],'PO','');
		
		% select target(s)
		PG          = spm_select(1,'image', ['Target image, subj ' num2str(i)]);
		mireg(i).VG = spm_vol(PG);
		
		% select source(s)
		PF          = spm_select(1,'image', ['Source image, subj ' num2str(i)]);
		mireg(i).VF = spm_vol(PF);

		PO = spm_select(Inf,'image', ['Other images, subj ' num2str(i)]);
		if isempty(PO),
			mireg(i).PO = PF;
		else,
			mireg(i).PO = strvcat(PF,PO);
		end;
	end;
end;

if p==2,
	for i = 1:nsubjects,
		mireg(i) = struct('VG',[],'VF',[],'PO',[]);

		% select target space
		PG          = spm_select(1,'image', ['Space defining image, subj ' num2str(i)]);
		mireg(i).VG = spm_vol(PG);

		PO          = spm_select(Inf,'image', ['Images to reslice, subj ' num2str(i)]);
		mireg(i).PO = PO;
	end;
end;

% For each subject, call the program to perform the registration.
%-----------------------------------------------------------------------
spm('Pointer','Watch')
for i=1:nsubjects,

	if p == 1 | p == 3,
		spm('FigName',['Coregister subj ' num2str(i)],Finter,CmdLine);
		x  = spm_coreg(mireg(i).VG, mireg(i).VF,flags.estimate);
		M  = inv(spm_matrix(x));
		MM = zeros(4,4,size(mireg(i).PO,1));
		for j=1:size(mireg(i).PO,1),
			MM(:,:,j) = spm_get_space(deblank(mireg(i).PO(j,:)));
		end;
		for j=1:size(mireg(i).PO,1),
			spm_get_space(deblank(mireg(i).PO(j,:)), M*MM(:,:,j));
		end;
	end;
	if p == 2 | p == 3,
		 spm('FigName',['Reslice subj ' num2str(i)],Finter,CmdLine);
		fprintf('Reslicing Subject %d\n', i);
		P         = strvcat(mireg(i).VG.fname,mireg(i).PO);
		flg       = flags.write;
		flg.mean  = 0;
		flg.which = 1;
		flg.mean  = 0;
		spm_reslice(P,flg);
	end;
end;
spm('FigName','Coregister: done',Finter,CmdLine);
spm('Pointer');
return;

function defs = get_defs(defs)
fun = lower(defs.estimate.cost_fun);

funs={'mi','ecc','nmi','ncc'};
sel = 1;
for i=1:length(funs),
	if strcmp(funs{i},fun), sel = i; break; end;
end;
defs.estimate.cost_fun = spm_input('Cost Function?','+1','m',[...
	'Mutual Information|Entropy Correlation Coefficient|'...
	'Normalised Mutual Information|Normalised Cross Correlation'],...
	funs,sel);
defs.estimate.cost_fun = defs.estimate.cost_fun{1};


tmp2 = [0 1 2 3 4 5 6 7 Inf];
tmp = find(defs.write.interp == tmp2);
if ~finite(defs.write.interp), tmp = 9; end;
if isempty(tmp), tmp = 2; end;
defs.write.interp = spm_input('Reslice interpolation method?','+1','m',...
	['Nearest Neighbour|Trilinear|2nd Degree B-Spline|'...
	'3rd Degree B-Spline|4th Degree B-Spline|5th Degree B-Spline|'...
	'6th Degree B-Spline|7th Degree B-Spline|Fourier Interpolation'],...
	tmp2,tmp);

wraps = [0 0 0 ; 1 0 0; 0 1 0; 1 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 1];
t     = find(all(repmat(defs.write.wrap(:)',8,1) == wraps, 2));
if isempty(t), t = 1; end;
p     = spm_input('Way to wrap images?','+1','m',...
	['No wrap|Wrap X|Wrap Y|Wrap X & Y|Wrap Z|Wrap X & Z|Wrap Y & Z|Wrap X, Y & Z'],...
	[1 2 3 4 5 6 7 8], t);
defs.write.wrap    = wraps(p,:);
defs.estimate.wrap = defs.write.wrap;
 
tmp = 2;
if defs.write.mask == 1, tmp = 1; end;
defs.write.mask  = spm_input(['Mask images?'], '+1', 'm',...
	'  Mask images|Dont mask images', [1 0], tmp);
return;


