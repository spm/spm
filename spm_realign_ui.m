function spm_realign_ui(opt)
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
%_______________________________________________________________________
%
% In addition to realigning (transforming the images according to a rigid
% body model) the time series there is the option of modelling the residual 
% movement related variance in an EPI (fMRI) time series that can be 
% explained by a model for suceptibility-by-movement interactions. This
% is given by the "Realign & Unwarp" option.
%
% Susceptibility artefacts in EPI time series is a consequence of the
% "field disturbances" caused by the presence of an object in the field.
% Susceptibility is a property of a material, and can be thought of as
% the "resistance" put up by that material against being magnetised.
% At the interface between materials with different susceptibility
% rather severe field disturbancies will ensue. I.e. the field is
% no longer nice and homogenous. Since varying field strenght (gradients)
% is used to encode position in MRI, these disturbances cause spatial
% misplacement of RF signal, i.e. geometric distortions. Because of
% the low bandwidth in the phase encode direction these distortions
% will be appreciable mainly in that direction (typically the Ant-Post
% direction). The distortions are noticable mainly near air-tissue
% interfaces inside the scull, i.e. near the sinuses and the auditory
% canals. 
%
% In these areas in particular the observed image is a severly warped 
% version of reality, much like a funny mirror at a fair ground. When 
% one moves in front of such a mirror ones image will distort in 
% different ways and ones head may change from very elongated to 
% seriously flattened. If we were to take digital snapshots of the 
% reflection at these different positions it is rather obvious that 
% realignment alone will not suffice to bring them into a common space. 
%
% The situation is similar with EPI images, and an image collected 
% for a given subject position will not be identical to that collected 
% at another. We call this effect suscebtibility-by-movement interaction. 
% The "Unwarp" toolbox is predicated on the assumption that the 
% suscebtibility-by-movement interaction is responsible for a sizeable 
% part of residual movement related variance.
%
% Assume that we know how the deformations change when the subject changes 
% position (i.e. we know the derivatives of the deformations with respect 
% to subject position). That means that for a given time series and a 
% given set of subject movements we should be able to predict the "shape 
% changes" in the object and the ensuing variance in the time series. 
% It also means that, in principle, we should be able to formulate the 
% inverse problem, i.e. given the observed variance (after realignment) 
% and known (estimated) movements we should be able to estimate how 
% deformations change with subject movement. 
%
% We have made an attempt at formulating such an inverse model, and at 
% solving for the "derivative fields". A deformation field can be thought 
% of as little vectors at each position in space showing how that particular 
% location has been deflected. A "derivative field" is then the rate of 
% change of those vectors with respect to subject movement. Given these 
% "derivative fields" we should be able to remove the variance caused by 
% the suscebtibility-by-movement interaction. Since the underlying model 
% is so restricted we would also expect experimentally induced variance 
% to be preserved. 
%
% Refs
%
% Background reading:
%
% Friston KJ, Williams SR, Howard R, Frackowiak RSJ and Turner R (1995)
% Movement-related effect in fMRI time-series.  Mag. Res. Med. 35:346-355
%
% Jezzard P and Balaban RS (1995) Correction for geometric distortions
% in echoplanar images from B0 field variations. Magn Reson Med 34:65-73
%
% Wu DH, Lewin JS and Duerk JL (1997) Inadequacy of motion correction
% algorithms in functional MRI: Role of susceptibility-induced artefacts.
% J Magn Reson Imag 7:365-370
%
% About this particular method:
%
% Andersson JLR, Hutton C, Ashburner J, Turner R, Friston K (2001)
% Modelling geometric deformations in EPI time series. NeuroImage
% 13:903-919. doi:10.1006/nimg.2001.0746
%
%_______________________________________________________________________
%
%                        The Prompts Explained
%_______________________________________________________________________
%
% 'Number of subjects'
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
% 'Images, subject # ..'
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
%_______________________________________________________________________
%
%                           Defaults Options
%_______________________________________________________________________
%
% 'Registration Quality?'
% Quality versus speed trade-off.  Highest quality (1) gives most
% precise results, whereas lower qualities gives faster realignment.
% The idea is that some voxels contribute little to the estimation of
% the realignment parameters. This parameter is involved in selecting
% the number of voxels that are used.
% [defaults.realign.estimate.quality]
%
% 'Allow weighting of reference image?'
% Give the option of providing a weighting image to weight each voxel
% of the reference image differently when estimating the realignment
% parameters.  The weights are proportional to the inverses of the
% standard deviations.
% For example, when there is a lot of extra-brain motion - e.g., during
% speech, or when there are serious artifacts in a particular region of
% the images.
% [defaults.realign.estimate.weight]
%
%
% 'Reslice interpolation Method?'
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
%       'Fourier space Interpolation'
%               - Rigid body rotations are executed as a series of shears,
%                 which are performed in Fourier space (Eddy et. al. 1996).
%                 Unfortunately, this method can only be applied to images
%                 with cubic voxels (since zooms can not be done by
%                 convolution in Fourier space). 
% [defaults.realign.write.interp]
%
% These are typically:
%       'No wrapping' - for PET or images that have already
%                       been spatially transformed.
%       'Wrap in  Y'  - for (un-resliced) MRI where phase encoding
%                       is in the Y direction (voxel space).
% [defaults.realign.write.wrap]
%
% 'Mask images?'
% Because of subject motion, different images are likely to have different
% patterns of zeros from where it was not possible to sample data.
% With masking enabled, the program searches through the whole time series
% looking for voxels which need to be sampled from outside the original
% images. Where this occurs, that voxel is set to zero for the whole set
% of images (unless the image format can represent NaN, in which case
% NaNs are used where possible). This is in order to avoid artifactual
% movement-related variance the realigned images.
% [defaults.realign.write.mask]
%
% Other Defaults are useful for making spm_realign_ui.m more flexible:
% [defaults.realign.estimate.interp] - interpolation method for parameter
% estimation.
% [defaults.realign.estimate.wrap] - wrapping used for parameter estimation.
% This should probably also influence the way that smoothing is done, but
% the smoothing code has not yet been modified to include wrapping.
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
%
% The *uw.mat file.
%
% The *uw.mat file contains information about input to the Unwarp module 
% (i.e. how it was run) and about the result. It contains all the information
% needed to unwarp the appurtenant time-series, and to assess exactly how the
% estimation was done.
%__________________________________________________________________________
% %W% John Ashburner - with input from Oliver Josephs %E%
% and Jesper Andersson

global defaults

if nargin==0 | strcmp(lower(opt),'ui'),
   if get(gcbo,'Value') < 2
	run_ui(defaults.realign, defaults.modality);
   else
        run_ui(defaults.realign, defaults.modality, defaults.unwarp);
   end
elseif nargin>0 & strcmp(lower(opt),'defaults'),
	defaults.realign = get_defs(defaults.realign);
elseif nargin>0 & strcmp(lower(opt),'unwarpdefaults'),
        defaults.unwarp = get_unwarp_defs(defaults.realign);
end;
return;

function run_ui(defs, modality, unwarp)
% User interface.
%_______________________________________________________________________
SPMid                   = spm('FnBanner',mfilename,'2.10');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Realign');
spm_help('!ContextHelp',mfilename);

n     = spm_input('Num subjects', '+1', 'e', 1);
if n<1, spm_figure('Clear','Interactive'); return; end

%
% Possibility to include measured field-map temporarily removed.
%
pcpm = 0;
%if strcmp(lower(modality),'fmri') & nargin > 2
%	pcpm = spm_input('Include pre-calculated phase maps?','+1','m',...
%		[' Yes| No'],[1 0],0);
%end
   

P = cell(n,1);
for i = 1:n,
	if strcmp(lower(modality),'fmri'),
		ns = spm_input(['Num sessions, subj ' num2str(i)], '+1', 'e', 1);
		pp = cell(1,ns);
		for s=1:ns,
			p = '';
                        pm = '';
			while size(p,1)<1,
				p = spm_get(Inf,'IMAGE',...
					['Images, subj ' num2str(i) ', sess' num2str(s)]);
				if pcpm == 1
					pm = spm_get(1,'IMAGE',...
					['Phasemap, subj ' num2str(i) ', sess' num2str(s)]);
				end
			end;
			pp{s} = p;
                        ppm{s} = pm;
		end;
		P{i} = pp;
                pP{i} = ppm;
	else,
		p  = cell(1,1);
		p{1} = '';
		while size(p{1},1)<1,
		      p{1} = spm_get(Inf,'IMAGE',...
			  ['Images, subj ' num2str(i)]);
		end;
		P{i} = p;
	end;
end;


if strcmp(lower(modality),'pet'),
	FlagsC = struct('quality',defs.estimate.quality,'fwhm',7,'rtm',1);
else,
	FlagsC = struct('quality',defs.estimate.quality,'fwhm',5,'rtm',0);
end;

if nargin < 3
	WhchPtn = spm_input('Which option?', '+1', 'm',...
		'Coregister only|Reslice Only|Coregister & Reslice',...
		[1 2 3],3);
else
	WhchPtn = 3;
end

PW = '';
if (WhchPtn == 1 | WhchPtn == 3) & defs.estimate.weight,
	if spm_input(...
		['Weight the reference image(s)?'],...
		'+1', 'm',...
		['Dont weight registration|'...
		 'Weight registration'], [0 1], 1);

		PW = spm_get(n,'IMAGE', 'Weight images for each subj');
	end;
end;

if nargin > 2
	foe = spm_input('Model field changes w.r.t.','+1','m',...
		['Pitch & Roll|All|Customise'],...
		[1:2 0],1);
	if foe == 1
		foe = [4 5];
	elseif foe == 2
		foe = 1:6;
	else
		foe = spm_input('First order effects','+0','e',...
			'4 5',Inf)';
	end
	if unwarp.estimate.soe == 2
		tmp = spm_input('Second order effects','+1','m',...
		['None|All first order effects|Customise'],...
		[1:3],1);
		if tmp == 1
			soe = [];
		else
			cnt = 1;
			for i=1:size(foe,2)
				for j=i:size(foe,2)
					soe(cnt,1) = foe(i);
					soe(cnt,2) = foe(j);
					cnt = cnt+1;
				end
			end
			if tmp == 3
				string = '[';
				for i=1:cnt-1
					string = [string sprintf('%d %d; ',soe(i,1),soe(i,2))];
				end
				string = [string(1:end-2) ']'];
				soe = spm_input('Second order effects','+0','e',...
                 			string,[Inf 2]);
			end
		end
	else
		soe = [];
	end
		
end


% Reslicing options
%-----------------------------------------------------------------------
if WhchPtn == 2 | WhchPtn == 3,
	FlagsR = struct('interp',defs.write.interp,...
		'wrap',defs.write.wrap,...
		'mask',defs.write.mask,...
		'which',2,'mean',1);

	if strcmp(lower(modality),'pet'), FlagsR.wrap = [0 0 0]; end;

	if nargin < 3
		p = spm_input('Create what?','+1','m',...
			[' All Images (1..n)| Images 2..n|'...
			 ' All Images + Mean Image| Mean Image Only'],...
			[1 2 3 4],3);
		if p==1, FlagsR.which = 2; FlagsR.mean = 0; end
		if p==2, FlagsR.which = 1; FlagsR.mean = 0; end
		if p==3, FlagsR.which = 2; FlagsR.mean = 1; end
		if p==4, FlagsR.which = 0; FlagsR.mean = 1; end
	else
		p = spm_input('Create what?','+1','m',...
			[' All Images (1..n)| All Images + Mean Image'],[1 3],3);
		if p==1, FlagsR.which = 2; FlagsR.mean = 0; end
		if p==3, FlagsR.which = 2; FlagsR.mean = 1; end
	end

end;

%
% If doing realignment only
%
if nargin < 3
	spm('Pointer','Watch');
	for i = 1:n
		if WhchPtn==1 | WhchPtn==3,
			spm('FigName',['Realigning subj ' num2str(i)],Finter,CmdLine);
			flagsC = FlagsC;
			if ~isempty(PW), flagsC.PW = deblank(PW(i,:)); end;
			spm_realign(P{i},flagsC);
		end
		if WhchPtn==2 | WhchPtn==3,
			spm('FigName',['Reslicing subj ' num2str(i)],Finter,CmdLine);
			spm_reslice(P{i},FlagsR);
		end;
	end;
	spm('FigName','Realign: done',Finter,CmdLine);
	spm('Pointer');
%
% If doing realignment and unwarping
%
else
	uwe_flags = struct('order',    unwarp.estimate.basfcn,...
			'sfield',      [],...
			'regorder',    unwarp.estimate.regorder,...
			'lambda',      unwarp.estimate.regwgt,...
			'jm',          unwarp.estimate.jm,...
			'fot',         foe,...
			'sot',         soe,...
			'fwhm',        unwarp.estimate.fwhm,...
			'rem',         unwarp.estimate.rem,...
			'noi',         unwarp.estimate.noi,...
			'exp_round',   unwarp.estimate.expround);


	uwr_flags = struct('interp',   defs.write.interp,...
			'wrap',        defs.write.wrap,...
			'mask',        defs.write.mask,...
			'which',       FlagsR.which,...
			'mean',        FlagsR.mean);
	if unwarp.estimate.jm == 1
		uwr_flags.udc = 2;
	else
		uwr_flags.udc = 1;
	end

	for i = 1:n
		spm('Pointer','Watch');
		spm('FigName',['Realigning subj ' num2str(i)],Finter,CmdLine);
		flagsC = FlagsC;
		if ~isempty(PW), flagsC.PW = deblank(PW(i,:)); end;
		spm_realign(P{i},flagsC);
		spm('Pointer');
                clear ads
                tmpP = spm_vol(P{i}{1}(1,:));
                uwe_flags.M = tmpP.mat;
		for j=1:length(P{i})
			if pcpm > 0
				uwe_flags.sfield = pP{i}{j};
			end
			ds = spm_uw_estimate(P{i}{j},uwe_flags);
                        ads(j) = ds;
			[path,name,ext,ver] = fileparts(P{i}{j}(1,:));
			pefile = fullfile(path,[name '_uw.mat']);
			save(pefile,'ds');
		end		
		spm_uw_apply(ads,uwr_flags);
	end
end


return;
%_______________________________________________________________________

%_______________________________________________________________________
function defs = get_defs(defs)

tmp2 = [1.00 0.90 0.75 0.50 0.25 0.10 0.05 0.01 0.005 0.001];
tmp = find(defs.estimate.quality == tmp2);
if isempty(tmp) tmp = 1; end
defs.estimate.quality = spm_input('Registration Quality?','+1','m',...
	['Quality 1.00  (slowest/most accurate) |Quality 0.90|' ...
	 'Quality 0.75|Quality 0.50|Quality 0.25|Quality 0.10|' ...
	 'Quality 0.05|Quality 0.01|' ...
	 'Quality 0.005|Quality 0.001 (fastest/poorest)'],tmp2, tmp);

tmp = 1;
if defs.estimate.weight == 1, tmp = 2; end;
defs.estimate.weight = spm_input(...
	['Allow weighting of reference image?'],...
	'+1', 'm',...
	['Dont allow weighting|'...
	 'Allow weighting'], [0 1], tmp);

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

tmp = 1;
if ~defs.write.mask, tmp = 2; end;
defs.write.mask  = spm_input(['Mask images?'], '+1', 'm',...
		'  Mask images|Dont mask images', [1 0], tmp);

return;
%_______________________________________________________________________

%_______________________________________________________________________
function defs = get_unwarp_defs(defs)

orders = [6 6; 8 8; 10 10; 12 12];
lambdas = [1e4 1e5 1e6];

defs.estimate.fwhm = spm_input('Filter width (mm)','+1','e',...
                               '4',1);

defs.estimate.basfcn = spm_input('No. of basis functions','+1','m',...
                                 ['6x6x*|8x8x*|10x10x*|12x12x*'],...
                                 [1:4],3);
defs.estimate.basfcn = orders(defs.estimate.basfcn,:);

defs.estimate.regwgt = spm_input('Amount of regularisation','+1','m',...
                                 ['A little|Medium|A lot'],...
                                 [1:3],2);
defs.estimate.regwgt = lambdas(defs.estimate.regwgt);

defs.estimate.rem = spm_input('Re-estimation of movement parameters?','+1','m',...
                              ['Yes|No'],[1 0],1);

defs.estimate.jm = spm_input('Include Jacobian intensity modulation?','+1','m',...
                             ['Yes|No'],[1 0],2);

defs.estimate.soe = spm_input('Second order effects','+1','m',...
                              ['None|Customise'],...
                              [1:2],1);

return

