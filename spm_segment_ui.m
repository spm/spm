function spm_segment_ui(P)
% Segment an MR image into Grey, White & CSF.
%_______________________________________________________________________
%
%                      The algorithm is three step:
%
% 1) Determine the affine transform which best matches the image with a
%    template image. If the name of more than one image is passed, then
%    the first image is used in this step. This step is not performed if
%    no template images are specified.
%
% 2) Perform Cluster Analysis with a modified Mixture Model and a-priori
%    information about the likelihoods of each voxel being one of a
%    number of different tissue types. If more than one image is passed,
%    then they they are all assumed to be in register, and the voxel
%    values are fitted to multi-normal distributions.
%
% 3) Do a "cleanup" of the partitions, analagous to what <Xtract Brain>
%    followed by <ImCalc> with "i1.*i4./(i1+i2+i3+eps)" did with SPM99.
%
% 4) Write the segmented image. The names of these images have
%    "_seg1", "_seg2" & "_seg3" appended to the name of the
%    first image passed.
%
%_______________________________________________________________________
% Refs:
%
% Ashburner J & Friston KJ (1997) Multimodal Image Coregistration and
% Partitioning - a Unified Framework. NeuroImage 6:209-217
%
% Ashburner J & Friston KJ (2000) Voxel-Based Morphometry - The Methods.
% NeuroImage 11(6):805-821
%
% Ashburner J (2002) "Another MRI Bias Correction Approach" [abstract].
% Presented at the 8th International Conference on Functional Mapping of
% the Human Brain, June 2-6, 2002, Sendai, Japan. Available on CD-Rom
% in NeuroImage, Vol. 16, No. 2.
%
%_______________________________________________________________________
%
%                        The Prompts Explained
%_______________________________________________________________________
%
% 'Scan(s) for subject #'
% If more than one volume is specified (eg T1 & T2), then they must be
% in register (same position, size, voxel dims etc..).  If no images
% are selected, then the interface continues on to the next step.
%
% 'Are they spatially normalised?'
% If not, then an affine registration needs to be done.
% 'Modality?'
% The first image of each subject is matched to a template, so this
% bit prompts for which one to use.
% 'Template(s) for affine matching'
% If none of the default templates are any use, then select one or more
% template images, which will be used for affine matching. Note that the
% affine transform is only determined from the first image specified
% for segmentation. 
%
% 'Attempt to correct intensity inhomogeneities?'
% This uses a Bayesian framework (again) to model intensity
% inhomogeneities in the image(s).  The variance associated with each
% tissue class is assumed to be multiplicative (with the
% inhomogeneities).  The low frequency intensity variability is
% modelled by a linear combination of three dimensional DCT basis
% functions (again), using a fast algorithm (again) to generate the
% curvature matrix.  The regularization is based upon minimizing the
% integral of square of the fourth derivatives of the modulation field
% (the integral of the squares of the first and second derivs give the
% membrane and bending energies respectively).
%
%_______________________________________________________________________
%
%                           Defaults Options
%_______________________________________________________________________
%
% 'Bias regularisation?'
% The importance of smoothness for the estimated bias field. Without
% any regularisation, the algorithm will attempt to correct for
% different grey levels arising from different tissue types, rather than
% just correcting bias artifact.
% Bias correction uses a Bayesian framework (again) to model intensity
% inhomogeneities in the image(s).  The variance associated with each
% tissue class is assumed to be multiplicative (with the
% inhomogeneities).  The low frequency intensity variability is
% modelled by a linear combination of three dimensional DCT basis
% functions (again), using a fast algorithm (again) to generate the
% curvature matrix.  The regularization is based upon minimizing the
% integral of square of the fourth derivatives of the modulation field
% (the integral of the squares of the first and second derivs give the
% membrane and bending energies respectively).
% [defaults.segment.estimate.reg]
%
% 'Bias cutoff?'
% Cutoff of DCT bases.  Only DCT bases of periods longer than the
% cutoff are used to describe the warps. The number used will
% depend on the cutoff and the field of view of the image.
% [defaults.segment.estimate.cutoff]
%
% 'Clean up the partitions?'
% This uses a crude routine for extracting the brain from segmented
% images.  It begins by taking the white matter, and eroding it a
% couple of times to get rid of any odd voxels.  The algorithm
% continues on to do conditional dilations for several iterations,
% where the condition is based upon gray or white matter being present.
% This identified region is then used to clean up the grey and white
% matter partitions, and has a slight influences on the CSF partition.
% [defaults.segment.write.cleanup]
%
% 'Write bias correted image?'
% Pretty self explanitory.
% [defaults.segment.write.wrt_cor]
%
% Other Defaults are useful for making spm_segment.m more flexible: 
% [defaults.segment.estimate.priors] - prior probability images for
% grey, white and CSF.
% [defaults.segment.estimate.samp] - distance to sample between points
% duting the EM algorithm for segmenting and bias correction.
% [defaults.segment.estimate.bb]
% The bounding box enclosing the region to use for the EM algorithm.
% [defaults.segment.estimate.affreg.smosrc]
% Amount of smoothing to use for the source image during the affine
% registration step.
% [defaults.segment.estimate.affreg.regtype] - specifies the way that affine
% registration is regularised (see spm_affreg.m).
% [defaults.segment.estimate.affreg.weight] -  Template weighting image
% applied during affine registration.
%
%_______________________________________________________________________
% @(#)spm_segment_ui.m	2.4 John Ashburner 03/02/17

global defaults
if nargin==1 & strcmp(lower(P),'defaults');
	defaults.segment = edit_defaults(defaults.segment);
	return;
end;
segment_ui(defaults.segment);
return;
%=======================================================================
 
%=======================================================================
function segment_ui(flags)

SPMid = spm('FnBanner',mfilename,'2.4');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Segment');
spm_help('!ContextHelp',mfilename);

don = 0;
for i = 1:1000,
	ok = 0;
	while ~ok,
		PF = spm_get(Inf,'IMAGE',...
			['Image(s), subj ' num2str(i)]);
		if isempty(PF), don = 1; break; end;
		VF{i} = spm_vol(PF);
		ok    = dims_ok(VF{i});
	end;
	if don, break; end;
end;

if spm_input('Already spatially normalised?', 1, 'y/n') == 'n',
	% Get template
	%-----------------------------------------------------------------------
	templates = str2mat(	fullfile(spm('Dir'),'templates','T1.mnc'),...
				fullfile(spm('Dir'),'templates','T2.mnc'),...
				fullfile(spm('Dir'),'templates','PD.mnc'),...
				fullfile(spm('Dir'),'templates','EPI.mnc'));

	% Get modality of target
	respt = spm_input('Modality?','+1','m',...
		['modality - T1 MRI|modality - T2 MRI|'...
		 'modality - PD MRI|modality - EPI MR|--other--'],...
		[1 2 3 4 0],1);
	if respt > 0,
		PG = deblank(templates(respt,:));
		VG = spm_vol(PG);
	else,
		ok = 0;
		while ~ok,
			PG = spm_get(Inf,'IMAGE',['Template(s) for affine matching'],...
				fullfile(spm('Dir'),'templates'));
			VG = spm_vol(PG);
			ok = dims_ok(VG);
		end;
	end;
else,
	VG = eye(4);
end;

spm('Pointer','Watch');
for i = 1:length(VF),
	spm('FigName',['Segment: working on subj ' num2str(i)],Finter,CmdLine);
	spm_segment(VF{i},VG,flags);
end;

spm_figure('Clear','Interactive');
spm('FigName','Segment: done',Finter,CmdLine);
spm('Pointer');
return;
%=======================================================================
 
%=======================================================================
function ok = dims_ok(vv)
ok = 0;
if isempty(vv), return; end;
if prod(size(vv))==1,
	ok = 1;
else,
	tmp1 = cat(1,vv.dim);
	tmp2 = cat(3,vv.mat);
	if ~any(any(diff(tmp1(:,1:3)))) & ~any(any(any(diff(tmp2,1,3)))),
		ok=1;
	end;
end;
return;
%=======================================================================
 
%=======================================================================
function defs = edit_defaults(defs)

rg  = [0 0.00001 0.0001 0.001 0.01 0.1 1.0 10];
tmp = find(rg == defs.estimate.reg);
if isempty(tmp), tmp = 4; end;
defs.estimate.reg = spm_input('Bias regularisation?','+1','m',...
	['no regularisation (0)|extremely light regularisation (0.00001)|'...
	 'very light regularisation (0.0001)|light regularisation (0.001)|',...
	 'medium regularisation (0.01)|heavy regularisation (0.1)|'...
	 'very heavy regularisation (1)|extremely heavy regularisation (10)'],...
	rg, tmp);
 
co  = [20 25 30 35 40 45 50 60 70 80 90 100 Inf];
if isinf(defs.estimate.cutoff),
	% not sure if Inf==Inf holds for all Matlab versions
	tmp = length(co);
else,
	tmp = find(co == defs.estimate.cutoff);
	if isempty(tmp), tmp = 4; end;
end;
defs.estimate.cutoff =  spm_input('Bias cutoff?','+1','m',...
	[' 20mm cutoff| 25mm cutoff| 30mm cutoff| 35mm cutoff| 40mm cutoff|'...
	 ' 45mm cutoff| 50mm cutoff| 60mm cutoff| 70mm cutoff| 80mm cutoff|'...
	 ' 90mm cutoff|100mm cutoff|No correction'],...
	co, tmp);

tmp = find(defs.write.cleanup == [1 0]);
if isempty(tmp), tmp = 1; end;
defs.write.cleanup =  spm_input('Clean up the partitions?','+1','m',...
	['Clean up the partitions|Dont do cleanup'],[1 0],tmp); 

tmp = find(defs.write.wrt_cor == [1 0]);
if isempty(tmp), tmp = 1; end;
defs.write.wrt_cor = spm_input('Write bias correted image?','+1','m',...
	['Write bias corrected|Dont write bias corrected'],[1 0],tmp);
%=======================================================================
 
%=======================================================================
