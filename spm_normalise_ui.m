function spm_normalise_ui(opt)
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
% parameters are saved in the "*_sn.mat" file.
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
%                        The Prompts Explained
%_______________________________________________________________________
%
% 'Which option?'
% 	'Determine Parameters Only'
% 	Computes the parameters which best fit the image to the
% 	template(s), and saves them to the file imagename'_sn.mat'.
%
% 	'Write Normalised Only'
% 	Select the appropriate '_sn.mat' file, and the images you
% 	wish to have normalised. The normalised images will be written
% 	prefixed with 'r'.
%
% 	'Determine Parameters & Write Normalised'
% 	Combines the above two steps into one.
%
% Options for Determine Parameters:
%
% 'Template image(s) '
% Select one or more templates. The image(s) will be fitted (in the least 
% squares sense) to the optimum linear combination of these templates.
%
% 'Source image, subj #'
% This is the image that is matched to the template by minimising the
% (weighted) mean squared difference between them.
%
% If source weighting is enabled, then:
% 'Source weighting image (or Done for none)'
% Applies a weighting image to the source image during the parameter
% estimation.  Weights are as for template weighting (0-1).  Used
% (usually) to prevent unusual or abnormal areas of brain (e.g.
% stroke, tumour) influencing normalisation to normal brain.
%
% Options for Write Normalised:
%
% 'Parameters (or Done), subj #'
% Select the '_sn.mat' file containing the spatial normalisation
% parameters for that subject. If nothing is selected, then the routine
% will assume that no more subjects need to be selected.
%
% 'Images to write, subj #'
% These are the images for warping according to the estimated parameters. 
%
%_______________________________________________________________________
%
%                           Defaults Options
%_______________________________________________________________________
%[   things in square brackets indicate corresponding defaults field   ]
%
% Defaults for Parameter Estimation:
%
% 'Weight template when registering?'
%       'No Weighting'       - no template weighting
%       'Default Brain Mask' - weighting with .../apriori/brainmask.mnc
%       'Specified Weighting' - user selected weighting file
% Applies a weighting mask to the template(s) during the parameter
% estimation.  With the default brain mask, weights in and around the
% brain have values of one whereas those clearly outside the brain are
% zero.  This is an attempt to base the normalization purely upon
% the shape of the brain, rather than the shape of the head (since
% low frequency basis functions can not really cope with variations
% in skull thickness).
% The option is now available for a user specified weighting image.
% This should have the same dimensions and mat file as the template
% images, with values in the range of zero to one.
% [defaults.normalise.estimate.weight]
%
% 'Weight source images when registering?'
% 	'Weight sources'      - enables source weighting
% 	'Dont weight sources' - disables source weighting
% Applies a weighting mask to the object image(s) during the parameter
% estimation.  Weights are as for the template mask (0-1).  Used
% (usually) to prevent unusual or abnormal areas of brain (e.g.
% stroke, tumour) influencing normalisation of normal brain.
% [defaults.normalise.estimate.wtsrc]
%
% 'Nonlinear Regularization'
% Allows the amount of regularisation to be changed.
% Pick a medium value.  However, if your normalized images appear
% distorted, then it may be an idea to increase the amount of
% regularization - or even just use an affine normalization.
% The regularization influences the smoothness of the deformation
% fields.
% [defaults.normalise.estimate.reg] 
%
% 'Cutoff'
% Cutoff of DCT bases.  Only DCT bases of periods longer than the
% cutoff are used to describe the warps. The number used will
% depend on the cutoff and the field of view of the template image(s).
% [defaults.normalise.estimate.cutoff]
%
% '# Nonlinear Iterations?'
% Number of iterations of nonlinear warping performed.
% [defaults.normalise.estimate.nits]
%
%
% Defaults for Writing Normalised:
%
% 'Preserve what?'
% 	'Preserve Concentrations'
% 	Spatially normalised images are not "modulated". The warped images
% 	preserve the intensities of the original images.
% 	'Preserve Total'
% 	Spatially normalised images are "modulated" in order to preserve the
% 	total amount of signal in the images. Areas that are expanded during
% 	warping are correspondingly reduced in intensity.
% [defaults.normalise.write.preserve]
%
% 'Bounding Box?'
% The bounding box (in mm) of the volume which is to be written
% (relative to the anterior commissure).
% [defaults.normalise.write.bb]
%
% 'Voxel Sizes '
% The voxel sizes (x, y & z, in mm) of the written normalised images.
% [defaults.normalise.write.vox]
%
% 'Interpolation Method?'
% The method by which the images are sampled when being written in a
% different space.
% 	'Nearest Neighbour'
% 		- Fastest, but not normally recommended.
% 	'Bilinear Interpolation'
% 		- OK for PET, or realigned fMRI.
%
% 	'B-spline Interpolation'
%		- Better quality (but slower) interpolation, especially
% 		  with higher degree splines.  Don't use B-splines when
% 		  there is any region of NaN or Inf in the images.
% [defaults.normalise.write.interp]
%
% These are typically:
%       'No wrapping' - for PET or images that have already
%                       been spatially transformed.
%       'Wrap in  Y'  - for (un-resliced) MRI where phase encoding
%                       is in the Y direction (voxel space).
% [defaults.normalise.write.wrap]
%
%
% Other Defaults are useful for making spm_normalise.m more flexible:
% [defaults.normalise.estimate.smoref] - additional smoothing of
% template image.
% [defaults.normalise.estimate.regtype] - specifies the way that affine
% registration is regularised (see spm_affreg.m).
%
%_______________________________________________________________________
% %W% John Ashburner %E%

% Programmers notes
%-----------------------------------------------------------------------
%
% FORMAT spm_normalise('Defaults')
% acts as a user interface for setting normalisation defaults.
%
%-----------------------------------------------------------------------

global defaults
defs = defaults.normalise;

if nargin==0 | strcmp(lower(opt),'ui'),
	run_ui(defs);
elseif nargin>0 & strcmp(lower(opt),'defaults'),
	defaults.normalise = get_defs(defs);
end;
return;
%_______________________________________________________________________
 
%_______________________________________________________________________
function run_ui(defs)
% Run spatial normalisation ui
% FORMAT run_ui(defs)

SPMid = spm('FnBanner',mfilename,'%I%');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Normalise');
spm_help('!ContextHelp',mfilename);

a1 = spm_input('Which option?',1,'m',...
        ['Determine Parameters Only|'...
        'Write Normalised Only|'...
        'Determine Parameters & Write Normalised'],...
        [1 2 3],3);

if a1 == 1 | a1 == 3,
        % Get template(s)
        ok = 0;
        while ~ok,
                Template = spm_get(Inf,'IMAGE',['Template image(s)'],...
                        fullfile(spm('Dir'),'templates'));
                vv = spm_vol(Template);
                if prod(size(vv))==1,
                        ok = 1;
                elseif prod(size(vv)) ~= 0,
                        tmp1 = cat(1,vv.dim);
                        tmp2 = reshape(cat(3,vv.mat),4*4,prod(size(vv)));
                        if ~any(any(diff(tmp1(:,1:3)))) &...
                           ~any(any(diff(tmp2,1,2))),
                                ok=1;
                        end;
		end;
	end;
end;

% Select images..
%-----------------------------------------------------------------------
for i=1:1000,
	if a1 == 1 | a1 == 3,

		P = spm_get([0,1],'IMAGE',['Source image, subj ' num2str(i)]);
		if isempty(P), break; end;

		subj(i).P = P;
		% source weight
		if defs.estimate.wtsrc,
			subj(i).objmask = spm_get([0,1],'IMAGE',...
			['Source weighting image (or Done for none)']);
		else,
			subj(i).objmask = '';
		end;
		subj(i).matname = [spm_str_manip(subj(i).P,'sd') '_sn.mat'];
	else,
		matname = spm_get([0,1],'_sn.mat',['Parameters (or Done), subj ' num2str(i)]);
		if isempty(matname), break; end;
		subj(i).matname = matname;
	end;
	if a1 == 2 | a1 == 3,
		subj(i).PP = spm_get(Inf,'IMAGE',['Images to write, subj ' num2str(i)]);
	end;
end;
if ~exist('subj'),
        spm_figure('Clear','Interactive');
        return;
end;

% Go and do the work
%-----------------------------------------------------------------------
spm('Pointer','Watch')
if a1 == 1 | a1 == 3,
	for i=1:length(subj),
		spm('FigName',['Normalising (est) subj ' num2str(i)],...
			Finter,CmdLine);
		spm_normalise(Template, subj(i).P, subj(i).matname,...
			defs.estimate.weight, subj(i).objmask,defs.estimate);
	end;
end;

if a1 == 2 | a1 == 3,
	for i=1:length(subj),
		spm('FigName',['Normalising (write) subj ' num2str(i)],...
			Finter,CmdLine);
		spm_write_sn(subj(i).PP,subj(i).matname,defs.write);
	end;
end;

fprintf('\n\n');
spm('FigName','Normalise: done',Finter,CmdLine);
spm('Pointer');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function defs = get_defs(defs)
% Change spatial normalisation defaults
% FORMAT defs = get_defs(defs)
if spm_input(['Defaults for..?'],1,'m',...
	['Defaults for Parameter Estimation|'...
	 'Defaults for Writing Normalised'],[1 0]),
	defs.estimate.weight = get_weight(defs.estimate.weight);
	defs.estimate.wtsrc  = get_wtsrc(defs.estimate.wtsrc);
	defs.estimate.cutoff = get_cutoff(defs.estimate.cutoff);
	if ~isinf(defs.estimate.cutoff),
		defs.estimate.reg    = get_reg(defs.estimate.reg);
		defs.estimate.nits   = get_nits(defs.estimate.nits);
	end;
else,
	defs.write.preserve = get_preserve(defs.write.preserve);
	defs.write.bb       = get_bb(defs.write.bb);
	defs.write.vox      = get_vox(defs.write.vox);
	defs.write.interp   = get_interp(defs.write.interp);
	defs.write.wrap     = get_wrap(defs.write.wrap);
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function wtsrc = get_wtsrc(wtsrc)
% ask for source image weighting
% FORMAT wtsrc = get_wtsrc(wtsrc)

wtsrc = spm_input('Weight source images when registering?',...
	'+1', 'm','Dont weight sources|Weight sources',[0 1],...
	find([0 1] == wtsrc));
%_______________________________________________________________________

%_______________________________________________________________________
function nits = get_nits(nits)
% Get number of nonlinear iterations
% FORMAT nits = get_nits(nits)

if prod(nits) > 0,
        tmp2 = [1 3 5 8 12 16];
        tmp = find(tmp2 == nits);
        if isempty(tmp) tmp = length(tmp2); end;
        nits = spm_input(['# Nonlinear Iterations?'],'+1','m',...
        	['1  nonlinear iteration |3  nonlinear iterations'...
        	'|5  nonlinear iterations|8  nonlinear iterations'...
        	'|12 nonlinear iterations|16 nonlinear iterations'],tmp2, tmp);
else, nits = 0; end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function weight = get_weight(weight)
% Get an image to weight the registration with
% FORMAT weight = get_weight(weight)

def_brainmask = fullfile(spm('Dir'),'apriori','brainmask.img');
tmp = ~isempty(weight);
if tmp, tmp = tmp + 1 - strcmp(weight,def_brainmask); end;
tmp = spm_input('Weight template when registering?', '+1', 'm',...
	'No Weighting|Default Brain Mask|Specified Weighting',[0 1 2],...
		tmp+1);
if ~tmp,   weight = ''; end;
if tmp==1, weight = def_brainmask; end;
if tmp==2, weight = spm_get(1,'IMAGE','Specify Weighting Image'); end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function cutoff = get_cutoff(cutoff)
% Get cutoff frequency of DCT bases
 
tmp2 = [15 20 25 30 35 40 45 50 60 70 80 Inf];
tmp = find(tmp2 == cutoff);
if isempty(tmp) tmp = length(tmp2); end;
cutoff = spm_input('Cutoff','+1','m',...
        ['15mm cutoff|20mm cutoff|25mm cutoff|30mm cutoff|'...
	 '35mm cutoff|40mm cutoff|45mm cutoff|50mm cutoff|'...
	 '60mm cutoff|70mm cutoff|80mm cutoff|Affine only'], tmp2, tmp);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function reg = get_reg(reg)
% Get amount of regularisation
% FORMAT reg = get_reg(reg)

tmp2 = [100 10 1 0.1 0.01];
tmp = find(tmp2 == reg);
if isempty(tmp) tmp = length(tmp2); end;
reg = spm_input('Nonlinear Regularization','+1','m',...
        ['Extremely Heavy regularization (100)|Heavy regularization (10)|'...
         'Medium regularization (1)|Light regularization (0.1)|'...
         'Very Light regularization (0.01)'], tmp2, tmp);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function vx = get_vox(vx)
% Allow bounding boxes to be modified
% format vx = get_vox(vx)

voxdims    = [ 1   1   1 ; 1.5 1.5 1.5 ; 2   2   2 ; 3   3   3 ; 4   4   4 ; 1   1   2 ; 2   2   4];
voxprompts = ' 1   1   1 | 1.5 1.5 1.5 | 2   2   2 | 3   3   3 | 4   4   4 | 1   1   2 | 2   2   4';
if prod(size(vx)) == 3,
        tmp = find(voxdims(:,1)==vx(1) & voxdims(:,2)==vx(2) & voxdims(:,3)==vx(3));
        if isempty(tmp), tmp = size(voxdims,1)+1; end;
else,
        tmp = size(voxdims,1)+2;
end;

tmp = spm_input(...
        ['Voxel Sizes?'], '+1','m', [ voxprompts '|Customise'],...
        [1:size(voxdims,1) 0], tmp);
if tmp>0, vx = voxdims(tmp,:);
elseif tmp == 0,
        Vox = [];
        if (prod(size(vx)) ~= 3) vx = [2 2 2]; end
                vx = spm_input('Voxel Sizes ','+0', 'e', ...
                        sprintf('%d %d %d', vx(1), vx(2), vx(3)), 3)';
                vx = reshape(vx,1,3);
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function preserve = get_preserve(preserve)
% Preserve amounts or concentrations during warping.

preserve = spm_input('Preserve what?','+1','m',...
        ['Preserve Concentrations|Preserve Total'],...
         [0 1], preserve+1);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function bb = get_bb(bb)
% Allow bounding box to be modified.
% FORMAT bb = get_bb(bb)

bboxes  = [   -78 78 -112 76  -50 85
	      -64 64 -104 68  -28 72
	      -90 91 -126 91  -72 109
              -95 95 -112 76 -50 95];
bbprompt =  [' -78:78 -112:76  -50:85  (Default)|'...
	     ' -64:64 -104:68  -28:72  (SPM95)   |'...
	     ' -90:91 -126:91  -72:109 (Template)|'...
             ' -95:95 -112:76  -50:95 '];
if prod(size(bb)) == 6,
        tmp = find(     bb(1) == bboxes(:,1) & bb(2) == bboxes(:,2) & ...
                        bb(3) == bboxes(:,3) & bb(4) == bboxes(:,4) & ...
                        bb(5) == bboxes(:,5) & bb(6) == bboxes(:,6));
        if isempty(tmp), tmp = size(bboxes,1)+1; end;
else,
        tmp = size(bboxes,1)+2;
        bb  = reshape(bboxes(1,:),2,3);
end;

tmp = spm_input('Bounding Box?','+1','m',...
        [ bbprompt '|Customise'], [1:size(bboxes,1) 0],tmp);
if tmp>0,
        bb =reshape(bboxes(tmp,:),2,3);
elseif tmp == 0,
        if prod(size(bb)) ~= 6, bb = reshape(bboxes(1,:),2,3); end;
        directions = 'XYZ';
        nbb = zeros(2,3);
        for d=1:3,
                str = sprintf('%d %d', bb(1,d), bb(2,d));
                nbb(:,d) = spm_input(['Bounding Box ' directions(d) ],....
                        '+1', 'e',str, 2);
        end;
	bb = nbb;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function interp = get_interp(interp)
% Get interpolation method (for writing images)
% FORMAT interp = get_interp(interp)

interp = spm_input('Interpolation Method?','+1','m',...
	['Nearest Neighbour|Trilinear Interpolation|',...
	 '2nd Degree B-spline|3rd Degree B-spline|4th Degree B-spline|',...
	 '5th Degree B-spline|6th Degree B-spline|7th Degree B-spline'],...
	 [0 1 2 3 4 5 6 7], interp+1);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function wrap = get_wrap(wrap)
% Get image wrapping information
% FORMAT wrap = get_wrap(wrap)
wraps = [0 0 0 ; 1 0 0; 0 1 0; 1 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 1];
t     = find(all(repmat(wrap(:)',8,1) == wraps, 2));
if isempty(t), t = 1; end;
p     = spm_input('Way to wrap images?','+1','m',...
	['No wrap|Wrap X|Wrap Y|Wrap X & Y|Wrap Z|Wrap X & Z|Wrap Y & Z|Wrap X, Y & Z'],...
	[1 2 3 4 5 6 7 8], t);
wrap = wraps(p,:);
return;

