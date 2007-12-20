function opts = spm_config_coreg
% Configuration file for coregister jobs
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_config_coreg.m 1032 2007-12-20 14:45:55Z john $

ref.type = 'files';
ref.name = 'Reference Image';
ref.tag  = 'ref';
ref.filter = 'image';
ref.num  = 1;
ref.help = {[...
'This is the image that is assumed to remain stationary (sometimes ',...
'known as the target or template image), while the source image ',...
'is moved to match it.']};

%------------------------------------------------------------------------

source.type = 'files';
source.name = 'Source Image';
source.tag  = 'source';
source.filter = 'image';
source.num  = 1;
source.help = {...
'This is the image that is jiggled about to best match the reference.'};

%------------------------------------------------------------------------

other.type = 'files';
other.name = 'Other Images';
other.tag  = 'other';
other.filter = 'image';
other.num  = [0 Inf];
other.val  = {''};
other.help = {[...
'These are any images that need to remain in alignment with the ',...
'source image.']};

%------------------------------------------------------------------------

cost_fun.type = 'menu';
cost_fun.name = 'Objective Function';
cost_fun.tag  = 'cost_fun';
cost_fun.labels = {'Mutual Information','Normalised Mutual Information',...
'Entropy Correlation Coefficient','Normalised Cross Correlation'};
cost_fun.values = {'mi','nmi','ecc','ncc'};
cost_fun.def  = 'coreg.estimate.cost_fun';
cost_fun.help = {[...
'Registration involves finding parameters that either maximise or ',...
'minimise some objective function. ',...
'For inter-modal registration, use Mutual Information\* \cite{collignon95,wells96}*/, Normalised Mutual Information/* \cite{studholme99}*/, or ',...
'Entropy Correlation Coefficient/* \cite{maes97}*/.',...
'For within modality, you could also use Normalised Cross Correlation.']};

%------------------------------------------------------------------------

sep.type = 'entry';
sep.name = 'Separation';
sep.tag  = 'sep';
sep.num  = [1 Inf];
sep.strtype = 'e';
sep.def  = 'coreg.estimate.sep';
sep.help = {[...
'The average distance between sampled points (in mm).  Can be a vector ',...
'to allow a coarse registration followed by increasingly fine ones.']};

%------------------------------------------------------------------------

tol.type = 'entry';
tol.name = 'Tolerances';
tol.tag = 'tol';
tol.num = [1 12];
tol.strtype = 'e';
tol.def = 'coreg.estimate.tol';
tol.help = {[...
'The accuracy for each parameter.  Iterations stop when differences ',...
'between successive estimates are less than the required tolerance.']};

%------------------------------------------------------------------------

fwhm.type = 'entry';
fwhm.name = 'Histogram Smoothing';
fwhm.tag  = 'fwhm';
fwhm.num  = [1 2];
fwhm.strtype = 'e';
fwhm.def = 'coreg.estimate.fwhm';
fwhm.help = {[...
'Gaussian smoothing to apply to the 256x256 joint histogram. Other ',...
'information theoretic coregistration methods use fewer bins, but ',...
'Gaussian smoothing seems to be more elegant.']};

%------------------------------------------------------------------------

eoptions.type = 'branch';
eoptions.name = 'Estimation Options';
eoptions.tag  = 'eoptions';
eoptions.val  = {cost_fun,sep,tol,fwhm};
eoptions.help = {...
['Various registration options, which are passed to the ',...
'Powell optimisation algorithm/* \cite{press92}*/.']};

%------------------------------------------------------------------------

est.type = 'branch';
est.name = 'Coreg: Estimate';
est.tag  = 'estimate';
est.val  = {ref,source,other,eoptions};
est.prog = @estimate;
p1 = [...
'The registration method used here is based on work by Collignon et al/* \cite{collignon95}*/. ',...
'The original interpolation method described in this paper has been ',...
'changed in order to give a smoother cost function.  The images are ',...
'also smoothed slightly, as is the histogram.  This is all in order to ',...
'make the cost function as smooth as possible, to give faster ',...
'convergence and less chance of local minima.'];
p2 = [...
'At the end of coregistration, the voxel-to-voxel affine transformation ',...
'matrix is displayed, along with the histograms for the images in the ',...
'original orientations, and the final orientations.  The registered ',...
'images are displayed at the bottom.'];
p3 = [...
'Registration parameters are stored in the headers of the "source" ',...
'and the "other" images.'];

est.help = {p1,'',p2,'',p3};

%------------------------------------------------------------------------

interp.type = 'menu';
interp.name = 'Interpolation';
interp.tag  = 'interp';
interp.labels = {'Nearest neighbour','Trilinear','2nd Degree B-Spline',...
'3rd Degree B-Spline','4th Degree B-Spline','5th Degree B-Spline',...
'6th Degree B-Spline','7th Degree B-Spline'};
interp.values = {0,1,2,3,4,5,6,7};
interp.def  = 'coreg.write.interp';
interp.help = {[...
'The method by which the images are sampled when being written in a ',...
'different space. ',...
'Nearest Neighbour is fastest, but not normally recommended. ',...
'It can be useful for re-orienting images while preserving the original ',...
'intensities (e.g. an image consisting of labels). ',...
'Bilinear Interpolation is OK for PET, or realigned and re-sliced fMRI. ',...
'If subject movement (from an fMRI time series) is included in the transformations ',...
'then it may be better to use a higher degree approach. ',...
'Note that higher degree B-spline interpolation/* \cite{thevenaz00a,unser93a,unser93b}*/ is slower because it uses more ',...
'neighbours.']};

%------------------------------------------------------------------------

wrap.type = 'menu';
wrap.name = 'Wrapping';
wrap.tag  = 'wrap';
wrap.labels = {'No wrap','Wrap X','Wrap Y','Wrap X & Y','Wrap Z',...
'Wrap X & Z','Wrap Y & Z','Wrap X, Y & Z'};
wrap.values = {[0 0 0],[1 0 0],[0 1 0],[1 1 0],[0 0 1],[1 0 1],[0 1 1],[1 1 1]};
wrap.def    = 'coreg.write.wrap';
wrap.help = {...
'These are typically:',...
['    No wrapping - for PET or images that have already',...
 '                  been spatially transformed.'],...
['    Wrap in  Y  - for (un-resliced) MRI where phase encoding',...
 '                  is in the Y direction (voxel space).']};

%------------------------------------------------------------------------

mask.type = 'menu';
mask.name = 'Masking';
mask.tag  = 'mask';
mask.labels = {'Mask images','Dont mask images'};
mask.values = {1,0};
mask.def    = 'coreg.write.mask';
mask.help = {[...
'Because of subject motion, different images are likely to have different ',...
'patterns of zeros from where it was not possible to sample data. ',...
'With masking enabled, the program searches through the whole time series ',...
'looking for voxels which need to be sampled from outside the original ',...
'images. Where this occurs, that voxel is set to zero for the whole set ',...
'of images (unless the image format can represent NaN, in which case ',...
'NaNs are used where possible).']};

%------------------------------------------------------------------------

roptions.type = 'branch';
roptions.name = 'Reslice Options';
roptions.tag  = 'roptions';
roptions.val  = {interp,wrap,mask};
roptions.help = {'Various reslicing options.'};

%------------------------------------------------------------------------

estwrite.type = 'branch';
estwrite.name = 'Coreg: Estimate & Reslice';
estwrite.tag  = 'estwrite';
estwrite.val  = {ref,source,other,eoptions,roptions};
estwrite.prog = @estimate_reslice;
estwrite.vfiles = @vfiles_estwrite;
p1 = [...
'The registration method used here is based on work by Collignon et al/* \cite{collignon95}*/. ',...
'The original interpolation method described in this paper has been ',...
'changed in order to give a smoother cost function.  The images are ',...
'also smoothed slightly, as is the histogram.  This is all in order to ',...
'make the cost function as smooth as possible, to give faster ',...
'convergence and less chance of local minima.'];
p2 = [...
'At the end of coregistration, the voxel-to-voxel affine transformation ',...
'matrix is displayed, along with the histograms for the images in the ',...
'original orientations, and the final orientations.  The registered ',...
'images are displayed at the bottom.'];
p3 = [...
'Registration parameters are stored in the headers of the "source" ',...
'and the "other" images. These images are also resliced to match the ',...
'source image voxel-for-voxel. The resliced images are named the same as ',...
'the originals except that they are prefixed by ''r''.'];
estwrite.help = {p1,'',p2,'',p3};

%------------------------------------------------------------------------

ref.type = 'files';
ref.name = 'Image Defining Space';
ref.tag  = 'ref';
ref.filter = 'image';
ref.num  = 1;
ref.help = {[...
'This is analogous to the reference image.  Images are resliced to match ',...
'this image (providing they have been coregistered first).']};

%------------------------------------------------------------------------
 
source.type = 'files';
source.name = 'Images to Reslice';
source.tag  = 'source';
source.filter = 'image';
source.num  = Inf;
source.help = {[...
'These images are resliced to the same dimensions, voxel sizes, ',...
'orientation etc as the space defining image.']};

%------------------------------------------------------------------------

write.type = 'branch';
write.name = 'Coreg: Reslice';
write.tag  = 'write';
write.val  = {ref,source,roptions};
write.prog = @reslice;
write.vfiles = @vfiles_write;
write.help = {[...
'Reslice images to match voxel-for-voxel with an image defining ',...
'some space. The resliced images are named the same as the originals ',...
'except that they are prefixed by ''r''.']};

%------------------------------------------------------------------------

opts.type = 'repeat';
opts.name = 'Coreg';
opts.tag  = 'coreg';
opts.values = {est,write,estwrite};
opts.num  = [1 Inf];
opts.modality = {'PET','FMRI','VBM'};
p1 = [...
'Within-subject registration using a rigid-body model. ',...
'A rigid-body transformation (in 3D) can be parameterised by three ',...
'translations and three rotations about the different axes.'];
p2 = [...
'You get the options of estimating the transformation, reslicing images ',...
'according to some rigid-body transformations, or estimating and ',...
'applying rigid-body transformations.'];
opts.help = {p1,'',p2};

return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function estimate(varargin)
job = varargin{1};
%disp(job);
%disp(job.eoptions);

x  = spm_coreg(strvcat(job.ref), strvcat(job.source),job.eoptions);
M  = inv(spm_matrix(x));
PO = strvcat(strvcat(job.source),strvcat(job.other));
MM = zeros(4,4,size(PO,1));
for j=1:size(PO,1),
	MM(:,:,j) = spm_get_space(deblank(PO(j,:)));
end;
for j=1:size(PO,1),
	spm_get_space(deblank(PO(j,:)), M*MM(:,:,j));
end;

return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function reslice(varargin)
job = varargin{1};

P            = strvcat(strvcat(job.ref),strvcat(job.source));
flags.mask   = job.roptions.mask;
flags.mean   = 0;
flags.interp = job.roptions.interp;
flags.which  = 1;
flags.wrap   = job.roptions.wrap;

spm_reslice(P,flags);

return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function estimate_reslice(varargin)
job = varargin{1};
%disp(job);
%disp(job.eoptions);
%disp(job.roptions);

job.ref    = strvcat(job.ref);
job.source = strvcat(job.source);
job.other  = strvcat(job.other);

x  = spm_coreg(job.ref, job.source,job.eoptions);
M  = inv(spm_matrix(x));
PO = strvcat(job.source,job.other);
MM = zeros(4,4,size(PO,1));
for j=1:size(PO,1),
        MM(:,:,j) = spm_get_space(deblank(PO(j,:)));
end;
for j=1:size(PO,1),
        spm_get_space(deblank(PO(j,:)), M*MM(:,:,j));
end;

P            = strvcat(job.ref,job.source,job.other);
flags.mask   = job.roptions.mask;
flags.mean   = 0;
flags.interp = job.roptions.interp;
flags.which  = 1;
flags.wrap   = job.roptions.wrap;

spm_reslice(P,flags);

return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function vf = vfiles_write(varargin)
job = varargin{1};
vf  = cell(size(job.source));
for i=1:numel(job.source),
    [pth,nam,ext,num] = spm_fileparts(job.source{i});
    vf{i} = fullfile(pth,['r', nam, ext, num]);
end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function vf = vfiles_estwrite(varargin)
job = varargin{1};
P   = {job.source{:},job.other{:}};
vf  = cell(size(P));
for i=1:numel(P),
    [pth,nam,ext,num] = spm_fileparts(P{i});
    vf{i} = fullfile(pth,['r', nam, ext, num]);
end;


