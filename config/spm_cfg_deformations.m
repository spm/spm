function conf = spm_cfg_deformations
% Configuration file for deformation jobs
%_______________________________________________________________________
% Copyright (C) 2008-2016 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_cfg_deformations.m 7700 2019-11-21 17:09:15Z john $

hsummary = {
'Utility for working with deformation fields.',...
['They can be loaded, inverted, combined etc, and the results ',...
'either saved to disk, or applied to some image or surface file. ',...
'This utility was intended for imaging experts and may therefore ',...
'be a bit difficult for naive users. ',...
'It provides a great deal of flexibility, which may be confusing to some.']};

hinv = {[...
'Creates the inverse of a deformation field. ',...
'Deformations are assumed to be one-to-one, in which case they ',...
'have a unique inverse.  If y'':A->B is the inverse of y:B->A, then ',...
'y'' o y = y o y'' = Id, where Id is the identity transform.'],...
'',...
'Deformations are inverted using the method described in the appendix of:',...
['    * Ashburner J, Andersson JLR & Friston KJ (2000) ',...
 '"Image Registration using a Symmetric Prior - in Three-Dimensions." ',...
 'Human Brain Mapping 9(4):212-225']};

hcomp = {[...
'Deformation fields can be thought of as mappings. ',...
'These can be combined by the operation of "composition", which is ',...
'usually denoted by a circle "o". ',...
'Suppose x:A->B and y:B->C are two mappings, where A, B and C refer ',...
'to domains in 3 dimensions. ',...
'Each element a in A points to element x(a) in B. ',...
'This in turn points to element y(x(a)) in C, so we have a mapping ',...
'from A to C. ',...
'The composition of these mappings is denoted by yox:A->C. ',...
'Compositions can be combined in an associative way, such that zo(yox) = (zoy)ox.'],...
'',[...
'In this utility, the right-to-left order of the compositions is ',...
'from top to bottom (note that the rightmost deformation would ',...
'actually be applied first).']};

hsn = {[...
'Spatial normalisation, and the unified segmentation model of ',...
'SPM5 save a parameterisation of deformation fields.  These consist ',...
'of a combination of an affine transform, and nonlinear warps that ',...
'are parameterised by a linear combination of cosine transform ',...
'basis functions.  These are saved in *_sn.mat files, which can be ',...
'converted to deformation fields.']};

hvox = {[...
'Specify the voxel sizes of the deformation field to be produced. ',...
'Non-finite values will default to the voxel sizes of the template image',...
'that was originally used to estimate the deformation.']};

hbb = {[...
'Specify the bounding box of the deformation field to be produced. ',...
'Non-finite values will default to the bounding box of the template image',...
'that was originally used to estimate the deformation.']};

himgr = {[...
'Deformations can be thought of as vector fields, and represented ',...
'by three-volume images.  In SPM, deformation fields are saved in ',...
'NIfTI format, with dimensions xdim x ydim x zdim x 1 x 3. ',...
'Each voxel contains the x, y and z mm coordinates of where the deformation points.']};

himgw = {[...
'Save the result as a three-volume image.  "y_" will be prepended to the ',...
'filename.']};

hdetw = {[...
'Save the Jacobian determinants as an image.  "j_" will be prepended to the ',...
'filename.']};

happly = {[...
'Apply the resulting deformation field to some images. ',...
'The filenames will be prepended by "w".']};

hmatname = {...
'Specify the _sn.mat to be used.'};

himg = {...
'Specify the image file on which to base the dimensions, orientation etc.'};

hid = {[...
'This option generates an identity transform, but this can be useful for ',...
'changing the dimensions of the resulting deformation (and any images that ',...
'are generated from it).  Dimensions, orientation etc are derived from ',...
'an image.']};

hidbbvox = {[...
'This option generates an identity transform, but this can be useful for ',...
'changing the dimensions of the resulting deformation (and any images that ',...
'are generated from it).  Dimensions, orientation etc are derived from ',...
'a specified bounding box and voxel dimensions.']};

def          = files('Deformation Field','def','.*y_.*\.nii$',[1 1]);
def.help     = himgr;

matname      = files('Parameter File','matname','.*_sn\.mat$',[1 1]);
matname.help = hmatname;

vox          = entry('Voxel sizes','vox','r',[1 3]);
vox.val      = {[NaN NaN NaN]};
vox.help     = hvox;

bb           = entry('Bounding box','bb','r',[2 3]);
bb.val       = {[NaN NaN NaN;NaN NaN NaN]};
bb.help      = hbb;

sn2def       = branch('Imported _sn.mat','sn2def',{matname,vox,bb});
sn2def.help  = hsn;

img          = files('Image to base Id on','space','nifti',[1 1]);
img.help     = himg;

id           = branch('Identity (Reference Image)','id',{img});
id.help      = hid;

idbbvox      = branch('Identity (Bounding Box and Voxel Size)','idbbvox',{vox, bb});
idbbvox.help = hidbbvox;

ffield = files('Flow field','flowfield','nifti',[1 1]);
ffield.ufilter = '^u_.*';
ffield.help = {...
    ['The flow field stores the deformation information. '...
     'The same field can be used for both forward or backward deformations '...
     '(or even, in principle, half way or exaggerated deformations).']};
%------------------------------------------------------------------------
forbak = mnu('Forward/Backwards','times',{'Backward','Forward'},{[1 0],[0 1]});
forbak.val  = {[1 0]};
forbak.help = {[...
    'The direction of the Dartel flow.  '...
    'Note that a backward transform will warp an individual subject''s '...
    'to match the template (ie maps from template to individual). '...
    'A forward transform will warp the template image to the individual.']};
%------------------------------------------------------------------------
K = mnu('Time Steps','K',...
        {'1','2','4','8','16','32','64','128','256','512'},...
        {0,1,2,3,4,5,6,7,8,9});
K.val  = {6};
K.help = {...
    ['The number of time points used for solving the '...
     'partial differential equations.  A single time point would be '...
     'equivalent to a small deformation model. '...
     'Smaller values allow faster computations, '...
     'but are less accurate in terms '...
     'of inverse consistency and may result in the one-to-one mapping '...
     'breaking down.']};
% ---------------------------------------------------------------------
template        = cfg_files;
template.tag    = 'template';
template.name   = 'Dartel Template';
template.filter = 'nifti';
template.num    = [0 1];
template.val    = {{''}};
template.help   = {...
['Select the final Template file generated by Dartel. This will be affine '...
 'registered with a TPM file, such that the resulting spatially normalised '...
 'images are closer aligned to MNI space. Leave empty if you do not wish to '...
 'incorporate a transform to MNI space '...
 '(ie just click ``done'' on the file selector, without selecting any images).']};
%------------------------------------------------------------------------
drtl = branch('Dartel flow','dartel',{ffield,forbak,K,template});
drtl.help = {'Imported Dartel flow field.'};
%------------------------------------------------------------------------
other = {drtl,def,id,idbbvox,sn2def};

img          = files('Image to base inverse on','space','nifti',[1 1]);
img.help     = himg;

comp0        = repeat('Composition','comp',other);
comp0.help   = hcomp;

iv0          = branch('Inverse','inv',{comp0,img});
iv0.help     = hinv;

comp1        = repeat('Composition','comp',[other,{iv0},{comp0}]);
comp1.num    = [1 Inf];
comp1.help   = hcomp;

iv1          = branch('Inverse','inv',{comp1,img});
iv1.help     = hinv;

comp2        = repeat('Composition','comp',[other,{iv1},{comp1}]);
comp2.num    = [1 Inf];
comp2.help   = hcomp;

iv2          = branch('Inverse','inv',{comp2,img});
iv2.help     = hinv;

comp         = repeat('Composition','comp',[other,{iv2},{comp2}]);
comp.num     = [1 Inf];
comp.help    = hcomp;

saveas       = entry('Save as','ofname','s',[0 Inf]);
saveas.help  = himgw;

savedas       = entry('Save as','ofname','s',[0 Inf]);
savedas.help  = hdetw;

applyto      = files('Apply to','fnames','nifti',[0 Inf]);
applyto.help = happly;

savepwd      = cfg_const;
savepwd.name = 'Current directory';
savepwd.tag  = 'savepwd';
savepwd.val  = {1};
savepwd.help = {['All created files (deformation fields and warped images) ' ...
                 'are written to the current directory.']};

savesrc      = cfg_const;
savesrc.name = 'Source directories';
savesrc.tag  = 'savesrc';
savesrc.val  = {1};
savesrc.help = {['The combined deformation field is written into the ' ...
                 'directory of the first deformation field, warped images ' ...
                 'are written to the same directories as the source ' ...
                 'images.']};

savedef      = cfg_const;
savedef.name = 'Source directory (deformation)';
savedef.tag  = 'savedef';
savedef.val  = {1};
savedef.help = {['The combined deformation field and the warped images ' ...
                 'are written into the directory of the first deformation ' ...
                 'field.']};

saveusr      = files('Output directory','saveusr','dir',[1 1]);
saveusr.help = {['The combined deformation field and the warped images ' ...
                 'are written into the specified directory.']};

savedir      = cfg_choice;
savedir.name = 'Output destination';
savedir.tag  = 'savedir';
savedir.values = {savepwd savesrc saveusr};
savedir.val  = {savepwd};

savedir1      = cfg_choice;
savedir1.name = 'Output destination';
savedir1.tag  = 'savedir';
savedir1.values = {savepwd saveusr};
savedir1.val  = {savepwd};

interp      = cfg_menu;
interp.name = 'Interpolation';
interp.tag  = 'interp';
interp.labels = {'Nearest neighbour','Trilinear','2nd Degree B-spline',...
'3rd Degree B-Spline ','4th Degree B-Spline ','5th Degree B-Spline',...
'6th Degree B-Spline','7th Degree B-Spline','Categorical'};
interp.values = {0,1,2,3,4,5,6,7,-1};
interp.def  = @(val)spm_get_defaults('normalise.write.interp',val{:});
interp.help    = {
                  ['The method by which the images are sampled when ' ...
                  'being written in a different space. ' ...
                  '(Note that Inf or NaN values are treated as zero, ' ...
                  'rather than as missing data)']
                  '    Nearest Neighbour:'
                  '      - Fastest, but not normally recommended.'
                  '    Trilinear Interpolation:'
                  '      - OK for PET, realigned fMRI, or segmentations'
                  '    B-spline Interpolation:'
                  ['      - Better quality (but slower) interpolation' ...
                  '/* \cite{thevenaz00a}*/, especially with higher ' ...
                  'degree splines. Can produce values outside the ' ...
                  'original range (e.g. small negative values from an ' ...
                  'originally all positive image).']
                  '    Categorical:'
                  ['       - Slow (particularly when there are lots of '...
                  'categories). This is intended to warp categorical images ' ...
                  'such as label maps.']
}';

% ---------------------------------------------------------------------
fwhm         = cfg_entry;
fwhm.tag     = 'fwhm';
fwhm.name    = 'Gaussian FWHM';
fwhm.val     = {[0 0 0]};
fwhm.strtype = 'r';
fwhm.num     = [1 3];
fwhm.help    = {'Specify the full-width at half maximum (FWHM) of the Gaussian blurring kernel in mm. Three values should be entered, denoting the FWHM in the x, y and z directions.'};
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.val     = {''};
prefix.strtype = 's';
prefix.num     = [0 Inf];
prefix.help    = {'The name of the output file(s) will be the name of the input file(s) prefixed with this prefix. Leave empty to use SPM default prefixes.'};
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
mask         = cfg_menu;
mask.tag     = 'mask';
mask.name    = 'Masking';
mask.help    = {'Because of subject motion, different images are likely to have different patterns of zeros from where it was not possible to sample data. With masking enabled, the program searches through the whole time series looking for voxels which need to be sampled from outside the original images. Where this occurs, that voxel is set to zero for the whole set of images (unless the image format can represent NaN, in which case NaNs are used where possible).'};
mask.labels = {
               'Mask images'
               'Dont mask images'
}';
mask.values = {1 0};
mask.val    = {1};
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
preserve         = cfg_menu;
preserve.tag     = 'preserve';
preserve.name    = 'Preserve';
preserve.help    = {
'Preserve Concentrations: Smoothed spatially normalised images (sw*) represent weighted averages of the signal under the smoothing kernel, approximately preserving the intensities of the original images. This option is currently suggested for eg fMRI.'
''
'Preserve Amount: Smoothed and spatially normalised images preserve the total amount of signal from each region in the images (smw*). Areas that are expanded during warping are correspondingly reduced in intensity. This option is suggested for VBM.'
''
'Preserve Labels: This is intended for warping label images. While it is quite slow to run, it is intended to give more accurately warped categorical data.'
}';
preserve.labels = {
                   'Preserve Concentrations (no "modulation")'
                   'Preserve Amount ("modulation")'
                   'Preserve Labels (categorical data)'
}';
preserve.values = {0 1 2};
preserve.val    = {0};
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
fromimage       = cfg_files;
fromimage.name   = 'Image Defined';
fromimage.tag    = 'file';
fromimage.filter = 'nifti';
fromimage.num    = [1 1];
fromimage.help   = {'Use the dimensions, orientation etc of some pre-existing image.'};
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
surfa        = cfg_files;
surfa.name   = 'Surface';
surfa.tag    = 'surface';
surfa.filter = 'mesh';
surfa.num    = [1 Inf];
surfa.help   = {'Select a GIFTI file to warp.'};
% ---------------------------------------------------------------------

bbvox         = cfg_branch;
bbvox.name    = 'User Defined';
bbvox.tag     = 'bbvox';
bbvox.val     = {bb,vox};
bbvox.help   = {[...
'The part of the deformation to use is specified by defining the bounding box and ',...
'voxel sizes that you would like to have. This is probably stating the obvious to many ',...
'but smaller voxels and a broader bounding box will take up more disk space, but may ',...
'give a little more accuracy.']};

deffov        = cfg_choice;
deffov.name   = 'Field of View';
deffov.tag    = 'fov';
deffov.values = {fromimage,bbvox};
deffov.help   = {[...
'The dimensions and voxel size of the resulting deformation may be defined from some image, ',...
'or by specifying voxel sizes and a bounding box.']};

savedef       = cfg_branch;
savedef.name  = 'Save Deformation';
savedef.tag   = 'savedef';
savedef.val   ={saveas,savedir1};
savedef.help  = {'The deformation may be saved to disk as a ``y_*.nii'''' file.'};

savedet       = cfg_branch;
savedet.name  = 'Save Jacobian Determinants';
savedet.tag   = 'savejac';
savedet.val   ={savedas,savedir1};
savedet.help  = {'The Jacobian determinants may be saved to disk as a ``j_*.nii'''' file.'};

pullback      = cfg_branch;
pullback.name = 'Pullback';
pullback.tag  = 'pull';
pullback.val  = {applyto,savedir,interp,mask,fwhm,prefix};
pullback.help = {[...
'This is the old way of warping images, which involves resampling images based on a mapping from ',...
'the new (warped) image space back to the original image.  ',...
'The deformation should be the inverse of the deformation that would be used for the pushforward procedure.']};

weight        = cfg_files;
weight.name   = 'Weight Image';
weight.tag    = 'weight';
weight.filter = 'nifti';
weight.num    = [0 1];
weight.help   = {'Select an image file to weight the warped data with.  This is optional, but the idea is the same as was used by JE Lee et al (2009) in their ``A study of diffusion tensor imaging by tissue-specific, smoothing-compensated voxel-based analysis'''' paper.  In principle, a mask of (eg) white matter could be supplied, such that the warped images contain average signal intensities in WM.'};
weight.val    = {{''}};

% add note on aliasing to fwhm.help for Pushforward
fwhm.help     = {[fwhm.help{1} ' Note that you can specify [0 0 0], ',...
    'but any "modulated" data will show aliasing, which occurs because of ',...
    'the way the warped images are generated.']};

pushfo        = cfg_branch;
pushfo.name   = 'Pushforward';
pushfo.tag    = 'push';
pushfo.val    = {applyto,weight,savedir,deffov,preserve,fwhm,prefix};
pushfo.help   = {[...
'This is a newer way of warping images (for SPM at least), and involves the ',...
'forward pushing of voxel values from the original image into the appropriate place in the warped image. ',...
'The deformation field should be the inverse of the one used for the pullback procedure.'],...
'',...
[...
'``Smoothed'''' (blurred) spatially normalised images are generated in such a ',...
'way that the original signal is preserved. Normalised images are ',...
'generated by a ``pushing'''' rather than a ``pulling'''' (the usual) procedure. ',...
'Note that a procedure related to trilinear interpolation is used, and no masking is done.  It ',...
'is therefore recommended that the images are realigned and resliced ',...
'before they are spatially normalised, in order to benefit from motion correction using higher order interpolation.  Alternatively, contrast images ',...
'generated from unsmoothed native-space fMRI/PET data can be spatially ',...
'normalised for a 2nd level analysis.'],[...
'Two ``preserve'''' options are provided.  One of them should do the ',...
'equavalent of generating smoothed ``modulated'''' spatially normalised ',...
'images.  The other does the equivalent of smoothing the modulated ',...
'normalised fMRI/PET, and dividing by the smoothed Jacobian determinants.']};

pushsurf      = cfg_branch;
pushsurf.name = 'Surface';
pushsurf.tag  = 'surf';
pushsurf.val  = {surfa,savedir};
pushsurf.help = {[...
'Surfaces may be warped using the resulting deformation. ',...
'Note that a procedure similar to the pushforward is used, so the deformation should ',...
'be the inverse of the one that would be used for spatially normalising images via the pullback procedure.']}; 

output        = cfg_repeat;
output.name   = 'Output';
output.tag    = 'out';
output.values = {savedef,pullback, pushfo,pushsurf,savedet};
output.help = {[...
'Various output options are available.  ',...
'The deformation may be saved to disk as a ``y_*.nii'''' file.',...
'Images may be warped using the resulting deformation, either using a ``pullback'''' procedure, or a ``pushforward''''.',...
'The old style of spatial normalisation involved the pullback, whereas the pushforward requires ',...
'the inverse of the deformation used by the pullback.  ',...
'Finally, the deformation may be used to warp a GIFTI surface file.']};

conf         = exbranch('Deformations','defs',{comp,output});
conf.prog    = @spm_deformations;
conf.vout    = @vout;
conf.help    = hsummary;


%==========================================================================
function vo = vout(job)
vo = [];
savedef   = false;
saveimage = false;
savesurf  = false;
savejac   = false;
for i=1:numel(job.out)
    out = job.out{i};
    if isfield(out,'savedef') && ~savedef
        savedef = true;
        if isempty(vo), vo = cfg_dep; else vo(end+1) = cfg_dep; end
        vo(end).sname      = 'Deformation';
        vo(end).src_output = substruct('.','def');
        vo(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
    end
    if (isfield(out,'pull') || isfield(out,'push')) && ~saveimage
        saveimage = true;
        if isempty(vo), vo = cfg_dep; else vo(end+1) = cfg_dep; end
        vo(end).sname      = 'Warped Images';
        vo(end).src_output = substruct('.','warped');
        vo(end).tgt_spec   = cfg_findspec({{'filter','image'}});
    end
    if isfield(out,'surf') && ~savesurf
        savesurf = true;
        if isempty(vo), vo = cfg_dep; else vo(end+1) = cfg_dep; end
        vo(end).sname      = 'Warped Surfaces';
        vo(end).src_output = substruct('.','surf');
        vo(end).tgt_spec   = cfg_findspec({{'filter','mesh'}});
    end
    if isfield(out,'savejac') && ~savejac
        savejac = true;
        if isempty(vo), vo = cfg_dep; else vo(end+1) = cfg_dep; end
        vo(end).sname      = 'Jacobian';
        vo(end).src_output = substruct('.','jac');
        vo(end).tgt_spec   = cfg_findspec({{'filter','image'}});
    end
end


%==========================================================================
function entry_item = entry(name, tag, strtype, num)
entry_item         = cfg_entry;
entry_item.name    = name;
entry_item.tag     = tag;
entry_item.strtype = strtype;
entry_item.num     = num;

function files_item = files(name, tag, fltr, num)
files_item        = cfg_files;
files_item.name   = name;
files_item.tag    = tag;
files_item.filter = fltr;
files_item.num    = num;

function branch_item = branch(name, tag, val)
branch_item      = cfg_branch;
branch_item.name = name;
branch_item.tag  = tag;
branch_item.val  = val;

function exbranch_item = exbranch(name, tag, val)
exbranch_item      = cfg_exbranch;
exbranch_item.name = name;
exbranch_item.tag  = tag;
exbranch_item.val  = val;

function repeat_item = repeat(name, tag, values)
repeat_item        = cfg_repeat;
repeat_item.name   = name;
repeat_item.tag    = tag;
repeat_item.values = values;

function menu_item = mnu(name, tag, labels, values)
menu_item        = cfg_menu;
menu_item.name   = name;
menu_item.tag    = tag;
menu_item.labels = labels;
menu_item.values = values;

