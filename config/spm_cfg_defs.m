function conf = spm_cfg_defs
% Configuration file for deformation jobs.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_cfg_defs.m 4136 2010-12-09 22:22:28Z guillaume $

hsummary = {[...
'This is a utility for working with deformation fields. ',...
'They can be loaded, inverted, combined etc, and the results ',...
'either saved to disk, or applied to some image.']};

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
'In this utility, the left-to-right order of the compositions is ',...
'from top to bottom (note that the rightmost deformation would ',...
'actually be applied first). ',...
'i.e. ...((first o second) o third)...o last. The resulting deformation field will ',...
'have the same domain as the first deformation specified, and will map ',...
'to voxels in the codomain of the last specified deformation field.']};

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
'Deformations can be thought of as vector fields. These can be represented ',...
'by three-volume images.']};

himgw = {[...
'Save the result as a three-volume image.  "y_" will be prepended to the ',...
'filename.  The result will be written to the current directory.']};

happly = {[...
'Apply the resulting deformation field to some images. ',...
'The warped images will be written to the current directory, and the ',...
'filenames prepended by "w".  Note that trilinear interpolation is used ',...
'to resample the data, so the original values in the images will ',...
'not be preserved.']};

hmatname = {...
'Specify the _sn.mat to be used.'};

himg = {...
'Specify the image file on which to base the dimensions, orientation etc.'};

hid = {[...
'This option generates an identity transform, but this can be useful for ',...
'changing the dimensions of the resulting deformation (and any images that ',...
'are generated from it).  Dimensions, orientation etc are derived from ',...
'an image.']};

def          = files('Deformation Field','def','.*y_.*\.nii$',[1 1]);
def.help     = himgr;

matname      = files('Parameter File','matname','.*_sn\.mat$',[1 1]);
matname.help = hmatname;

vox          = entry('Voxel sizes','vox','e',[1 3]);
vox.val      = {[NaN NaN NaN]};
vox.help     = hvox;

bb           = entry('Bounding box','bb','e',[2 3]);
bb.val       = {[NaN NaN NaN;NaN NaN NaN]};
bb.help      = hbb;

sn2def       = branch('Imported _sn.mat','sn2def',{matname,vox,bb});
sn2def.help  = hsn;

img          = files('Image to base Id on','space','image',[1 1]);
img.help     = himg;

id           = branch('Identity (Reference Image)','id',{img});
id.help      = hid;

voxid        = entry('Voxel sizes','vox','e',[1 3]);
bbid         = entry('Bounding box','bb','e',[2 3]);

idbbvox      = branch('Identity (Bounding Box and Voxel Size)','idbbvox',{voxid, bbid});
id.help      = hid;

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
    'The direction of the DARTEL flow.  '...
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
%------------------------------------------------------------------------
drtl = branch('DARTEL flow','dartel',{ffield,forbak,K});
drtl.help = {'Imported DARTEL flow field.'};
%------------------------------------------------------------------------
other = {sn2def,drtl,def,id,idbbvox};

img          = files('Image to base inverse on','space','image',[1 1]);
img.help     = himg;

comp0        = repeat('Composition','comp',other);
comp0.help   = hcomp;

iv0          = branch('Inverse','inv',{comp0,img});
iv0.help     = hinv;

comp1        = repeat('Composition','comp',{other{:},iv0,comp0});
comp1.num    = [1 Inf];
comp1.help   = hcomp;

iv1          = branch('Inverse','inv',{comp1,img});
iv1.help     = hinv;

comp2        = repeat('Composition','comp',{other{:},iv1,comp1});
comp2.num    = [1 Inf];
comp2.help   = hcomp;

iv2          = branch('Inverse','inv',{comp2,img});
iv2.help     = hinv;

comp         = repeat('Composition','comp',{other{:},iv2,comp2});
comp.num     = [1 Inf];
comp.help    = hcomp;

saveas       = entry('Save as','ofname','s',[0 Inf]);
saveas.val   = {''};
saveas.help  = himgw;

applyto      = files('Apply to','fnames','image',[0 Inf]);
applyto.val  = {''};
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
savedir.values = {savepwd savesrc savedef saveusr};
savedir.val  = {savepwd};

interp      = cfg_menu;
interp.name = 'Interpolation';
interp.tag  = 'interp';
interp.labels = {'Nearest neighbour','Trilinear','2nd Degree B-spline',...
'3rd Degree B-Spline ','4th Degree B-Spline ','5th Degree B-Spline',...
'6th Degree B-Spline','7th Degree B-Spline'};
interp.values = {0,1,2,3,4,5,6,7};
interp.def  = @(val)spm_get_defaults('normalise.write.interp',val{:});
interp.help    = {
                  ['The method by which the images are sampled when ' ...
                  'being written in a different space. ' ...
                  '(Note that Inf or NaN values are treated as zero, ' ...
                  'rather than as missing data)']
                  '    Nearest Neighbour:'
                  '      - Fastest, but not normally recommended.'
                  '    Bilinear Interpolation:'
                  '      - OK for PET, realigned fMRI, or segmentations'
                  '    B-spline Interpolation:'
                  ['      - Better quality (but slower) interpolation' ...
                  '/* \cite{thevenaz00a}*/, especially with higher ' ...
                  'degree splines. Can produce values outside the ' ...
                  'original range (e.g. small negative values from an ' ...
                  'originally all positive image).']
}';

conf         = exbranch('Deformations','defs',{comp,saveas,applyto,savedir,interp});
conf.prog    = @spm_defs;
conf.vout    = @vout;
conf.help    = hsummary;
return;
%_______________________________________________________________________

%_______________________________________________________________________

function vo = vout(job)
vo = [];
if ~isempty(job.ofname) && ~isequal(job.ofname,'<UNDEFINED>') 
    vo            = cfg_dep;
    vo.sname      = 'Combined deformation';
    vo.src_output = substruct('.','def');
    vo.tgt_spec   = cfg_findspec({{'filter','image','filter','nifti'}});
end
if ~isempty(job.fnames) && ~isequal(job.fnames, {''})
    if isempty(vo), vo = cfg_dep; else vo(end+1) = cfg_dep; end
    vo(end).sname      = 'Warped images';
    vo(end).src_output = substruct('.','warped');
    vo(end).tgt_spec   = cfg_findspec({{'filter','image'}});
end
return;
%_______________________________________________________________________

%_______________________________________________________________________
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
