function voi = spm_cfg_voi
% SPM Configuration file for VOIs
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_cfg_voi.m 6925 2016-11-09 17:23:40Z guillaume $

% -------------------------------------------------------------------------
% spmmat Select SPM.mat
% -------------------------------------------------------------------------
spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
spmmat.name    = 'Select SPM.mat';
spmmat.help    = {'Select SPM.mat file. If empty, use the SPM.mat selected above.'};
spmmat.filter  = 'mat';
spmmat.ufilter = '^SPM\.mat$';
spmmat.num     = [0 1];
spmmat.val     = {{''}};

% -------------------------------------------------------------------------
% contrast Contrast
% -------------------------------------------------------------------------
contrast         = cfg_entry;
contrast.tag     = 'contrast';
contrast.name    = 'Contrast';
contrast.help    = {'Index of contrast. If more than one index is entered, a conjunction analysis is performed.'};
contrast.strtype = 'n';
contrast.num     = [1 Inf];

% -------------------------------------------------------------------------
% conjunction Conjunction Number
% -------------------------------------------------------------------------
conjunction         = cfg_entry;
conjunction.tag     = 'conjunction';
conjunction.name    = 'Conjunction number';
conjunction.help    = {'Conjunction number. Unused if a simple contrast is entered.'
    'For Conjunction Null, enter 1.'
    'For Global Null, enter the number of selected contrasts.'
    'For Intermediate, enter the number of selected contrasts minus the number of effects under the Null.'}';
conjunction.strtype = 'n';
conjunction.num     = [1 1];
conjunction.val     = {1};

% -------------------------------------------------------------------------
% threshdesc Threshold type
% -------------------------------------------------------------------------
threshdesc         = cfg_menu;
threshdesc.tag     = 'threshdesc';
threshdesc.name    = 'Threshold type';
threshdesc.help    = {''};
threshdesc.labels  = {'FWE' 'none'};
threshdesc.values  = {'FWE' 'none'};
threshdesc.val     = {'none'};

% -------------------------------------------------------------------------
% thresh Threshold
% -------------------------------------------------------------------------
thresh         = cfg_entry;
thresh.tag     = 'thresh';
thresh.name    = 'Threshold';
thresh.help    = {''};
thresh.strtype = 'r';
thresh.num     = [1 1];
thresh.val     = {0.001};

% -------------------------------------------------------------------------
% extent Extent (voxels)
% -------------------------------------------------------------------------
extent         = cfg_entry;
extent.tag     = 'extent';
extent.name    = 'Extent (voxels)';
extent.help    = {''};
extent.strtype = 'w';
extent.num     = [1 1];
extent.val     = {0};

% -------------------------------------------------------------------------
% contrast Contrast
% -------------------------------------------------------------------------
contrastm         = cfg_entry;
contrastm.tag     = 'contrast';
contrastm.name    = 'Contrast';
contrastm.help    = {'Indices of contrast(s).'};
contrastm.strtype = 'n';
contrastm.num     = [1 Inf];

% -------------------------------------------------------------------------
% threshm Threshold
% -------------------------------------------------------------------------
threshm         = cfg_entry;
threshm.tag     = 'thresh';
threshm.name    = 'Uncorrected mask p-value';
threshm.help    = {''};
threshm.strtype = 'r';
threshm.num     = [1 1];
threshm.val     = {0.05};

% -------------------------------------------------------------------------
% mtype Nature of mask
% -------------------------------------------------------------------------
mtype         = cfg_menu;
mtype.tag     = 'mtype';
mtype.name    = 'Nature of mask';
mtype.help    = {''};
mtype.labels  = {'Inclusive' 'Exclusive'};
mtype.values  = {0 1};

% -------------------------------------------------------------------------
% mask Mask definition
% -------------------------------------------------------------------------
mask         = cfg_branch;
mask.tag     = 'mask';
mask.name    = 'Mask definition';
mask.val     = {contrastm threshm mtype};
mask.help    = {''};

% -------------------------------------------------------------------------
% generic Masking
% -------------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Masking';
generic.help    = {''};
generic.values  = {mask};
generic.num     = [0 1];

% -------------------------------------------------------------------------
% map Thresholded SPM
% -------------------------------------------------------------------------
map         = cfg_branch;
map.tag     = 'spm';
map.name    = 'Thresholded SPM';
map.val     = {spmmat contrast conjunction threshdesc thresh extent generic};
map.help    = {'Thresholded SPM'}';

% -------------------------------------------------------------------------
% fix Movement of centre: fixed
% -------------------------------------------------------------------------
fix      = cfg_const;
fix.tag  = 'fixed';
fix.name = 'Fixed';
fix.val  = { 1 };
fix.help = {'Fixed centre'};

% -------------------------------------------------------------------------
% mp SPM index
% -------------------------------------------------------------------------
mp         = cfg_entry;
mp.tag     = 'spm';
mp.name    = 'SPM index';
mp.help    = {'SPM index'};
mp.strtype = 'n';
mp.num     = [1 1];

% -------------------------------------------------------------------------
% mskexp Expression
% -------------------------------------------------------------------------
mskexp         = cfg_entry;
mskexp.tag     = 'mask';
mskexp.name    = 'Mask expression';
mskexp.help    = {'Example expressions (f):'};
mskexp.strtype = 's';
mskexp.num     = [0 Inf];
mskexp.val     = {''};

% -------------------------------------------------------------------------
% glob Movement of centre: Global maximum
% -------------------------------------------------------------------------
glob      = cfg_branch;
glob.tag  = 'global';
glob.name = 'Global maximum';
glob.val  = { mp mskexp };
glob.help = {'Global maximum'};

% -------------------------------------------------------------------------
% loc Movement of centre: Nearest local maximum
% -------------------------------------------------------------------------
loc      = cfg_branch;
loc.tag  = 'local';
loc.name = 'Nearest local maximum';
loc.val  = { mp mskexp };
loc.help = {'Nearest local maximum'};

% -------------------------------------------------------------------------
% supra Movement of centre: Nearest suprathreshold voxel
% -------------------------------------------------------------------------
supra      = cfg_branch;
supra.tag  = 'supra';
supra.name = 'Nearest suprathreshold voxel';
supra.val  = { mp mskexp };
supra.help = {'Nearest suprathreshold voxel.'};

% -------------------------------------------------------------------------
% mvt Movement of centre
% -------------------------------------------------------------------------
mvt         = cfg_choice;
mvt.tag     = 'move';
mvt.name    = 'Movement of centre';
mvt.help    = {'Movement of centre.'};
mvt.values  = {fix glob loc supra};
mvt.val     = { fix };

% -------------------------------------------------------------------------
% centre Centre of Sphere/Box
% -------------------------------------------------------------------------
centre         = cfg_entry;
centre.tag     = 'centre';
centre.name    = 'Centre';
centre.help    = {'Centre [x y z] {mm}.'};
centre.strtype = 'r';
centre.num     = [1 3];

% -------------------------------------------------------------------------
% radius Radius of Sphere
% -------------------------------------------------------------------------
radius         = cfg_entry;
radius.tag     = 'radius';
radius.name    = 'Radius';
radius.help    = {'Sphere radius (mm).'};
radius.strtype = 'r';
radius.num     = [1 1];

% -------------------------------------------------------------------------
% sphere Sphere
% -------------------------------------------------------------------------
sphere         = cfg_branch;
sphere.tag     = 'sphere';
sphere.name    = 'Sphere';
sphere.val     = {centre radius mvt};
sphere.help    = {'Sphere.'}';

% -------------------------------------------------------------------------
% dim Box Dimension
% -------------------------------------------------------------------------
dim         = cfg_entry;
dim.tag     = 'dim';
dim.name    = 'Dimensions';
dim.help    = {'Box dimensions [x y z] {mm}.'};
dim.strtype = 'r';
dim.num     = [1 3];

% -------------------------------------------------------------------------
% box Box
% -------------------------------------------------------------------------
box         = cfg_branch;
box.tag     = 'box';
box.name    = 'Box';
box.val     = {centre dim mvt};
box.help    = {'Box.'}';

% -------------------------------------------------------------------------
% image Image
% -------------------------------------------------------------------------
image         = cfg_files;
image.tag     = 'image';
image.name    = 'Image file';
image.help    = {'Select image.'};
image.filter  = 'image';
image.ufilter = '.*';
image.num     = [1 1];

% -------------------------------------------------------------------------
% threshold Threshold
% -------------------------------------------------------------------------
threshold         = cfg_entry;
threshold.tag     = 'threshold';
threshold.name    = 'Threshold';
threshold.help    = {'Threshold.'};
threshold.strtype = 'r';
threshold.num     = [1 1];
threshold.val     = {0.5};

% -------------------------------------------------------------------------
% mask Mask
% -------------------------------------------------------------------------
mask         = cfg_branch;
mask.tag     = 'mask';
mask.name    = 'Mask Image';
mask.val     = {image threshold};
mask.help    = {'Mask Image.'}';

% -------------------------------------------------------------------------
% list List of labels
% -------------------------------------------------------------------------
list         = cfg_entry;
list.tag     = 'list';
list.name    = 'List of labels';
list.help    = {'List of labels.'};
list.strtype = 'r';
list.num     = [1 Inf];

% -------------------------------------------------------------------------
% label Label Image
% -------------------------------------------------------------------------
label         = cfg_branch;
label.tag     = 'label';
label.name    = 'Label Image';
label.val     = {image list};
label.help    = {'Label Image.'}';

% -------------------------------------------------------------------------
% roi ROI
% -------------------------------------------------------------------------
roi         = cfg_repeat;
roi.tag     = 'roi';
roi.name    = 'Region(s) of Interest';
roi.help    = {'Region(s) of Interest'};
roi.values  = {map sphere box mask label};
roi.num     = [1 Inf];

% -------------------------------------------------------------------------
% expression Expression
% -------------------------------------------------------------------------
expression         = cfg_entry;
expression.tag     = 'expression';
expression.name    = 'Expression';
expression.help    = {['Expression to be evaluated, using i1, i2, ...',...
    'and logical operators (&, |, ~)']};
expression.strtype = 's';
expression.num     = [2 Inf];

% -------------------------------------------------------------------------
% spmmat Select SPM.mat
% -------------------------------------------------------------------------
spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
spmmat.name    = 'Select SPM.mat';
spmmat.help    = {'Select SPM.mat file'};
spmmat.filter  = 'mat';
spmmat.ufilter = '^SPM\.mat$';
spmmat.num     = [1 1];

% -------------------------------------------------------------------------
% adjust Contrast used to adjust data
% -------------------------------------------------------------------------
adjust         = cfg_entry;
adjust.tag     = 'adjust';
adjust.name    = 'Adjust data';
adjust.help    = {'Index of F-contrast used to adjust data. Enter ''0'' for no adjustment. Enter ''NaN'' for adjusting for everything.'}';
adjust.strtype = 'w';
adjust.num     = [1 1];

% -------------------------------------------------------------------------
% session Session index
% -------------------------------------------------------------------------
session         = cfg_entry;
session.tag     = 'session';
session.name    = 'Which session';
session.help    = {'Enter the session number from which you want to extract data.'}';
session.strtype = 'n';
session.num     = [1 1];

% -------------------------------------------------------------------------
% name Name of VOI
% -------------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Name of VOI';
name.help    = {['Name of the VOI mat file that will be saved in the ' ...
'same directory than the SPM.mat file. A ''VOI_'' prefix will be added. ' ...
'' ...
'A binary NIfTI image of the VOI will also be saved.']};
name.strtype = 's';
name.num     = [1 Inf];

% ---------------------------------------------------------------------
% dir Directory
% ---------------------------------------------------------------------
% dir         = cfg_files;
% dir.tag     = 'dir';
% dir.name    = 'Directory';
% dir.help    = {'Output directory.'};
% dir.filter  = 'dir';
% dir.ufilter = '.*';
% dir.num     = [1 1];

% -------------------------------------------------------------------------
% voi VOI
% -------------------------------------------------------------------------
voi         = cfg_exbranch;
voi.tag     = 'voi';
voi.name    = 'Volume of Interest';
voi.val     = {spmmat adjust session name roi expression};
voi.help    = {' VOI time-series extraction of adjusted data (& local eigenimage analysis).'};
voi.prog    = @spm_run_voi;
voi.vout    = @vout_voi;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function dep = vout_voi(varargin)
dep(1)            = cfg_dep;
dep(1).sname      = ' VOI mat File';
dep(1).src_output = substruct('.','voimat');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
dep(2)            = cfg_dep;
dep(2).sname      = ' VOI Image File';
dep(2).src_output = substruct('.','voiimg');
dep(2).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
