function imcalc = spm_cfg_imcalc
% SPM Configuration file for ImCalc
%__________________________________________________________________________
% Copyright (C) 2008-2015 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_imcalc.m 6536 2015-08-26 13:53:45Z guillaume $

%--------------------------------------------------------------------------
% input Input Images
%--------------------------------------------------------------------------
input         = cfg_files;
input.tag     = 'input';
input.name    = 'Input Images';
input.help    = {'These are the images that are used by the calculator.  They are referred to as i1, i2, i3, etc in the order that they are specified.'};
input.filter  = 'image';
input.ufilter = '.*';
input.num     = [1 Inf];

%--------------------------------------------------------------------------
% output Output Filename
%--------------------------------------------------------------------------
output         = cfg_entry;
output.tag     = 'output';
output.name    = 'Output Filename';
output.help    = {'The output image is written to current working directory unless a valid full pathname is given. If a path name is given here, the output directory setting will be ignored.'
    'If the field is left empty, i.e. set to '''', then the name of the 1st input image, preprended with ''i'', is used (change this letter in the spm_defaults if necessary).'};
output.strtype = 's';
output.num     = [0 Inf];
output.val     = {'output'};

%--------------------------------------------------------------------------
% outdir Output Directory
%--------------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output Directory';
outdir.val{1}  = {''};
outdir.help    = {'Files produced by this function will be written into this output directory. If no directory is given, images will be written to current working directory. If both output filename and output directory contain a directory, then output filename takes precedence.'};
outdir.filter  = 'dir';
outdir.ufilter = '.*';
outdir.num     = [0 1];

%--------------------------------------------------------------------------
% expression Expression
%--------------------------------------------------------------------------
expression         = cfg_entry;
expression.tag     = 'expression';
expression.name    = 'Expression';
expression.help    = {
                      'Example expressions (f):'
                      '    * Mean of six images (select six images)'
                      '       f = ''(i1+i2+i3+i4+i5+i6)/6'''
                      '    * Make a binary mask image at threshold of 100'
                      '       f = ''i1>100'''
                      '    * Make a mask from one image and apply to another'
                      '       f = ''i2.*(i1>100)'''
                      '             - here the first image is used to make the mask, which is applied to the second image'
                      '    * Sum of n images'
                      '       f = ''i1 + i2 + i3 + i4 + i5 + ...'''
                      '    * Sum of n images (when reading data into a data-matrix - use dmtx arg)'
                      '       f = ''sum(X)'''
}';
expression.strtype = 's';
expression.num     = [2 Inf];

%--------------------------------------------------------------------------
% name Name
%--------------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Name';
name.help    = {'Variable name used in expression.'};
name.strtype = 's';
name.num     = [1 Inf];

%--------------------------------------------------------------------------
% value Value
%--------------------------------------------------------------------------
value         = cfg_entry;
value.tag     = 'value';
value.name    = 'Value';
value.help    = {'Value of the variable.'};
value.strtype = 'e';
value.num     = [Inf Inf];

%--------------------------------------------------------------------------
% var Variable
%--------------------------------------------------------------------------
var      = cfg_branch;
var.tag  = 'var';
var.name = 'Variable';
var.val  = { name value };
var.help = {'Additional variable which can be used in expression.'};

%--------------------------------------------------------------------------
% generic Additional Variables
%--------------------------------------------------------------------------
generic        = cfg_repeat;
generic.tag    = 'generic';
generic.name   = 'Additional Variables';
generic.help   = {'Additional variables which can be used in expression.'};
generic.values = {var};
generic.num    = [0 Inf];

%--------------------------------------------------------------------------
% dmtx Data Matrix
%--------------------------------------------------------------------------
dmtx         = cfg_menu;
dmtx.tag     = 'dmtx';
dmtx.name    = 'Data Matrix';
dmtx.help    = {'If the dmtx flag is set, then images are read into a data matrix X (rather than into separate variables i1, i2, i3,...). The data matrix  should be referred to as X, and contains images in rows. Computation is plane by plane, so in data-matrix mode, X is a NxK matrix, where N is the number of input images [prod(size(Vi))], and K is the number of voxels per plane [prod(Vi(1).dim(1:2))].'};
dmtx.labels  = {'No - don''t read images into data matrix'
                'Yes -  read images into data matrix'}';
dmtx.values  = {0 1};
dmtx.val     = {0};

%--------------------------------------------------------------------------
% mask Masking
%--------------------------------------------------------------------------
mask         = cfg_menu;
mask.tag     = 'mask';
mask.name    = 'Masking';
mask.help    = {'For data types without a representation of NaN, implicit zero masking assumes that all zero voxels are to be treated as missing, and treats them as NaN. NaN''s are written as zero (by spm_write_plane), for data types without a representation of NaN.'};
mask.labels  = {'No implicit zero mask'
                'Implicit zero mask'
                'NaNs should be zeroed'}';
mask.values  = {0 1 -1};
mask.val     = {0};

%--------------------------------------------------------------------------
% interp Interpolation
%--------------------------------------------------------------------------
interp         = cfg_menu;
interp.tag     = 'interp';
interp.name    = 'Interpolation';
interp.help    = {
                  'With images of different sizes and orientations, the size and orientation of the first is used for the output image. A warning is given in this situation. Images are sampled into this orientation using the interpolation specified by the hold parameter.'
                  ''
                  'The method by which the images are sampled when being written in a different space.'
                  '    Nearest Neighbour'
                  '    - Fastest, but not normally recommended.'
                  '    Trilinear Interpolation'
                  '    - OK for PET, or realigned fMRI.'
                  '    Sinc Interpolation'
                  '    - Better quality (but slower) interpolation, especially'
                  '      with higher degrees.'
}';
interp.labels  = {
                  'Nearest neighbour'
                  'Trilinear'
                  '2nd Degree Sinc'
                  '3rd Degree Sinc'
                  '4th Degree Sinc'
                  '5th Degree Sinc'
                  '6th Degree Sinc'
                  '7th Degree Sinc'
}';
interp.values  = {0 1 -2 -3 -4 -5 -6 -7};
interp.val     = {1};

%--------------------------------------------------------------------------
% dtype Data Type
%--------------------------------------------------------------------------
dtype         = cfg_menu;
dtype.tag     = 'dtype';
dtype.name    = 'Data Type';
dtype.help    = {'Data-type of output image'};
dtype.labels  = {
                 'UINT8   - unsigned char'
                 'INT16   - signed short'
                 'INT32   - signed int'
                 'FLOAT32 - single prec. float'
                 'FLOAT64 - double prec. float'
}';
dtype.values  = {spm_type('uint8') spm_type('int16') spm_type('int32') ...
                 spm_type('float32') spm_type('float64')};
dtype.val     = {spm_type('int16')};

%--------------------------------------------------------------------------
% options Options
%--------------------------------------------------------------------------
options      = cfg_branch;
options.tag  = 'options';
options.name = 'Options';
options.val  = {dmtx mask interp dtype };
options.help = {'Options for image calculator'};

%--------------------------------------------------------------------------
% imcalc Image Calculator
%--------------------------------------------------------------------------
imcalc      = cfg_exbranch;
imcalc.tag  = 'imcalc';
imcalc.name = 'Image Calculator';
imcalc.val  = {input output outdir expression generic options };
imcalc.help = {'The image calculator is for performing user-specified algebraic manipulations on a set of images, with the result being written out as an image. The user is prompted to supply images to work on, a filename for the output image, and the expression to evaluate. The expression should be a standard MATLAB expression, within which the images should be referred to as i1, i2, i3,... etc.'};
imcalc.prog = @my_spm_imcalc;
imcalc.vout = @vout;


%==========================================================================
% function out = my_spm_imcalc(job)
%==========================================================================
function out = my_spm_imcalc(job)
[p,nam,ext] = spm_fileparts(job.output);
if isempty(p)
    if isempty(job.outdir{1})
        p = pwd;
    else
        p = job.outdir{1};
    end
end
if isempty(nam)
    nam = [spm_get_defaults('imcalc.prefix') spm_file(job.input{1},'basename')];
    ext = ['.' spm_file(job.input{1},'ext')];
end
if isempty(ext)
    ext = spm_file_ext;
end
out.files = { fullfile(p,[nam ext]) };
extra_vars = {};
if numel(job.var)
    extra_vars = { job.var };
end

spm_imcalc(char(job.input), out.files{1}, job.expression, job.options, extra_vars{:});

cmd = 'spm_image(''display'',''%s'')';
fprintf('ImCalc Image: %s\n',spm_file(out.files{1},'link',cmd));


%==========================================================================
% function dep = vout(job)
%==========================================================================
function dep = vout(job)
dep = cfg_dep;
if ~ischar(job.output) || strcmp(job.output, '<UNDEFINED>')
    dep.sname  = 'ImCalc Computed Image';
else
    dep.sname  = sprintf('ImCalc Computed Image: %s', job.output);
end
dep.src_output = substruct('.','files');
dep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
