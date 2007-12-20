function opts = spm_config_reorient
% Configuration file for reorient images
% Apply a given transformation matrix or reorientation parameters by
% left-multiplying the original image orientation with it.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Volkmar Glauche
% $Id: spm_config_movefile.m 549 2006-06-07 12:37:29Z volkmar $

%_______________________________________________________________________


srcfiles.type = 'files';
srcfiles.name = 'Images to reorient';
srcfiles.tag  = 'srcfiles';
srcfiles.filter = 'image';
srcfiles.num  = [0 Inf];
srcfiles.help = {'Select images to reorient.'};

transM.type = 'entry';
transM.name = 'Reorientation matrix';
transM.tag  = 'transM';
transM.strtype = 'e';
transM.num  = [4 4];
p1 = 'Enter a valid 4x4 matrix for reorientation.';
p2 = 'Example: This will L-R flip the images.';
p3 = '   -1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1';
transM.help = {p1,'',p2,'',p3};

transprm.type = 'entry';
transprm.name = 'Reorientation parameters';
transprm.tag  = 'transprm';
transprm.strtype = 'e';
transprm.num  = [1 12];
p0  = 'Enter 12 reorientation parameters.';
p1  = 'P(1)  - x translation';
p2  = 'P(2)  - y translation';
p3  = 'P(3)  - z translation';
p4  = 'P(4)  - x rotation about - {pitch} (radians)';
p5  = 'P(5)  - y rotation about - {roll}  (radians)';
p6  = 'P(6)  - z rotation about - {yaw}   (radians)';
p7  = 'P(7)  - x scaling';
p8  = 'P(8)  - y scaling';
p9  = 'P(9)  - z scaling';
p10 = 'P(10) - x affine';
p11 = 'P(11) - y affine';
p12 = 'P(12) - z affine';
p13 = 'Parameters are entered as listed above and then processed by spm_matrix.';
p14 = ['Example: This will L-R flip the images (extra spaces are inserted between ',...
       'each group for illustration purposes).'];
p15 = '   0 0 0   0 0 0   -1 0 0   0 0 0';

transprm.help = {p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,'',p14,'',p15,''};

transform.type = 'choice';
transform.name = 'Reorient by';
transform.tag  = 'transform';
transform.values = {transM, transprm};
transform.help = {'Specify reorientation method.'};

opts.type = 'branch';
opts.name = 'Reorient images';
opts.tag  = 'reorient';
opts.val  = {srcfiles,transform};
opts.prog = @my_reorient;
opts.vfiles = @vfiles_reorient;
opts.help = {[...
    'This facility allows to reorient images in a batch. The reorientation parameters ' ...
    'can be given either as a 4x4 matrix or as parameters as defined for spm_matrix.m. ' ...
    'The new image orientation will be computed by PRE-multiplying the original '...
    'orientation matrix with the supplied matrix.']};

%------------------------------------------------------------------------

%------------------------------------------------------------------------
function my_reorient(varargin)
job = varargin{1};
if isfield(job.transform,'transprm')
    job.transform.transM = spm_matrix(job.transform.transprm);
end;
spm_progress_bar('Init', numel(job.srcfiles), 'Reorient', 'Images completed');
for k = 1:numel(job.srcfiles)
    M = spm_get_space(job.srcfiles{k});
    spm_get_space(job.srcfiles{k},job.transform.transM*M);
    spm_progress_bar('Set',k);
end;
spm_progress_bar('Clear');
return;
%-------------------------------------------------------------------------

function vf = vfiles_reorient(job)
srcfiles = job.srcfiles{1};
vf    = {spm_select('CPath','image',srcfiles)};
