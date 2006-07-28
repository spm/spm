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
transM.help = {'Enter a valid 4x4 matrix for reorientation.'};

transprm.type = 'entry';
transprm.name = 'Reorientation parameters';
transprm.tag  = 'transprm';
transprm.strtype = 'e';
transprm.num  = [1 12];
transprm.help = {'Enter 12 reorientation parameters (see spm_matrix for details).'};

transform.type = 'choice';
transform.name = 'Reorient by';
transform.tag  = 'transform';
transform.values = {transM, transprm};
transform.help = {'Specify reorientation parameters.'};

opts.type = 'branch';
opts.name = 'Reorient images';
opts.tag  = 'reorient';
opts.val  = {srcfiles,transform};
opts.prog = @my_reorient;
opts.vfiles = @vfiles_reorient;
opts.help = {[...
    'This facilty allows to reorient images in a batch. The reorientation parameters ' ...
    'can be given either as a 4x4 matrix or as parameters as defined for spm_matrix.m. ' ...
    'The new image orientation will be computed by left-multiplying the original '...
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
