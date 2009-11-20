function exp_frames = spm_cfg_exp_frames
% SPM Configuration file for Expand image frames
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_exp_frames.m 3589 2009-11-20 17:17:41Z guillaume $

% ---------------------------------------------------------------------
% files NIfTI file(s)
% ---------------------------------------------------------------------
files         = cfg_files;
files.tag     = 'files';
files.name    = 'NIfTI file(s)';
files.help    = {'Files to read. If the same multi-frame image is specified more than once, it will be expanded as often as it is listed.'};
files.filter  = 'image';
files.ufilter = '.*';
files.num     = [1 Inf];
% ---------------------------------------------------------------------
% frames Frames
% ---------------------------------------------------------------------
frames         = cfg_entry;
frames.tag     = 'frames';
frames.name    = 'Frames';
frames.help    = {'Frame number(s) requested. Only frames that are actually present in the image file(s) will be listed. Enter ''Inf'' to list all frames.'};
frames.strtype = 'n';
frames.num     = [1  Inf];
% ---------------------------------------------------------------------
% exp_frames Expand image frames
% ---------------------------------------------------------------------
exp_frames      = cfg_exbranch;
exp_frames.tag  = 'exp_frames';
exp_frames.name = 'Expand image frames';
exp_frames.val  = {files frames };
exp_frames.help = {'Returns a list of image filenames with appended frame numbers.'};
exp_frames.prog = @(job)spm_run_exp_frames('run',job);
exp_frames.vout = @(job)spm_run_exp_frames('vout',job);
