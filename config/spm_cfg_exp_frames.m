function exp_frames = spm_cfg_exp_frames
% 'Expand image frames' - MATLABBATCH configuration
% This MATLABBATCH configuration file has been generated automatically
% by MATLABBATCH using ConfGUI. It describes menu structure, validity
% constraints and links to run time code.
% Changes to this file will be overwritten if the ConfGUI batch is executed again.
% Created at 2009-03-27 21:35:11.
% ---------------------------------------------------------------------
% files NIfTI file(s)
% ---------------------------------------------------------------------
files         = cfg_files;
files.tag     = 'files';
files.name    = 'NIfTI file(s)';
files.help    = {'Files to read. If the same multi-frame image is specified more than once, it will be expanded as often as it is listed.'};
files.filter = 'image';
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
exp_frames         = cfg_exbranch;
exp_frames.tag     = 'exp_frames';
exp_frames.name    = 'Expand image frames';
exp_frames.val     = {files frames };
exp_frames.help    = {'Returns a list of image filenames with appended frame numbers.'};
exp_frames.prog = @(job)spm_run_exp_frames('run',job);
exp_frames.vout = @(job)spm_run_exp_frames('vout',job);
