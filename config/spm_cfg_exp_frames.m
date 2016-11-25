function exp_frames = spm_cfg_exp_frames
% SPM Configuration file for Expand image frames
%__________________________________________________________________________
% Copyright (C) 2009-2016 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_exp_frames.m 6952 2016-11-25 16:03:13Z guillaume $


%--------------------------------------------------------------------------
% files NIfTI file(s)
%--------------------------------------------------------------------------
files         = cfg_files;
files.tag     = 'files';
files.name    = 'NIfTI file(s)';
files.help    = {'Files to read. If the same multi-frame image is specified more than once, it will be expanded as often as it is listed.'};
files.filter  = 'image';
files.ufilter = '.*';
files.num     = [1 Inf];

%--------------------------------------------------------------------------
% frames Frames
%--------------------------------------------------------------------------
frames         = cfg_entry;
frames.tag     = 'frames';
frames.name    = 'Frames';
frames.help    = {'Frame number(s) requested. Only frames that are actually present in the image file(s) will be listed. Enter ''Inf'' to list all frames.'};
frames.strtype = 'n';
frames.num     = [1 Inf];

%--------------------------------------------------------------------------
% exp_frames Expand image frames
%--------------------------------------------------------------------------
exp_frames      = cfg_exbranch;
exp_frames.tag  = 'exp_frames';
exp_frames.name = 'Expand image frames';
exp_frames.val  = {files frames };
exp_frames.help = {'Return a list of image filenames with appended frame numbers.'};
exp_frames.prog = @run_frames;
exp_frames.vout = @vout_frames;


%==========================================================================
function out = run_frames(job)
out.files = {};
for k = 1:numel(job.files)
    F = spm_select('Expand',spm_file(job.files(k),'number',''));
    if all(isfinite(job.frames))
        F = F(job.frames(job.frames <= numel(F)));
    end
    out.files = [out.files; F];
end


%==========================================================================
function out = vout_frames(job)
out            = cfg_dep;
out.sname      = 'Expanded filename list.';
out.src_output = substruct('.','files');
out.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
