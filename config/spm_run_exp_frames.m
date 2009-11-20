function out = spm_run_exp_frames(cmd, job)
% SPM job execution function for Expand image frames
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_exp_frames.m 3589 2009-11-20 17:17:41Z guillaume $

switch lower(cmd)
    case 'run'
        out.files = {};
        for k = 1:numel(job.files)
            % Strip existing frame spec, if any
            [p n e v] = spm_fileparts(job.files{k});
            V = spm_vol(fullfile(p, [n e]));
            if all(isfinite(job.frames))
                Vframes = job.frames(job.frames <= numel(V));
            else
                Vframes = 1:numel(V);
            end
            Vfiles = cell(numel(Vframes),1);
            for l = 1:numel(Vfiles)
                Vfiles{l} = fullfile(p, [n e ',' num2str(Vframes(l))]);
            end
            out.files = [out.files; Vfiles];
        end
    case 'vout'
        out = cfg_dep;
        out.sname = 'Expanded filename list.';
        out.src_output = substruct('.','files');
        out.tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
end