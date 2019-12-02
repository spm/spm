function out = spm_run_fmri_data(job)
% Set up the design matrix and run a design
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2005-2019 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_fmri_data.m 7739 2019-12-02 14:00:18Z guillaume $


original_dir = pwd;
cd(spm_file(job.spmmat{1},'fpath'));

load(fullfile(pwd,'SPM.mat'));

%-Image filenames
%--------------------------------------------------------------------------
SPM.xY.P = char(job.scans);

%-Let SPM configure the design
%--------------------------------------------------------------------------
SPM = spm_fmri_spm_ui(SPM);

if ~isempty(job.mask{1})
    SPM.xM.VM         = spm_data_hdr_read(job.mask{:});
    SPM.xM.xs.Masking = [SPM.xM.xs.Masking, '+explicit mask'];
end

%-Save SPM.mat
%--------------------------------------------------------------------------
fprintf('%-40s: ','Saving SPM configuration')                           %-#
fmt = spm_get_defaults('mat.format');
s = whos('SPM');
if s.bytes > 2147483647, fmt = '-v7.3'; end
save('SPM.mat','SPM', fmt);
fprintf('%30s\n','...SPM.mat saved')                                    %-#

out.spmmat{1} = fullfile(pwd, 'SPM.mat');

cd(original_dir);
