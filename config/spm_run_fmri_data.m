function out = spm_run_fmri_data(job)
% Set up the design matrix and run a design.
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_fmri_data.m 4470 2011-09-08 14:42:38Z guillaume $


original_dir = pwd;
p = spm_file(job.spmmat{1},'fpath');
my_cd(p);
load(job.spmmat{1});

%-Image filenames
%--------------------------------------------------------------------------
SPM.xY.P = strvcat(job.scans);

%-Let SPM configure the design
%--------------------------------------------------------------------------
SPM = spm_fmri_spm_ui(SPM);

if ~isempty(job.mask)&&~isempty(job.mask{1})
    SPM.xM.VM         = spm_vol(job.mask{:});
    SPM.xM.xs.Masking = [SPM.xM.xs.Masking, '+explicit mask'];
end

%-Save SPM.mat
%--------------------------------------------------------------------------
fprintf('%-40s: ','Saving SPM configuration')                           %-#
if spm_check_version('matlab','7') >= 0
    save('SPM.mat','-V6','SPM');
else
    save('SPM.mat','SPM');
end

fprintf('%30s\n','...SPM.mat saved')                                    %-#

out.spmmat{1} = fullfile(pwd, 'SPM.mat');
my_cd(original_dir);
fprintf('Done\n')

%==========================================================================
function my_cd(jobDir)
if ~isempty(jobDir)
    try
        cd(char(jobDir));
    catch
        error('Failed to change directory. Aborting run.')
    end
end
