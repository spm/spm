function out = spm_run_fmri_data(job)
% Set up the design matrix and run a design.
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_fmri_data.m 3756 2010-03-05 18:43:37Z guillaume $


spm('defaults','FMRI');

original_dir = pwd;
[p n e v] = spm_fileparts(job.spmmat{1});
my_cd(p);
load(job.spmmat{1});

% Image filenames
%-------------------------------------------------------------------------
SPM.xY.P = strvcat(job.scans);

% Let SPM configure the design
%-------------------------------------------------------------------------
SPM = spm_fmri_spm_ui(SPM);

if ~isempty(job.mask)&&~isempty(job.mask{1})
    SPM.xM.VM         = spm_vol(job.mask{:});
    SPM.xM.xs.Masking = [SPM.xM.xs.Masking, '+explicit mask'];
end

%-Save SPM.mat
%-------------------------------------------------------------------------
fprintf('%-40s: ','Saving SPM configuration')                          %-#
if spm_matlab_version_chk('7') >= 0
    save('SPM','-V6','SPM');
else
    save('SPM','SPM');
end;

fprintf('%30s\n','...SPM.mat saved')                                   %-#

out.spmmat{1} = fullfile(pwd, 'SPM.mat');
my_cd(original_dir); % Change back dir
fprintf('Done\n')
return
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function my_cd(varargin)
% jobDir must be the actual directory to change to, NOT the job structure.
jobDir = varargin{1};
if ~isempty(jobDir)
    try
        cd(char(jobDir));
        fprintf('Changing directory to: %s\n',char(jobDir));
    catch
        error('Failed to change directory. Aborting run.')
    end
end
return;
