function spmjobs = spm_cfg
% SPM Configuration file for MATLABBATCH
%__________________________________________________________________________
% Copyright (C) 2008-2016 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg.m 6952 2016-11-25 16:03:13Z guillaume $

%--------------------------------------------------------------------------
% Temporal
%--------------------------------------------------------------------------
temporal        = cfg_choice;
temporal.tag    = 'temporal';
temporal.name   = 'Temporal';
temporal.help   = {'Temporal pre-processing functions.'};
temporal.values = { spm_cfg_st };

%--------------------------------------------------------------------------
% Spatial
%--------------------------------------------------------------------------
spatial        = cfg_choice;
spatial.tag    = 'spatial';
spatial.name   = 'Spatial';
spatial.help   = {'Spatial pre-processing functions.'};
spatial.values = { spm_cfg_realign spm_cfg_realignunwarp spm_cfg_coreg spm_cfg_preproc8 spm_cfg_norm spm_cfg_smooth };

%--------------------------------------------------------------------------
% Stats
%--------------------------------------------------------------------------
stats        = cfg_choice;
stats.tag    = 'stats';
stats.name   = 'Stats';
stats.help   = {'Statistical modelling and inference functions.'};
stats.values = { spm_cfg_fmri_spec spm_cfg_fmri_design spm_cfg_fmri_data spm_cfg_factorial_design spm_cfg_model_review spm_cfg_fmri_est spm_cfg_con spm_cfg_results spm_cfg_mfx spm_cfg_bms_map spm_cfg_ppi spm_cfg_setlevel };

%--------------------------------------------------------------------------
% Dynamic Causal Modelling
%--------------------------------------------------------------------------
spm_cfg_dcm_spec        = cfg_choice;
spm_cfg_dcm_spec.tag    = 'spec';
spm_cfg_dcm_spec.name   = 'DCM specification';
spm_cfg_dcm_spec.values = { spm_cfg_dcm_fmri spm_cfg_dcm_meeg };

dcm        = cfg_choice;
dcm.tag    = 'dcm';
dcm.name   = 'DCM';
dcm.help   = {'Dynamic Causal Modelling.'};
dcm.values = { spm_cfg_dcm_spec spm_cfg_dcm_est spm_cfg_dcm_bms spm_cfg_dcm_peb };

%--------------------------------------------------------------------------
% Util
%--------------------------------------------------------------------------
spm_cfg_import        = cfg_choice;
spm_cfg_import.tag    = 'import';
spm_cfg_import.name   = 'Import';
spm_cfg_import.help   = {'Import.'};
spm_cfg_import.values = { spm_cfg_dicom spm_cfg_minc spm_cfg_ecat spm_cfg_parrec };

util        = cfg_choice;
util.tag    = 'util';
util.name   = 'Util';
util.help   = {'Utility tools.'};
util.values = { spm_cfg_disp spm_cfg_checkreg spm_cfg_render spm_cfg_import spm_cfg_imcalc spm_cfg_reorient spm_cfg_voi spm_cfg_cdir spm_cfg_md spm_cfg_bbox spm_cfg_deface spm_cfg_deformations spm_cfg_tissue_volumes spm_cfg_print spm_cfg_cat spm_cfg_split spm_cfg_exp_frames spm_cfg_sendmail };

%--------------------------------------------------------------------------
% Tools
%--------------------------------------------------------------------------
tools        = cfg_choice;
tools.tag    = 'tools';
tools.name   = 'Tools';
tools.help   = {'Other tools.', ...
                ['Toolbox configuration files should be placed in the ' ...
                 'toolbox directory, with their own *_cfg_*.m files. ' ...
                 'If you write a toolbox, then you can include it in ' ...
                 'this directory - but remember to try to keep the ' ...
                 'function names unique (to reduce  clashes with other ' ...
                 'toolboxes).'], ...
                ['See spm_cfg.m or MATLABBATCH documentation ' ...
                 'for information about the form of SPM''s configuration ' ...
                 'files.']};
tools.values = {};
if isdeployed
    %-In compiled mode, cfg_master will take care of toolbox detection
    % See spm_make_standalone.m
    tools.values = spm_cfg_static_tools;
else
    %-Toolbox directories autodetection
    tbxdir = spm_get_defaults('tbx.dir');
    for i=1:numel(tbxdir)
        [unused,d] = spm_select('FPList',tbxdir{i});
        if isempty(d), d = {}; else d = cellstr(d); end
        ft = {}; dt = {};
        %-Look for '*_cfg_*.m' files in these directories
        for j=1:numel(d)
            f = cellstr(spm_select('List',d{j},'^[^(._)].*_cfg_.*\.m$'));
            if ~isempty(f{1})
                ft = [ft f{:}];
                dt = [dt repmat(d(j),1,numel(f))];
            end
        end
        
        %-Path to the toolbox MUST be added to matlabpath in the
        % configuration or 'prog' file of said toolbox with, e.g.:
        % >> if ~isdeployed
        % >>   addpath(fileparts(mfilename('fullpath')));
        % >> end
        cwd = pwd;
        for j=1:numel(ft)
            try
                cd(dt{j});
                tools.values{end+1} = feval(strtok(ft{j},'.'));
            catch
                disp(['Loading of toolbox ' fullfile(dt{j},ft{j}) ' failed.']);
            end
        end
        cd(cwd);
    end
end

%==========================================================================
% spmjobs SPM
%==========================================================================
spmjobs        = cfg_choice;
spmjobs.tag    = 'spm';
spmjobs.name   = 'SPM';
spmjobs.help   = {
    '%* Statistical Parametric Mapping'
    ''
    'Statistical Parametric Mapping refers to the construction and assessment of spatially extended statistical processes used to test hypotheses about functional imaging data. These ideas have been instantiated in software that is called SPM.'
    ''
    'The SPM software package has been designed for the analysis of brain imaging data sequences. The sequences can be a series of images from different cohorts, or time-series from the same subject.'
    ''
    'The current release is designed for the analysis of fMRI, PET, SPECT, EEG and MEG.'
    ''
    }';
spmjobs.values = { temporal spatial stats dcm spm_cfg_eeg util tools };
spmjobs.rewrite_job = @spm_rewrite_job;
