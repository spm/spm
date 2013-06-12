function spmjobs = spm_cfg
% SPM Configuration file for MATLABBATCH
%__________________________________________________________________________
% Copyright (C) 2008-2011 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg.m 5545 2013-06-12 12:08:37Z gareth $

%--------------------------------------------------------------------------
% Temporal
%--------------------------------------------------------------------------
temporal         = cfg_choice;
temporal.tag     = 'temporal';
temporal.name    = 'Temporal';
temporal.help    = {'Temporal pre-processing functions.'};
temporal.values  = { spm_cfg_st };

%--------------------------------------------------------------------------
% Spatial
%--------------------------------------------------------------------------
spatial         = cfg_choice;
spatial.tag     = 'spatial';
spatial.name    = 'Spatial';
spatial.help    = {'Various spatial and other pre-processing functions.'};
spatial.values  = { spm_cfg_realign spm_cfg_realignunwarp spm_cfg_coreg spm_cfg_preproc8 spm_cfg_norm spm_cfg_smooth };

%--------------------------------------------------------------------------
% Stats
%--------------------------------------------------------------------------
stats         = cfg_choice;
stats.tag     = 'stats';
stats.name    = 'Stats';
stats.help    = {'Various analysis utilities.'};
stats.values  = { spm_cfg_fmri_spec spm_cfg_fmri_design spm_cfg_fmri_data spm_cfg_mfx spm_cfg_factorial_design spm_cfg_fmri_est spm_cfg_con spm_cfg_results spm_cfg_bms spm_cfg_ppi spm_cfg_setlevel };


%--------------------------------------------------------------------------
% Util
%--------------------------------------------------------------------------
util         = cfg_choice;
util.tag     = 'util';
util.name    = 'Util';
util.help    = {'Various useful tools.'};
util.values  = { spm_cfg_disp spm_cfg_checkreg spm_cfg_imcalc spm_cfg_reorient spm_cfg_voi spm_cfg_dicom spm_cfg_minc spm_cfg_ecat spm_cfg_spm_surf spm_cfg_cdir spm_cfg_md spm_cfg_bbox spm_cfg_deformations spm_cfg_tissue_volumes spm_cfg_print spm_cfg_cat spm_cfg_exp_frames spm_cfg_sendmail };

%--------------------------------------------------------------------------
% Tools
%--------------------------------------------------------------------------
tools         = cfg_choice;
tools.tag     = 'tools';
tools.name    = 'Tools';
tools.help    = {'Other tools', ...
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
    tools.values = spm_cfg_static_tools;
else
    %-Toolbox directories autodetection
    try
        tbxdir = spm_get_defaults('tbx.dir');
    catch
        tbxdir = { fullfile(spm('Dir'),'toolbox') };
    end
    for i=1:numel(tbxdir)
        d  = dir(tbxdir{i});
        d  = {d([d.isdir]).name};
        d(strncmp('.',d,1)) = [];
        d  = [{''} d];
        ft = {}; dt = {};
        %-Look for '*_cfg_*.m' files in these directories
        for j=1:length(d)
            d2 = fullfile(tbxdir{i},d{j});
            di = dir(d2); di = {di(~[di.isdir]).name};
            f2 = regexp(di,'.*_cfg_.*\.m$');
            if ~iscell(f2), f2 = {f2}; end
            fi = {di{~cellfun('isempty',f2)}};
            if ~isempty(fi)
                ft = {ft{:} fi{:}};
                dt(end+1:end+length(fi)) = deal({d2});
            end
        end
        if ~isempty(ft)
            % The toolbox developer MUST add path to his/her toolbox in his/her
            % 'prog' function, with a command line like:
            % >> if ~isdeployed,
            % >>   addpath(fileparts(mfilename('fullpath')),'-end');
            % >> end
            cwd = pwd;
            for j=1:length(ft)
                
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
end

%==========================================================================
% spmjobs SPM
%==========================================================================
spmjobs         = cfg_choice;
spmjobs.tag     = 'spm';
spmjobs.name    = 'SPM';
spmjobs.help    = {
                '%* Statistical Parametric Mapping'
                ''
                'Statistical Parametric Mapping refers to the construction and assessment of spatially extended statistical processes used to test hypotheses about functional imaging data. These ideas have been instantiated in software that is called SPM.'
                ''
                'The SPM software package has been designed for the analysis of brain imaging data sequences. The sequences can be a series of images from different cohorts, or time-series from the same subject.'
                ''
                'The current release is designed for the analysis of fMRI, PET, SPECT, EEG and MEG.'
                ''
}';
spmjobs.values  = { temporal spatial stats spm_cfg_eeg util tools};
spmjobs.rewrite_job = @spm_rewrite_job;
