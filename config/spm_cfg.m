function spmjobs = spm_cfg
% SPM Configuration file
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg.m 4257 2011-03-18 15:28:29Z vladimir $

%_______________________________________________________________________
% temporal Temporal
% ---------------------------------------------------------------------
temporal         = cfg_choice;
temporal.tag     = 'temporal';
temporal.name    = 'Temporal';
temporal.help    = {'Temporal pre-processing functions.'};
temporal.values  = { spm_cfg_st};
% ---------------------------------------------------------------------
% spatial Spatial
% ---------------------------------------------------------------------
spatial         = cfg_choice;
spatial.tag     = 'spatial';
spatial.name    = 'Spatial';
spatial.help    = {'Various spatial and other pre-processing functions.'};
spatial.values  = { spm_cfg_realign spm_cfg_realignunwarp spm_cfg_coreg spm_cfg_preproc spm_cfg_normalise spm_cfg_smooth};
% ---------------------------------------------------------------------
% stats Stats
% ---------------------------------------------------------------------
stats         = cfg_choice;
stats.tag     = 'stats';
stats.name    = 'Stats';
stats.help    = {'Various analysis utilities.'};
stats.values  = { spm_cfg_fmri_spec spm_cfg_fmri_design spm_cfg_fmri_data spm_cfg_factorial_design spm_cfg_fmri_est spm_cfg_con spm_cfg_results spm_cfg_bms spm_cfg_ppi};
% ---------------------------------------------------------------------
% meeg preprocessing
% ---------------------------------------------------------------------
meegprep        = cfg_choice;
meegprep.tag    = 'preproc';
meegprep.name   = 'M/EEG Preprocessing';
meegprep.help   = {'M/EEG preprocessing.'};
meegprep.values = {spm_cfg_eeg_montage spm_cfg_eeg_filter spm_cfg_eeg_bc spm_cfg_eeg_artefact spm_cfg_eeg_downsample spm_cfg_eeg_merge spm_cfg_eeg_fuse}; 
% ---------------------------------------------------------------------
% meeg time-frequency
% ---------------------------------------------------------------------
meegtf        = cfg_choice;
meegtf.tag    = 'tf';
meegtf.name   = 'M/EEG Time-frequency';
meegtf.help   = {'M/EEG time-frequency.'};
meegtf.values = {spm_cfg_eeg_tf spm_cfg_eeg_tf_rescale}; 
% ---------------------------------------------------------------------
% meeg source reconstruction
% ---------------------------------------------------------------------
source        = cfg_choice;
source.tag    = 'source';
source.name   = 'M/EEG Source reconstruction';
source.help   = {'M/EEG source reconstruction.'};
source.values = {spm_cfg_eeg_inv_headmodel, spm_cfg_eeg_inv_invert, spm_cfg_eeg_inv_results, spm_cfg_eeg_inv_extract}; 
% ---------------------------------------------------------------------
% meeg Meeg
% ---------------------------------------------------------------------
meeg         = cfg_choice;
meeg.tag     = 'meeg';
meeg.name    = 'M/EEG';
meeg.help    = {'M/EEG functions.'};
meeg.values  = { spm_cfg_eeg_convert spm_cfg_eeg_epochs meegprep spm_cfg_eeg_average spm_cfg_eeg_review spm_cfg_eeg_contrast spm_cfg_eeg_grandmean spm_cfg_eeg_convert2images meegtf source};
% ---------------------------------------------------------------------
% util Util
% ---------------------------------------------------------------------
util         = cfg_choice;
util.tag     = 'util';
util.name    = 'Util';
util.help    = {'Various useful tools.'};
util.values  = { spm_cfg_disp spm_cfg_checkreg spm_cfg_imcalc spm_cfg_reorient spm_cfg_voi spm_cfg_dicom spm_cfg_minc spm_cfg_ecat spm_cfg_spm_surf spm_cfg_cdir spm_cfg_md spm_cfg_movefile spm_cfg_deletefiles spm_cfg_defs spm_cfg_print spm_cfg_cat spm_cfg_exp_frames};
% ---------------------------------------------------------------------
% tools Tools
% ---------------------------------------------------------------------
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
if isdeployed
    tools.values = spm_cfg_static_tools;
else
    %-Toolbox autodetection
    % In compiled mode, cfg_master will take care of this
    % Disable warnings when converting SPM5 toolboxes - set this to 'on' to
    % debug problems with SPM5 toolboxes
    warning('off','matlabbatch:cfg_struct2cfg:verb');
    %-Get the list of toolbox directories
    tbxdir = fullfile(spm('Dir'),'toolbox');
    d  = dir(tbxdir); d = {d([d.isdir]).name};
    dd = regexp(d,'^\.');
    %(Beware, regexp returns an array if input cell array is of dim 0 or 1)
    if ~iscell(dd), dd = {dd}; end
    d  = {'' d{cellfun('isempty',dd)}};
    ft = {}; dt = {};
    ftc = {}; dtc = {};
    %-Look for '*_cfg_*.m' or '*_config_*.m' files in these directories
    for i=1:length(d)
        d2 = fullfile(tbxdir,d{i});
        di = dir(d2); di = {di(~[di.isdir]).name};
        f2 = regexp(di,'.*_cfg_.*\.m$');
        if ~iscell(f2), f2 = {f2}; end
        fi = {di{~cellfun('isempty',f2)}};
        if ~isempty(fi)
            ft = {ft{:} fi{:}};
            dt(end+1:end+length(fi)) = deal({d2});
        else
            % try *_config_*.m files, if toolbox does not have '*_cfg_*.m' files
            f2 = regexp(di,'.*_config_.*\.m$');
            if ~iscell(f2), f2 = {f2}; end
            fi = {di{~cellfun('isempty',f2)}};
            ftc = {ftc{:} fi{:}};
            dtc(end+1:end+length(fi)) = deal({d2});
        end;        
    end
    if ~isempty(ft)||~isempty(ftc)
        % The toolbox developer MUST add path to his/her toolbox in his/her 'prog'
        % function, with a command line like:
        % >> addpath(fullfile(spm('Dir'),'toolbox','mytoolbox'),'-end');
        cwd = pwd;
        j = 1;
        for i=1:length(ft)
            try
                cd(dt{i});
                tools.values{j} = feval(strtok(ft{i},'.'));
                j = j + 1;
            catch
                disp(['Loading of toolbox ' fullfile(dt{i},ft{i}) ' failed.']);
            end
        end
        for i=1:length(ftc)
            try
                cd(dtc{i});
                % use cfg_struct2cfg to convert from SPM5 to matlabbatch
                % configuration tree
                tools.values{j} = cfg_struct2cfg(feval(strtok(ftc{i},'.')));
                j = j + 1;
            catch
                disp(['Loading of toolbox ' fullfile(dtc{i},ftc{i}) ' failed.']);
            end
        end
        cd(cwd);
    end
end
%_______________________________________________________________________
% spmjobs SPM
%_______________________________________________________________________
spmjobs         = cfg_choice;
spmjobs.tag     = 'spm';
spmjobs.name    = 'SPM';
spmjobs.help    = {
                '%* Menu and Toolbar'
                '/*\subsection*{Menu and Toolbar}*/'
                'The "File" and "Edit" menu offer options to load, save and run a job and to modify the configuration of the batch system. For each application which is known to the batch system, a separate pulldown menu lists the available modules. Depending on the application, these modules may be grouped into submenus. Application specific defaults can be edited by choosing "Edit Defaults" from the application menu. The toolbar offers some shortcuts to frequently used operations (e.g. load, save, run).'
                'Jobs are saved as MATLAB .m files. These files contain a MATLAB script, which can be executed in MATLAB to recreate the job variable. Multiple jobs can be loaded at once. This allows to concatenate parts of a job.'
                ''
                '%* Top Left Panel'
                '/*\subsection*{Top Left Panel}*/'
                'The current job, which is represented as a list of executable modules. Modules marked with DEP depend on the successful execution of other modules in the job. Modules marked with X still require some values to be set before the job can be run, although an incompletely specified job can still be saved and loaded.'
                ''
                '%* Top Right Panel'
                '/*\subsection*{Top Right Panel}*/'
                'These are the configuration details for the currently selected module. Items marked with DEP depend on the successful execution of other modules in the job. Items marked with X still require some values to be set before the job can be run. Depending on the kind of detail currently selected, a choice of buttons appears below the Centre Right Panel to manipulate the current value.'
                ''
                '%* Centre Right Panel'
                '/*\subsection*{Centre Right Panel}*/'
                'This panel shows the current value of the highlighted item (where relevant).'
                ''
                '%* Edit Buttons'
                '/*\subsection*{Edit Buttons}*/'
                'Depending on the type of configuration item, different edit buttons appear.'
                '/*\begin{description}*/'
                '/*\item[Files]*/'
                '%* Files'
                '"Select Files" opens a file selection dialog box to select multiple files. "Edit Value" opens a generic value edit dialog to edit the list of files. "Dependencies" offers a list of outputs from other modules that can be used as an input to this item.'
                '/*\item[Generic Value]*/'
                '%* Generic Value'
                '"Edit Value" opens a generic value edit dialog to edit the list of files. "Dependencies" offers a list of outputs from other modules that can be used as an input to this item.'
                '%* Menu'
                '/*\item[Menu]*/'
                '"Edit Value" opens a selection dialog showing allowed menu options.'
                '%* Choice'
                '/*\item[Choice]*/'
                '"Edit Value" opens a selection dialog showing allowed menu options. Depending on the choosen option the module configuration may change.'
                '%* Repeat'
                '/*\item[Repeat]*/'
                '"Add Item", "Replicate Item", "Delete Item" allow to add new repeated items, to replicate or to delete items from the list. If more than one item or item type exists, a dialog popup will appear listing the available options. Multiple selections are allowed.'
                '/*\end{description}*/'
                ''
                '%* Bottom Panel'
                '/*\subsection*{Bottom Panel}*/'
                'This panel provides information about the meaning of the current item.'
                '/*\begin{figure} \begin{center} \includegraphics[width=70mm]{images/batch_ui1} \includegraphics[width=70mm]{images/batch_ui2} \includegraphics[width=70mm]{images/ui3} \includegraphics[width=70mm]{images/ui4}\end{center} \caption{The SPM5 user interface. \emph{Top left:} The usual user-interface.  \emph{Top right:} The Defaults user-interface. \emph{Bottom left:} The file selector (click the (?) button for more information about filtering filenames, or selecting individual volumes within a 4D file). \emph{Bottom right:} more online help can be obtained via the main help button.} \end{figure} */'
}';
spmjobs.values  = { temporal spatial stats meeg util tools};
