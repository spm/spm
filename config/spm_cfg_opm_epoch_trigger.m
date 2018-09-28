function epoch = spm_cfg_opm_epoch_trigger
% configuration file for epoching OPM data
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_cfg_opm_epoch_trigger.m 7429 2018-09-28 09:29:20Z spm $

%--------------------------------------------------------------------------
% Output Directory
%--------------------------------------------------------------------------
D = cfg_files;
D.tag = 'D';
D.name = 'File Name';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the M/EEG mat file.'};

%--------------------------------------------------------------------------
% Output filename
%--------------------------------------------------------------------------
timewin  = cfg_entry;
timewin.tag     = 'timewin';
timewin.name    = 'Time window';
timewin.help    = {'Time window around trigger which to epoch(ms). e.g [-200,300]'};
timewin.strtype = 'r';
timewin.num     = [Inf,2];

%--------------------------------------------------------------------------
% Condition Labels
%--------------------------------------------------------------------------
condLabels         = cfg_entry;
condLabels.tag     = 'condLabels';
condLabels.name    = 'Condition Labels';
condLabels.help    = {'Labels of conditions. Enter each label on a new line. If left empty the default behaviout is to label conditions accrding to number(e.g. Cond1,Cond2,...)'};
condLabels.strtype = 's+';
condLabels.num     = [1,100];
condLabels.val     = {{''}};


%--------------------------------------------------------------------------
% simulation parameters
%--------------------------------------------------------------------------
epoch          = cfg_exbranch;
epoch.tag      = 'epoch';
epoch.name     = 'Epoch M/EEG object on Trigger';
epoch.val      = {D,timewin,condLabels};
epoch.help     = {'Epoch M/EEG data at the rise of every trigger in the dataset'}';
epoch.prog     = @epoch_trigger;
epoch.vout     = @vout_epoch_trigger;
epoch.modality = {'EEG'};


%==========================================================================
function out = epoch_trigger(job)
% construct the S struct

% datset parameters
S=[];
S.D= spm_eeg_load(job.D{1});
S.timewin= job.timewin;
S.condLabels = job.condLabels;

% check if condLabels is empty
defaultLabels = strcmp(S.condLabels{1},'');

% remove the unused fields
if(defaultLabels)
S = rmfield(S,'condLabels');
end

% run the main function 
out.D= spm_opm_epoch_trigger(S);
out.Dfname = {fullfile(out.D.path, out.D.fname)};


%==========================================================================
function dep = vout_epoch_trigger(job)
% return dependencies
dep = cfg_dep;
dep.sname = 'Epoched Data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2) = cfg_dep;
dep(2).sname = 'Epoched Datafile';
% reference field "Dfname" from output
dep(2).src_output = substruct('.','Dfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
