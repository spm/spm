function grandmean = spm_cfg_eeg_grandmean
% Configuration file for averaging evoked responses
%__________________________________________________________________________

% Stefan Kiebel
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


%--------------------------------------------------------------------------
% S Grandmean
%--------------------------------------------------------------------------
grandmean          = cfg_exbranch;
grandmean.tag      = 'grandmean';
grandmean.name     = 'Grandmean';
grandmean.val      = @grandmean_cfg;
grandmean.help     = {'Average multiple evoked responses'};
grandmean.prog     = @eeg_grandmean;
grandmean.vout     = @vout_eeg_grandmean;
grandmean.modality = {'EEG'};


%==========================================================================
function varargout = grandmean_cfg

persistent cfg
if ~isempty(cfg), varargout = {cfg}; return; end

%--------------------------------------------------------------------------
% D File Names
%--------------------------------------------------------------------------
D        = cfg_files;
D.tag    = 'D';
D.name   = 'File Names';
D.filter = 'mat';
D.num    = [1 inf];
D.help   = {'Select the M/EEG mat file.'};

%--------------------------------------------------------------------------
% Output filename
%--------------------------------------------------------------------------
outfile         = cfg_entry;
outfile.tag     = 'outfile';
outfile.name    = 'Output filename';
outfile.strtype = 's';
outfile.num     = [1 inf];
outfile.help    = {'Choose filename'};

% -------------------------------------------------------------------------
% weighted Weighted average
% -------------------------------------------------------------------------
weighted        = cfg_menu;
weighted.tag    = 'weighted';
weighted.name   = 'Weighted average?';
weighted.help   = {'Average weighted by number of replications in input.'};
weighted.labels = {'Yes' 'No'};
weighted.values = {1 0};

[cfg,varargout{1}] = deal({D outfile weighted});


%==========================================================================
function out = eeg_grandmean(job)
% construct the S struct
S = job;
S.D = strvcat(S.D);

out.D = spm_eeg_grandmean(S);
out.Dfname = {fullfile(out.D)};

%==========================================================================
function dep = vout_eeg_grandmean(job)

dep(1)            = cfg_dep;
dep(1).sname      = 'Grandmean Data';
dep(1).src_output = substruct('.','D');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'Grandmean Datafile';
dep(2).src_output = substruct('.','Dfname');
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
