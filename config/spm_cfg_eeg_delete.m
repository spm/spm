function delete = spm_cfg_eeg_delete
% Configuration file for deleting M/EEG datasets
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2009-2022 Wellcome Centre for Human Neuroimaging


%--------------------------------------------------------------------------
% D
%--------------------------------------------------------------------------
D        = cfg_files;
D.tag    = 'D';
D.name   = 'File Name';
D.filter = 'mat';
D.num    = [1 Inf];
D.help   = {'Select the M/EEG mat file(s).'};

%--------------------------------------------------------------------------
% delete
%--------------------------------------------------------------------------
delete          = cfg_exbranch;
delete.tag      = 'delete';
delete.name     = 'Delete';
delete.val      = {D};
delete.help     = {'Deleting M/EEG datasets'}';
delete.prog     = @eeg_delete;
delete.modality = {'EEG'};


%==========================================================================
function eeg_delete(job)

for i = 1:numel(job.D)
    try
        delete(spm_eeg_load(job.D{i}));
    end
end
