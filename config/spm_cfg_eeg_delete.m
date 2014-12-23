function delete = spm_cfg_eeg_delete
% Configuration file for deleting M/EEG datasets
%__________________________________________________________________________
% Copyright (C) 2009-2014 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_delete.m 6293 2014-12-23 18:15:57Z guillaume $


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
delete.help     = {'Copying M/EEG datasets'}';
delete.prog     = @eeg_delete;
delete.modality = {'EEG'};


%==========================================================================
function eeg_delete(job)

for i = 1:numel(job.D)
    try
        delete(spm_eeg_load(job.D{i}));
    end
end
