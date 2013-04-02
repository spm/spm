function delete = spm_cfg_eeg_delete
% configuration file for deleting M/EEG datasets
%__________________________________________________________________________
% Copyright (C) 2009-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_delete.m 5377 2013-04-02 17:07:57Z vladimir $

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
    delete(spm_eeg_load(job.D{i}));
end

%==========================================================================
