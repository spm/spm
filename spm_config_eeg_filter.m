function S = spm_config_eeg_filter
% configuration file for EEG Filtering
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: spm_config_eeg_filter.m 112 2005-05-04 18:20:52Z john $


D = struct('type','files','name','File Name','tag','D',...
           'filter','mat','num',1,...
           'help',{{'Select the EEG mat file.'}});

typ = struct('type','menu','name','Filter type','tag','type',...
           'labels',{{'Butterworth'}},'values',{{1}},'val',{{1}},...
           'help', {{'Select the filter type.'}});

PHz = struct('type','entry','name','Cutoff','tag','cutoff',...
           'strtype','r','num',[1 1],...
           'help', {{'Enter the filter cutoff'}});
flt = struct('type','branch','name','Filter','tag','filter','val',{{typ,PHz}});

S   = struct('type','branch','name','EEG Filter','tag','eeg_filter',...
           'val',{{D,flt}},'prog',@eeg_filter,'modality',{{'EEG'}},...
           'help',{{'Low-pass filters EEG/MEG epoched data.'}});

function eeg_filter(S)
S.D = strvcat(S.D{:});
spm_eeg_filter(S)
