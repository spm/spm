function S = spm_cfg_eeg_review
% configuration file for M/EEG reviewing tool
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_review.m 5377 2013-04-02 17:07:57Z vladimir $

rev = '$Rev: 5377 $';

D = cfg_files;
D.tag = 'D';
D.name = 'File Name';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the EEG mat file.'};


S = cfg_exbranch;
S.tag = 'review';
S.name = 'Display';
S.val = {D};
S.help = {'Run the reviewing tool with the given dataset as input.'};
S.prog = @eeg_review;
S.modality = {'EEG'};

%==========================================================================
function eeg_review(job)

spm_eeg_review(spm_eeg_load(job.D{1}));


