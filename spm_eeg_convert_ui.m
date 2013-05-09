function spm_eeg_convert_ui
% User interface for M/EEG data conversion facility
% FORMAT spm_eeg_convert_ui
% 
% See spm_eeg_convert for details.
%__________________________________________________________________________
% Copyright (C) 2008-2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_convert_ui.m 5482 2013-05-09 10:11:51Z vladimir $

SVNrev = '$Rev: 5482 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FnUIsetup','M/EEG data conversion',0);

[fname, sts] = spm_select(1, '.*', 'Select M/EEG data file');
if ~sts, return; end

if spm_input('Define settings?','+1','yes|just read',[1 0], 0);
    matlabbatch{1}.spm.meeg.convert.dataset = {fname};
    spm_jobman('interactive', matlabbatch);
else
    D = spm_eeg_convert(fname);
    spm_eeg_review(D);
end

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','');
