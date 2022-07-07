function spm_eeg_convert_ui
% User interface for M/EEG data conversion facility
% FORMAT spm_eeg_convert_ui
% 
% See spm_eeg_convert for details.
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename);
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
