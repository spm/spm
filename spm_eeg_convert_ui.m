function spm_eeg_convert_ui
% User interface for M/EEG data conversion facility
% FORMAT spm_eeg_convert_ui
% 
% See spm_eeg_convert for details.
%__________________________________________________________________________
% Copyright (C) 2008-2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_convert_ui.m 5428 2013-04-18 17:34:49Z guillaume $

SVNrev = '$Rev: 5428 $';

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
    spm_eeg_convert(fname);
end

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','');
