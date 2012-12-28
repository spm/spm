function spm_eeg_convert_ui
% User interface for M/EEG data conversion facility
% FORMAT spm_eeg_convert_ui
%__________________________________________________________________________
% 
% See spm_eeg_convert for a description of input structure S.
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_convert_ui.m 5164 2012-12-28 16:40:06Z vladimir $

SVNrev = '$Rev: 5164 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FnUIsetup','M/EEG data conversion',0);


if spm_input('Define settings?','+1','yes|just read',[1 0], 0);
    spm_jobman('interactive','','spm.meeg.convert', dataset);
else
    [fname, sts] = spm_select(1, '.*', 'Select M/EEG data file');
    if sts,
        spm_eeg_convert(fname);
    end
end

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName',''); 