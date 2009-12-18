function spm_Beamforming
% GUI gateway to Beamforming toolbox
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_Beamforming.m 3652 2009-12-18 18:54:43Z guillaume $

funlist = {
    'Volumetric LCMV beamformer',             'spm_eeg_ft_beamformer_gui'
    'Fieldtrip DICS beamformer',              'spm_eeg_ft_beamformer_freq'
    'Fieldtrip beamformer source extraction', 'spm_eeg_ft_beamformer_source'
};

str = sprintf('%s|', funlist{:, 1});
str = str(1:(end-1));

fun = spm_input('Beamforming',1,'m', str, strvcat(funlist(:, 2)));
  
eval(fun);