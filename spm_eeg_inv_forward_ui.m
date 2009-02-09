function D = spm_eeg_inv_forward_ui(varargin)
% Forward Solution user-interface routine
% commands the forward computation for either EEG or MEG data
% and calls for various types of solutions using BrainStorm functions
% as well as a realistic sphere solution (for EEG).
%
% FORMAT D = spm_eeg_inv_forward_ui(D,val)
% Input:
% D         - input data struct (optional)
% Output:
% D         - same data struct including the forward solution files and variables
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_forward_ui.m 2720 2009-02-09 19:50:46Z vladimir $

% initialise
%--------------------------------------------------------------------------
[D, val] = spm_eeg_inv_check(varargin{:});


D.inv{val}.forward = struct([]);

for i = 1:numel(D.inv{val}.datareg)
    switch D.inv{val}.datareg(i).modality
        case 'EEG'
            if D.inv{val}.mesh.template
                models = {'EEG BEM', '3-Shell Sphere (experimental)'};
            else
                models = {'EEG BEM', '3-Shell Sphere (experimental)'};
            end
        case 'MEG'
            models = {'Single Sphere', 'MEG Local Spheres', 'Single Shell'};
        otherwise
            error('Unsupported modality');
    end
    str = sprintf('%s|', models{:});
    str = str(1:(end-1));

    D.inv{val}.forward(i).voltype  = spm_input(sprintf('Which %s head model?', D.inv{val}.datareg(i).modality), 1, 'm', str, strvcat(models));
end

% compute forward model
%==========================================================================
D = spm_eeg_inv_forward(D);

spm_eeg_inv_checkforward(D, val);


fprintf('Foward model complete - thank you\n')
