function [D] = spm_eeg_invert_ui(D)
% GUI for ReML inversion of forward model for EEG-EMG
% FORMAT [D] = spm_eeg_invert_ui(D)
% ReML estimation of regularisation hyperparameters using the
% spatio-temporal hierarchy implicit in EEG data
% sets:
%
%     D.inverse.con    - condition or trial type
%     D.inverse.smooth - smoothness of source priors (mm)
%     D.inverse.type   - 'MSP' multiple sparse priors
%                        'LOR' LORETA-like model
%                        'IID' LORETA and WMN
%__________________________________________________________________________


% load data if necessary
%--------------------------------------------------------------------------
try
    D;
catch
    D = load(spm_select(1, '.mat','Select EEG/MEG mat file'));
    D = D.D;
end

try, val = D.val; catch, val = length(D.inv); end

try
    D.inv{val}.forward.gainmat;
catch
    warndlg('Please create forward model first');
    return
end

% Get condition or trial type
%--------------------------------------------------------------------------
if D.events.Ntypes > 1
    con = spm_input(sprintf('1:%i Condition or trial',D.events.Ntypes),'+1');

    try
        D.events.types(con);
    catch
        warndlg('please select a single trial');
        return
    end
else
    con = 1;
end
D.inv{val}.inverse.con = con;

% Type of analysis
%--------------------------------------------------------------------------
D.inv{val}.inverse.type   = spm_input('Type of inversion','+1','MSP|LOR|IID');

% D.inverse.smooth - smoothness of source priors (mm)
%--------------------------------------------------------------------------
D.inv{val}.inverse.smooth = spm_input('Smoothness (0-1)','+1','r',0.5);

% Number of sparse priors
%--------------------------------------------------------------------------
if strcmp(D.inv{val}.inverse.type,'MSP')
    D.inv{val}.inverse.Np = ...
    spm_input('Number of MSPs','+1','128|256|512',[128 256 512]);
end

% invert
%--------------------------------------------------------------------------
D = spm_eeg_invert(D);



