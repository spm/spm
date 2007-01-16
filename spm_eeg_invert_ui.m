function [D] = spm_eeg_invert_ui(varargin)
% GUI for ReML inversion of forward model for EEG-EMG
% FORMAT [D] = spm_eeg_invert_ui(D,val)
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


% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});

% Get condition or trial type
%--------------------------------------------------------------------------
if D.events.Ntypes > 1
    con   = 1:D.events.Ntypes;
    for i = con;
        str{i} = num2str(i);
    end
    con = spm_input('Condition or trial','+1','b',str,con,1);
    try
        D.events.types(con);
    catch
        warndlg('please select a single trial');
        return
    end
else
    con = 1;
end
D.inv{val}.inverse.con    = con;

% Type of analysis
%--------------------------------------------------------------------------
D.inv{val}.inverse.type   = ...
        spm_input('Type of inversion','+1','MSP|LOR|IID');

% D.inverse.smooth - smoothness of source priors (mm)
%--------------------------------------------------------------------------
D.inv{val}.inverse.smooth = ...
        spm_input('Smoothness (0-1)','+1','0.2|0.4|0.6',[0.2 0.4 0.6],2);

% Number of sparse priors
%--------------------------------------------------------------------------
if strcmp(D.inv{val}.inverse.type,'MSP')
    D.inv{val}.inverse.Np = ...
        spm_input('Number of MSPs','+1','128|192|256',[128 192 256]);

    if spm_input('Restrict solutions','+1','yes|no',[1 0],2);
        D.inv{val}.inverse.xyzr = ...
        spm_input('x,y,z and radius of spheres (mm)','0','r',[0 0 0 64]);
    end
end

% invert
%--------------------------------------------------------------------------
D = spm_eeg_invert(D);



