function [D] = spm_eeg_invert_ui(varargin)
% GUI for ReML inversion of forward model for EEG-MEG
% FORMAT [D] = spm_eeg_invert_ui(D,val)
% ReML estimation of regularisation hyperparameters using the
% spatio-temporal hierarchy implicit in EEG data
% sets:
%
%     D.inv{i}.inverse.trials - trials (in D.events.types) to invert
%     D.inv{i}.inverse.smooth - smoothness of source priors (mm)
%     D.inv{i}.inverse.type   - 'MSP' multiple sparse priors
%                               'LOR' LORETA-like model
%                               'IID' LORETA and WMN
%     D.inv{i}.inverse.xyz    - (n x 3) locations of spherical VOIs
%     D.inv{i}.inverse.rad    - radius (mm) of VOIs
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_eeg_invert_ui.m 6849 2016-07-31 12:34:33Z karl $

% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});


% check whether to use conventional or DCM temporal priors
%--------------------------------------------------------------------------
q_rec = spm_input('Reconstruction','+1','b',{'Imaging|VB-ECD|DCM'},[0 1 2],1);
switch q_rec
    case 2
        % record type in D and DCM structures
        %==================================================================
        inverse.type = 'DCM';

        % exchange filenames
        %------------------------------------------------------------------
        DCMfile            = ['DCM_' D.fname];
        D.inv{val}.DCMfile = DCMfile;
        DCM.val            = val;
        DCM.xY.Dfile       = fullfile(D.path,D.fname);
        DCM.options.type   = 2;
        DCM.name           = DCMfile;
        
        % an call API to specify DCM
        %------------------------------------------------------------------
        spm_api_erp(DCM);
        D.inv{val}.inverse = inverse;
        
    case 1
        % Use Variational Bayes Equivalent Current Dipole reconstruction
        %==================================================================
        D = spm_eeg_inv_vb_ecd_gui(D,val);
        
    case 0
        % Conventional imaging reconstruction: get conditions or trials
        %==================================================================
        if D.nconditions > 1
            if spm_input('All conditions or trials','+1','b',{'yes|no'},[1 0],1)
                trials = D.condlist;
            else
                trials = {};
                condlabels = D.condlist;
                for  i = 1:D.nconditions
                    str = sprintf('invert %s', condlabels{i});
                    if spm_input(str,'+1','b',{'yes|no'},[1 0],1);
                        trials{end + 1} = condlabels{i};
                    end
                end
            end
        else
            trials = D.condlist;
        end
        
        % Inversion parameters
        %------------------------------------------------------------------
        inverse        = spm_eeg_inv_custom_ui(D);
        inverse.trials = trials;

                
        % invert
        %==================================================================
        D.con               = 1;
        D.inv{val}.inverse  = inverse;
        
        % Modality
        %------------------------------------------------------------------
        [mod, list] = modality(D, 1, 1);
        if strcmp(mod, 'Multimodal')
            [selection, ok]= listdlg('ListString', list, 'SelectionMode', 'multiple' ,...
            'Name', 'Select modalities' , 'InitialValue', 1:numel(list),  'ListSize', [400 300]);
            if ~ok
                return;
            end
            
            D.inv{val}.inverse.modality  = list(selection);
            
            if numel(D.inv{val}.inverse.modality) == 1
                D.inv{val}.inverse.modality = D.inv{val}.inverse.modality{1};
            end
        else
            D.inv{val}.inverse.modality = mod;
        end
        
        D                            = spm_eeg_invert(D);
end
