function [x,f] = spm_dcm_x_neural_NMDA(P,model)
% Return the state and equation of neural mass models
% FORMAT [x,f] = spm_dcm_x_neural_NMDA(P,'model')
%
%  P      - parameter structure
% 'model' - 'ERP','SEP','LFP','NNM' or 'MFM'
%
% x   - initial states
% f   - state euquation
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% initial state and equation
%--------------------------------------------------------------------------
switch lower(model)

    % linear David et al model (linear in states)
    %======================================================================
    case{'erp'}

        % inital states and equations of motion
        %------------------------------------------------------------------
        x  =  spm_x_erp(P);
        f  = 'spm_fx_erp';


    % linear David et al model (linear in states) - fast version for SEPs
    %======================================================================
    case{'sep'}

        % inital states
        %------------------------------------------------------------------
        x  =  spm_x_erp(P);
        f  = 'spm_fx_erp';
        
    % linear David et al model (linear in states) - with self-inhibition
    %======================================================================
    case{'lfp'}

        % inital states
        %------------------------------------------------------------------
        x  =  spm_x_lfp(P);
        f  = 'spm_fx_lfp';


    % Neural mass model (nonlinear in states)
    %======================================================================
    case{'nmm'}

        % inital states and model
        %------------------------------------------------------------------
        x  = spm_x_nmm_NMDA(P);
        f  = 'spm_fx_mfm_NMDA';


    % Mean field model (nonlinear in states) - with covariance
    %======================================================================
    case{'mfm'}

        % inital states and model
        %------------------------------------------------------------------
        x  = spm_x_mfm_NMDA(P);
        f  = 'spm_fx_mfm_NMDA';
        

    otherwise
        warndlg('Unknown model')
end
