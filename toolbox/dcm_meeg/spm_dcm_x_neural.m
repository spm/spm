function [x,f] = spm_dcm_x_neural(P,model)
% Returns the state and equation of neural mass models
% FORMAT [x,f] = spm_dcm_x_neural(P,'model')
%
%  P      - parameter structure
% 'model'   - 'ERP','SEP','CMC','LFP','NNM' or 'MFM'
%
% x   - initial states
% f   - state euquation
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_x_neural.m 4232 2011-03-07 21:01:16Z karl $
 
 
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
        f  = 'spm_fx_sep';
        
    % Linear in states – canonical microcircuit
    %======================================================================
    case{'cmc'}
 
        % inital states
        %------------------------------------------------------------------
        x  =  spm_x_cmc(P);
        f  = 'spm_fx_cmc';
        
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
        x  = spm_x_nmm(P);
        f  = 'spm_fx_mfm';
 
 
    % Mean field model (nonlinear in states) - with covariance
    %======================================================================
    case{'mfm'}
 
        % inital states and model
        %------------------------------------------------------------------
        x  = spm_x_mfm(P);
        f  = 'spm_fx_mfm';
        
 
    otherwise
        warndlg('Unknown model')
end
