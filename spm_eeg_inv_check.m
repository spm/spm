function [D,val] = spm_eeg_inv_check(varargin)

%==========================================================================
% Checks that the EEG/EMG .mat file structure is loaded properly and that
% the particular inversion of interat has been specified
%
% FORMAT [D,val] = spm_eeg_inv_check(D,[val])
% Input:
% S              - input data struct (optional)
% val            - Inversion of interest
% Output:
% D              - data struct including the new files and parameters
% val            - Inversion of interest D.val
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_check.m 621 2006-09-12 17:22:42Z karl $


% Check - prompt for file if necessary
%--------------------------------------------------------------------------
try
    D = varargin{1};
    D.inv;
catch
    D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
    D = spm_eeg_ldata(D);
end

% Check - inversion
%--------------------------------------------------------------------------
if ~isfield(D,'inv')
    warndlg('No inverse structure has been created yet');
    return
end

% Check - val
%--------------------------------------------------------------------------
try
    val = varargin{2};
catch
    try
        val    = D.val;
    catch
        prompt = sprintf('which model (1 to %i)',length(D.inv));
        val    = inputdlg(prompt,'Source reconstruction',1,{'1'});
        val    = eval(val{1});
    end
end

