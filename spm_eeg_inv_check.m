function [D,val] = spm_eeg_inv_check(varargin)

%==========================================================================
% Checks that the EEG/EMG .mat file structure is loaded properly and that
% the particular inversion of interest has been specified
%
% FORMAT [D,val] = spm_eeg_inv_check(D,[val])
% Input:
% S              - data structure or its filename
% val            - model of interest (usually 1)
% Output:
% D              - data structure
% val            - model of interest D.val
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout, Karl Friston
% $Id: spm_eeg_inv_check.m 621 2006-09-12 17:22:42Z karl $


% Check - prompt for file if necessary
%--------------------------------------------------------------------------
try
    D = varargin{1};
    D.inv;
catch
    try
        D = spm_eeg_ldata(D);
        D.inv;
    catch
        D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
        D = spm_eeg_ldata(D);
        
        % Check for inversion
        %------------------------------------------------------------------
        if ~isfield(D,'inv')
            warndlg('No inverse structure has been created yet');
            return
        end
    end
end

% set val = 1 if only one model
%--------------------------------------------------------------------------
if length(D.inv) == 1
    val   = 1;
    D.val = val;
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
D.val = val;

