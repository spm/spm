function [D,ival] = spm_eeg_inv_check(varargin)
% Checks that the EEG/EMG .mat file structure is loaded properly and that
% the particular inversion of interest has been specified
%
% FORMAT [D,ival] = spm_eeg_inv_check(D,[ival])
% Input:
% S              - data structure or its filename
% ival            - model of interest (usually 1)
% Output:
% D              - data structure
% ival            - model of interest D.ival
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout, Karl Friston
% $Id: spm_eeg_inv_check.m 1437 2008-04-17 10:34:39Z christophe $


% Check - prompt for file if necessary
%--------------------------------------------------------------------------
try
    D = varargin{1};
    D.inv;
catch
    try
        D = spm_eeg_load(D);
        D.inv;
    catch
        D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
        D = spm_eeg_load(D);
        
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
    ival   = 1;
    D.ival = ival;
    return
end

% Check - val
%--------------------------------------------------------------------------
try
    ival = varargin{2};
catch
    try
        ival    = D.ival;
    catch
        prompt = sprintf('which model (1 to %i)',length(D.inv));
        ival    = inputdlg(prompt,'Source reconstruction',1,{'1'});
        ival    = eval(ival{1});
    end
end
D.ival = ival;

