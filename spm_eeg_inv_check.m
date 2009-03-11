function [D,val] = spm_eeg_inv_check(varargin)
% Checks that the EEG/MEG .mat file structure is loaded properly and that
% the particular inversion of interest has been specified
%
% FORMAT [D,val] = spm_eeg_inv_check(D,[val])
% Input:
% S              - data structure or its filename
% val            - model of interest (usually 1)
% Output:
% D              - data structure
% val            - model of interest D.val
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout, Karl Friston
% $Id: spm_eeg_inv_check.m 2861 2009-03-11 18:41:03Z guillaume $


% Check - prompt for file if necessary
%--------------------------------------------------------------------------
if nargin == 0
    [D, sts] = spm_select(1, 'mat', 'Select EEG/MEG mat file');
    if ~sts, D = []; val = 0; return; end
else
    D = varargin{1};
end

D = spm_eeg_load(D);

% Check for inversion
%--------------------------------------------------------------------------
if ~isfield(D,'inv')
    val = 0;
    return
end

% Set val = 1 if only one model
%--------------------------------------------------------------------------
if length(D.inv) == 1
    val   = 1;
    D.val = val;
    return
end

% Check - val
%--------------------------------------------------------------------------
try
    val        = varargin{2};
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
