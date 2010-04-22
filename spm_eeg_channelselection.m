function S = spm_eeg_channelselection(S)
% Function for selection of channels
% FORMAT S = spm_eeg_channelselection(S)
% S - existing configuration struct (optional)
% Fields of S:
% S.channels - can be 'MEG', 'EEG', 'file', 'gui' or cell array of labels
% S.chanfile - filename (used in case S.channels = 'file')
% S.dataset - MEEG dataset name
% S.inputformat - data type (optional) to force the use of specific data
% reader
%
% OUTPUT:
%   S.channels - cell array of labels
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_channelselection.m 3833 2010-04-22 14:49:48Z vladimir $

SVNrev = '$Rev: 3833 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG channels selection');

%-Get parameters
%--------------------------------------------------------------------------
try
    S.channels;
catch
    S.channels = 'gui';
end

try
    S.dataset;
catch
    [S.dataset, sts] = spm_select(1, '.*', 'Select M/EEG data file');
    if ~sts, return; end
end

try
    S.inputformat;
catch
    S.inputformat = [];
end

hdr = ft_read_header(S.dataset, 'fallback', 'biosig', 'headerformat', S.inputformat);

if strcmp(S.channels, 'file')
    if ~isfield(S, 'chanfile')
        [S.chanfile, sts] = spm_select(1, 'mat', 'Select channel selection file');
        if ~sts, return; end
    end
    label = load(S.chanfile, 'label');
    if ~isfield(label, 'label')
        error('Channel selection file does not contain labels.');
    else
        S.channels = label.label;
    end
    S.save = 0;
else
    switch S.channels
        case 'eeg'
            S.channels = 'EEG';
        case 'meg'
            S.channels = 'MEG';
        otherwise
    end

    % Make sure the order is like in the file
    channels     = ft_channelselection(S.channels, hdr.label);
    [sel1, sel2] = spm_match_str(hdr.label, channels);
    S.channels   = channels(sel2);
end

%-Create trial definition file
%--------------------------------------------------------------------------
if ~isfield(S, 'save')
    S.save = spm_input('Save channel selection?','+1','yes|no',[1 0], 0);
end

if S.save
    [chanfilename, chanpathname] = uiputfile( ...
        {'*.mat', 'MATLAB File (*.mat)'}, 'Save channel selection as');

    label = S.channels;

    save(fullfile(chanpathname, chanfilename), 'label');
end
