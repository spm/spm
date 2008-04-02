function S = spm_eeg_channelselection(S)
% Function for selection of channels
% FORMAT S = spm_eeg_channelselection(S)
% S - existing configuration struct (optional)
% Fields of S:
% S.channels - can be 'MEG', 'EEG', 'file', 'gui' or cell array of labels
% S.chanfile - filename (used in case S.channels = 'file')
% S.dataset - MEEG dataset name
%
% OUTPUT:
%   S.channels - cell array of labels
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_channelselection.m 1291 2008-04-02 13:58:28Z vladimir $

if nargin == 0
    S = [];
end

% ------------- Check inputs

Fig = spm_figure('GetWin','Interactive');
clf(Fig);

if ~isfield(S, 'channels') S.channels = 'gui'; end

if ~isfield(S, 'dataset')
    S.dataset = spm_select(1, '\.*', 'Select M/EEG data file');
end
hdr = read_header(S.dataset);

if strcmp(S.channels, 'file')
    if ~isfield(S, 'chanfile')
        S.chanfile = spm_select(1, '\.mat', 'Select channel selection file');
        label = load(S.chanfile, 'label');
        if ~isfield(label, 'label')
            error('Channel selection file does not contain labels.');
        else
            S.channels = label.label;
        end
    end
elseif exist('channelselection') == 2
    S.channels = channelselection(S.channels, hdr.label);
else
    error('Fieldtrip channelselection not found. Use a file to select channels');
end

% ------------- Create trial definition file
if ~isfield(S, 'save')
    S.save = spm_input('Save channel selection?','+1','yes|no',[1 0], 0);
end

if S.save
    [chanfilename, chanpathname] = uiputfile( ...
        {'*.mat', 'MATLAB File (*.mat)'}, 'Save channel selection as');

    label = S.channels;

    save(fullfile(chanpathname, chanfilename), 'label');
end

