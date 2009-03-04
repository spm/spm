function Heeg = spm_eeg_display_ui(varargin)
% user interface for displaying EEG/MEG channel data.
% Heeg = spm_eeg_display_ui(varargin)
%
% optional argument:
% S         - struct
%    fields of S:
%     D       - MEEG object
%     Hfig    - Figure (or axes) to work in (Defaults to SPM graphics window)
%     rebuild - indicator variable: if defined spm_eeg_display_ui has been
%                                   called after channel selection
%
% output:
%     Heeg      - Handle of resulting figure
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_display_ui.m 2826 2009-03-04 17:24:49Z james $

if nargin == 1
    S = varargin{1};
    try
        D = S.D;
    end
end
if ~exist('D','var')
    D = spm_select(1, '\.mat$', 'Select EEG mat file');
    if ~isempty(D)
        try
            D = spm_eeg_load(D);
            sD = struct(D); % transform to struct for access to some fields
        catch
            error(sprintf('Trouble reading file %s', D));
        end
    end
end

if ~isempty(D)
    spm_eeg_review(D,5)
end


