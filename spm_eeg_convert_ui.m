function spm_eeg_convert_ui(S)
% User interface for M/EEG conversion function
% FORMAT spm_eeg_convert(S)
% S - existing configuration struct (optional)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_convert_ui.m 1240 2008-03-25 16:34:30Z vladimir $
if nargin == 0
    S=[];
end

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','MEEG data conversion ',0);

if ~isfield(S, 'dataset')
    S.dataset = spm_select(1, '\.*', 'Select M/EEG data file');
end

if spm_input('Define settings?','+1','yes|just read',[1 0], 0);

    if ~isfield(S, 'continuous')
        S.continuous = spm_input('How to read?','+1','continuous|trials',[1 0], 1);
    end

    if S.continuous
        readall = spm_input('Read everything?','+1','yes|no',[1 0]);

        if ~isfield(S, 'timewindow')
            S.timewindow = [];
        end

        if ~readall  &&  isempty(S.timewindow);
            S.timewindow = spm_input('Input time window ([start end] in sec)', '+1', 'r');
        end
    else
        S.usetrials = spm_input('Where to look for trials?','+1','data|file',[1 0], 1);

        if  ~S.usetrials && ~isfield(S, 'trlfile')
            S.trlfile = spm_select(1, '\.mat$', 'Select a trial definition file');
        end
    end

    if ~isfield(S, 'allchannels')
        S.allchannels = spm_input('Read all channels?','+1','yes|no',[1 0], 1);
    end

    if ~S.allchannels && ~isfield(S, 'chanfile')
        S.chanfile = spm_select(1, '\.mat$', 'Select a channel selection file');
    end

    if ~isfield(S, 'outfile')
        S.outfile = spm_input('SPM EEG file name', '+1', 's', spm_str_manip(S.dataset,'tr'));
    end
    
end

spm_eeg_convert(S);
