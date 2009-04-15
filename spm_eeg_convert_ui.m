function spm_eeg_convert_ui(S)
% User interface for M/EEG data conversion facility
% FORMAT spm_eeg_convert_ui(S)
% S       - existing configuration struct (optional)
%__________________________________________________________________________
% 
% See spm_eeg_convert for a description of input structure S.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_convert_ui.m 3059 2009-04-15 18:09:13Z guillaume $

SVNrev = '$Rev: 3059 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FnUIsetup','M/EEG data conversion',0);

%-Get parameters
%--------------------------------------------------------------------------
try
    S.dataset;
catch
    [S.dataset, sts] = spm_select(1, '.*', 'Select M/EEG data file');
    if ~sts, return; end
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
            S.timewindow = spm_input('Input time window ([start end] in sec)', '+1', 'r', '', 2);
        end
    else
        res = spm_input('Where to look for trials?','+1','data|define|file',[1 2 3], 1);

        switch res
            case 1
                S.usetrials = 1;
            case 2
                S.usetrials = 0;
                [trl, conditionlabel, S] = spm_eeg_definetrial(S);
                S.trl = trl;
                S.conditionlabel = conditionlabel;
            case 3
                S.usetrials = 0;
                if  ~isfield(S, 'trlfile')
                    S.trlfile = spm_select(1, 'mat', 'Select a trial definition file');
                end
        end
    end

    if S.continuous || ~S.usetrials
         S.checkboundary = spm_input('Read across trial borders?','+1','yes|no',[0 1]);
    end
    
    if ~isfield(S, 'channels')
        S.channels = spm_input('What channels?','+1','all|meg|eeg|gui|file');
    end

    S = spm_eeg_channelselection(S);

    if ~isfield(S, 'outfile')
        if S.continuous
            prefix = 'spm8_';
        else
            prefix = 'espm8_';
        end
        S.outfile = spm_input('SPM M/EEG file name', '+1', 's', [prefix spm_str_manip(S.dataset,'tr')]);
    end
    
end

%-Call the conversion routine
%--------------------------------------------------------------------------
spm('Pointer','Watch');
D = spm_eeg_convert(S);
spm('Pointer','Arrow');

%-Display the imported M/EEG data file
%--------------------------------------------------------------------------
if ~isfield(S, 'review') || S.review
    spm_eeg_review(D);
end
