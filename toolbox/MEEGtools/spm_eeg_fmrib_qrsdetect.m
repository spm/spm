function D = spm_eeg_fmrib_qrsdetect(S)
% Detect QRS complexes in the data using FMRIB plugin
% FORMAT  D = spm_eeg_fmrib_qrsdetect(S)
%
% S                    - input structure (optional)
% (optional) fields of S:
%   S.D                - MEEG object or filename of M/EEG mat-file
%   S.ecgchan          - name of index of the ECG channel (optional). If absent
%                        will try to look for a channel with type ECG.
%
% Output:
% D                    - MEEG object with added QRS events(also written on
% disk)
%__________________________________________________________________________
%
% see http://users.fmrib.ox.ac.uk/~rami/fmribplugin/
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_fmrib_qrsdetect.m 3435 2009-10-01 10:24:08Z vladimir $

SVNrev = '$Rev: 3435 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','QRS detect'); spm('Pointer','Watch');

if exist('fmrib_qrsdetect', 'file')~=2
    error('This tool requires FMRIB plugin, see http://users.fmrib.ox.ac.uk/~rami/fmribplugin/');
end

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

D = spm_eeg_load(D);

if ~strcmp(D.type, 'continuous')
    error('Only continuous data is supported at the moment.');
end

% Get the indicex for ECG channel
%--------------------------------------------------------------------------
if ~(isfield(S, 'ecgchan') && ~isempty(S.ecgchan))
   ecgchan = setdiff(D.ecgchannels, D.badchannels);
   if length(ecgchan)~=1
       [selection, ok]= listdlg('ListString', D.chanlabels(ecgchan), 'SelectionMode', 'single' ,'Name', 'Select ecg channel' , 'ListSize', [400 300]);
       if ~ok
           return;
       end
       if ~isempty(ecgchan)
           ecgchan = ecgchan(selection);
       else
           ecgchan = selection;
       end
       S.ecgchan = D.chanlabels(ecgchan);
   end
elseif ~isnumeric(S.ecgchan)
    ecgchan = D.indchannel(S.ecgchan);
else
    ecgchan = S.ecgchan;
end
   
% Detect QRS peaks using FMRIB plugin
%--------------------------------------------------------------------------
EEG = [];
EEG.data = D(ecgchan, :, :);
EEG.srate = D.fsample;
qrs = fmrib_qrsdetect(EEG,1);

if ~isempty(qrs)
    % Update the event structure
    %--------------------------------------------------------------------------
    for n = 1:D.ntrials
        cqrs   = qrs(qrs>(D.nsamples*(n-1)) & qrs<(D.nsamples*n));
        ctime  = D.trialonset(n)+(cqrs - D.nsamples*(n-1)-1)/D.fsample;
        ctime  = num2cell(ctime);
        
        ev = events(D, n);
        
        if iscell(ev)
            ev = ev{1};
        end
        
        Nevents = numel(ev);
        for i=1:numel(ctime)
            ev(Nevents+i).type     = 'artefact';
            ev(Nevents+i).value    = 'qrs';
            ev(Nevents+i).duration = [];
            ev(Nevents+i).time     = ctime{i};
        end
        
        if ~isempty(ev)
            [tevent, I] = sort([ev.time]);
            ev = ev(I);
            D = events(D, n, ev);
        end
    end    
else
    warning(['No QRS events detected in the selected channel ' D.chanlabels(ecgchannel)]);
end

%-Update history (not necessary; leave as call to spm_eeg_montage?) Save new evoked M/EEG dataset
%--------------------------------------------------------------------------
D = D.history(mfilename, S);

save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
spm('FigName','QRS detect: done'); spm('Pointer','Arrow');
