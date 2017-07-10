function res = spm_eeg_artefact_heartbeat(S)
% Detects heart beats in SPM continuous data file
% S            - input structure
% fields of S:
%    S.D       - M/EEG object
%    S.chanind - vector of indices of channels that this plugin will look at
%
%    Additional parameters can be defined specific for each plugin.
%
% Output:
% res -
%    If no input is provided the plugin returns a cfg branch for itself.
%
%    If input is provided the plugin returns a matrix of size D.nchannels x D.ntrials
%    with zeros for clean channel/trials and ones for artefacts.
%
% See http://fsl.fmrib.ox.ac.uk/eeglab/fmribplugin/
%__________________________________________________________________________
% Copyright (C) 2008-2017 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_artefact_heartbeat.m 7132 2017-07-10 16:22:58Z guillaume $


%-This part if for creating a config branch that plugs into spm_cfg_eeg_artefact
% Any parameters can be specified and they are then passed to the plugin
% when it's called.
%--------------------------------------------------------------------------
if nargin == 0           
    excwin         = cfg_entry;
    excwin.tag     = 'excwin';
    excwin.name    = 'Excision window';
    excwin.strtype = 'r';
    excwin.num     = [1 1];
    excwin.val     = {0};
    excwin.help    = {'Window (in ms) to mark as bad around each heart beat, 0 to not mark data as bad.'};
    
    heartbeat      = cfg_branch;
    heartbeat.tag  = 'heartbeat';
    heartbeat.name = 'Heart beats';
    heartbeat.val  = {excwin};
    heartbeat.help = {''};
    
    res = heartbeat;
    
    return
end

SVNrev = '$Rev: 7132 $';

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','M/EEG heartbeat detection');

if exist('fmrib_qrsdetect', 'file')~=2
    error('This tool requires FMRIB plugin, see http://fsl.fmrib.ox.ac.uk/eeglab/fmribplugin/');
end

if isequal(S.mode, 'reject')
    error('Only mark mode is supported by this plug-in, use event-based rejection to reject.');
end

D = spm_eeg_load(S.D);

chanind = S.chanind;

if length(chanind)~=1
    error('More than one channel - not currently supported.')
end

% Detect QRS peaks using FMRIB plugin
%--------------------------------------------------------------------------
EEG = [];
EEG.data =  reshape(squeeze(D(chanind,:,:)), 1, []);
EEG.srate = D.fsample;
spikes = fmrib_qrsdetect(EEG,1);

% Update the event structure
%----------------------------------------------------------------------
if ~isempty(spikes)  
    for n = 1:D.ntrials
        cspikes   = spikes(spikes>(D.nsamples*(n-1)) & spikes<(D.nsamples*n));
        ctime  = D.trialonset(n)+(cspikes - D.nsamples*(n-1)-1)/D.fsample;
        ctime  = num2cell(ctime);
        
        ev = events(D, n);
        
        if iscell(ev)
            ev = ev{1};
        end
        
        if ~isempty(ev) && ~S.append
            ind1 = strmatch('artefact_heartbeat', {ev.type}, 'exact');
            if ~isempty(ind1)
                ind2 = strmatch(D.chanlabels(chanind), {ev(ind1).value}, 'exact');
                if ~isempty(ind2)
                    ev(ind1(ind2)) = [];
                end
            end
        end        
        
        Nevents = numel(ev);
        for i=1:numel(ctime)
            if ctime{i} == 0
                %likely to be trial border falsely detected as heartbeat
                continue;
            end
            ev(Nevents+i).type  = 'artefact_heartbeat';
            ev(Nevents+i).value = char(D.chanlabels(chanind));
            if S.excwin == 0
                ev(Nevents+i).duration = [];
                ev(Nevents+i).time     = ctime{i};
            else
                ev(Nevents+i).time     = max(D.trialonset(n), ctime{i} - 5e-4*S.excwin);
                ev(Nevents+i).duration = min(1e-3*S.excwin, (D.time(end)-D.time(1))-(ev(Nevents+i).time-D.trialonset(n)))+...
                    min(ctime{i} - 5e-4*S.excwin, 0);               
            end                      
        end
        
        if ~isempty(ev)
            [dum, I] = sort([ev.time]);
            ev = ev(I);
            D = events(D, n, ev);
        end
    end    
else
    warning('No heartbeat events detected in the selected channel.');
end

res = D;

spm('FigName','M/EEG heartbeat detection: done');
