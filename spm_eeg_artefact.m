function D = spm_eeg_artefact(S)
% Simple artefact detection, optionally with robust averaging
% FORMAT D = spm_eeg_artefact(S)
%
% S                 - input structure
%
% fields of S:
%   S.mode            'reject' [default]: reject bad channels and trials
%                     'mark': scan the data and create events marking the
%                             artefacts
%   S.D             - MEEG object or filename of M/EEG mat-file
%   S.badchanthresh - fraction of trials (or time) with artefacts above
%                     which a channel is declared as bad [default: 0.2]
%
%   S.append        - 1 [default]: append new markings to existing ones
%                     0: overwrite existing markings
%   S.methods       - structure array with configuration parameters for
%                     artefact detection plugins
%   S.prefix        - prefix for the output file [default: 'a']
%
% Output:
% D                 - MEEG object (also written on disk)
%__________________________________________________________________________
%
% This is a modular function for which plugins can be developed to detect
% artefacts with any algorithm.
% The name of a plugin function should start with 'spm_eeg_artefact_'.
% Several plugins are already implemented annd they can be used as
% templates for new plugins:
%
% peak2peak         - thresholds peak-to-peak amplitude
% (spm_eeg_artefact_peak2peak)
%
% jump              - thresholds the difference between adjacent samples
% (spm_eeg_artefact_jump)
%
% flat              - detects flat segments in the data
% (spm_eeg_artefact_flat)
%__________________________________________________________________________
% Copyright (C) 2008-2017 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_artefact.m 7132 2017-07-10 16:22:58Z guillaume $

SVNrev = '$Rev: 7132 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG artefact detection'); spm('Pointer','Watch');

if ~isfield(S, 'mode'),            S.mode = 'reject';          end
if ~isfield(S, 'badchanthresh'),   S.badchanthresh = 0.2;      end
if ~isfield(S, 'append'),          S.append = true;            end
if ~isfield(S, 'prefix'),          S.prefix        = 'a';      end

%-Get MEEG object
%--------------------------------------------------------------------------
D = spm_eeg_load(S.D);

if isequal(S.mode, 'reject')
    
    if isequal(D.type, 'continuous')
        error('Artefact rejection can only be applied to epoched data.');
    end
    
    %-Create a copy of the dataset
    %----------------------------------------------------------------------
    D = copy(D, [S.prefix D.fname]);
    
    %-Run the artefact detection routines
    %----------------------------------------------------------------------
    bad = zeros(D.nchannels, D.ntrials);
    
    for i = 1:numel(S.methods)
        chanind = D.selectchannels(S.methods(i).channels);
        if S.append
            chanind = setdiff(chanind, D.badchannels);
        end
        
        if ~isempty(chanind)
            S1 =  S.methods(i).settings;
            
            if isempty(S1)
                S1 = [];
            end
            
            S1.mode = S.mode;
            S1.D    = D;            
            S1.chanind = chanind;
            
            bad = bad | feval(['spm_eeg_artefact_' S.methods(i).fun], S1);
        end
    end
    
    %-Classify MEEG channels as bad if the fraction of bad trials exceeds threshold
    %------------------------------------------------------------------------------
    badchanind  = intersect(find(mean(bad, 2)>S.badchanthresh), indchantype(D, 'MEEG'));
    badchanind  = union(badchanind, D.badchannels);
    goodchanind = setdiff(1:D.nchannels, badchanind);
    
    %-Classify trials as bad if they have artefacts in good M/EEG channels
    %-or in non-M/EEG channels
    %----------------------------------------------------------------------
    badtrialind = find(any(bad(goodchanind, :), 1));
    
    %-Update and save new dataset
    %----------------------------------------------------------------------
    if ~S.append
        D = badtrials(D, ':', 0);
        D = badchannels(D, ':', 0);
    end
    
    D = badtrials(D, badtrialind, 1);
    D = badchannels(D, badchanind, ones(size(badchanind)));
    
    %-Report on command line
    %----------------------------------------------------------------------
    fprintf('%d rejected trials: %s\n', length(badtrialind), num2str(badtrialind));
    
elseif isequal(S.mode, 'mark')
    
    %-Create a copy of the dataset
    %----------------------------------------------------------------------
    D = copy(D, [S.prefix D.fname]);
           
    %-Run the artefact detection routines
    %----------------------------------------------------------------------
    for i = 1:numel(S.methods)
        chanind = setdiff(D.selectchannels(S.methods(i).channels), D.badchannels);
       
        if ~isempty(chanind)
            S1 =  S.methods(i).settings;
            
            if isempty(S1)
                S1 = [];
            end
            
            S1.mode = S.mode;
            S1.D    = D;
            S1.badchanthresh = S.badchanthresh;
            S1.append =  S.append;
            S1.chanind = chanind;
            
            
            D =  feval(['spm_eeg_artefact_' S.methods(i).fun], S1);
        end
    end    
    
    meegind = D.indchantype({'MEEG', 'LFP'});
    
    bad = squeeze(mean(mean(badsamples(D, meegind, ':', ':'), 2), 3)) > S.badchanthresh;
    
    badchanind = meegind(bad);
    
    %-Update and save new dataset
    %----------------------------------------------------------------------
    D = badchannels(D, badchanind, 1);
else
    error('Invalid mode specification.');
end

if isempty(badchanind)
    fprintf('There are no bad channels.\n');
else
    lbl = D.chanlabels(badchanind);
    if ~iscell(lbl), lbl = {lbl}; end
    fprintf('%d bad channels: %s\n', numel(lbl), sprintf('%s ', lbl{:}));
end

D = D.history(mfilename, S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
spm('FigName','M/EEG artefact detection: done'); spm('Pointer','Arrow');
