function D = spm_eeg_artefact(S)
% Simple artefact detection, optionally with robust averaging
% FORMAT D = spm_eeg_artefact(S)
%
% S                     - input structure
%
% fields of S:
%   S.D                 - MEEG object or filename of M/EEG mat-file with
%   S.badchanthresh     - fraction of trials with artefacts above which a 
%                         channel is declared as bad (default: 0.2)
%
%   S.methods           - a structure array with configuration parameters
%                         for artefact detection plugins.
%   S.prefix     - prefix for the output file (default - 'a')
% Output:
% D                     - MEEG object (also written on disk)
%__________________________________________________________________________
% This is a modular function for which plugins can be developed to detect
% artefacts with any algorithm. There are 3 basic plugins presently
% implemented and they can be used as templates for new plugins.
% The name of a plugin function should start with 'spm_eeg_artefact_'
%
% peak2peak (spm_eeg_artefact_peak2peak) - thresholds peak-to-peak
%                                          amplitude
%
% jump (spm_eeg_artefact_jump)           - thresholds the difference
%                                          between adjacent samples.
%
% flat (spm_eeg_artefact_flat)           - detects flat segments in the
%                                          data
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_artefact.m 5076 2012-11-23 16:05:21Z vladimir $

SVNrev = '$Rev: 5076 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG artefact detection'); spm('Pointer','Watch');

if ~isfield(S, 'badchanthresh'),   S.badchanthresh = 0.2;      end
if ~isfield(S, 'prefix'),          S.prefix        = 'a';      end

%-Get MEEG object
%--------------------------------------------------------------------------
D = spm_eeg_load(S.D);

if isequal(D.type, 'continuous')
    error('Artefact detection can only be applied to epoched data');
end

%-Create a copy of the dataset
%--------------------------------------------------------------------------
D = copy(D, [S.prefix D.fname]);

%-Run the artefact detection routines
%--------------------------------------------------------------------------
bad = zeros(D.nchannels, D.ntrials);

for i = 1:numel(S.methods)
    chanind = setdiff(D.selectchannels(S.methods(i).channels), D.badchannels);
    
    if ~isempty(chanind)
        S1 =  S.methods(i).settings;

        if isempty(S1)
            S1 = [];
        end
        
        S1.D = D;
        S1.chanind = chanind;

        bad = bad | feval(['spm_eeg_artefact_' S.methods(i).fun], S1);
    end
end

%-Classify MEEG channels as bad if the fraction of bad trials exceeds threshold
%-------------------------------------------------------------------------------
badchanind  = intersect(find(mean(bad, 2)>S.badchanthresh), indchantype(D, 'MEEG'));
badchanind  = union(badchanind, D.badchannels);
goodchanind = setdiff(1:D.nchannels, badchanind);

%-Classify trials as bad if they have artefacts in good M/EEG channels
%-or in non-M/EEG channels
%--------------------------------------------------------------------------
badtrialind = find(any(bad(goodchanind, :)));

%-Update and save new dataset
%--------------------------------------------------------------------------
D = badtrials(D, badtrialind, 1);
D = badchannels(D, badchanind, ones(size(badchanind)));


D = D.history(mfilename, S);
save(D);

%-Report on command line
%--------------------------------------------------------------------------
if isempty(badchanind)
    fprintf('There isn''t a bad channel.\n');                                   
else
    lbl = D.chanlabels(badchanind);
    if ~iscell(lbl), lbl = {lbl}; end
    fprintf('%d bad channels: %s\n', numel(lbl), sprintf('%s ', lbl{:}));      
end
fprintf('%d rejected trials: %s\n', length(badtrialind), num2str(badtrialind));

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG artefact detection: done'); spm('Pointer','Arrow');
