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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_artefact.m 4350 2011-06-13 16:31:02Z vladimir $

SVNrev = '$Rev: 4350 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG artefact detection'); spm('Pointer','Watch');

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

if isequal(D.type, 'continuous')
    error('Artefact detection can only be applied to epoched data');
end

%-Backward compatibility
%--------------------------------------------------------------------------
persistent runonce
if isfield(S, 'artefact')
    if isempty(runonce)
        warning(['The old version of the artefact function will be deprecated in the future']);
        runonce = 1;
    end
    D = spm_eeg_artefact5(S);
    return;
end

if ~isfield(S, 'badchanthresh')
    S.badchanthresh = 0.2;
end

%-Create a copy of the dataset
%--------------------------------------------------------------------------
S1 =[];
S1.D = D;
S1.newname = ['a' D.fname];
S1.updatehistory = 0;
D = spm_eeg_copy(S1);


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
badchanind  = intersect(find(mean(bad, 2)>S.badchanthresh), meegchannels(D));
badchanind  = union(badchanind, D.badchannels);
goodchanind = setdiff(1:D.nchannels, badchanind);

%-Classify trials as bad if they have artefacts in good M/EEG channels
%-or in non-M/EEG channels
%--------------------------------------------------------------------------
badtrialind = find(any(bad(goodchanind, :)));

%-Update and save new dataset
%--------------------------------------------------------------------------
if ~isempty(badtrialind)
    D = reject(D, badtrialind, 1);
end

if ~isempty(badchanind)
    D = badchannels(D, badchanind, ones(size(badchanind)));
end

D = D.history('spm_eeg_artefact', S);
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
