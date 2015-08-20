function D = meeg(varargin)
% Function for creating meeg objects.
% FORMAT 
%         D = meeg;
%             returns an empty object
%         D = meeg(D);
%             converts a D struct to object or does nothing if already
%             object
%         D = meeg(nchannels, nsamples, ntrials)
%             return a time dataset with default settings
%         D = meeg(nchannels, nfrequencies, nsamples, ntrials)
%             return TF time dataset with default settings        
%         
% SPM MEEG format consists of a header object optionally linked to 
% binary data file. The object is usually saved in the header mat file
%
% The header file will contain a struct called D. All 
% information other than data is contained in this struct and access to the
% data is via methods of the object. Also, arbitrary data can be stored 
% inside the object if their field names do not conflict with existing 
% methods' names.
%
% The following is a description of the internal implementation of meeg.
%
% Fields of meeg:
% .type - type of data in the file: 'continuous', 'single', 'evoked'
% .Fsample - sampling rate
%
% .data - file_array object linking to the data or empty if unlinked
%
%
% .Nsamples - length of the trial (whole data if the file is continuous).
% .timeOnset - the peri-stimulus time of the first sample in the trial (in sec)
%
% .fname, .path - strings added by spm_eeg_load to keep track of where a 
%                 header struct was loaded from.
%
% .trials - this describes the segments of the epoched file and is also a 
%           structure array.
%
%   Subfields of .trials
%
%       .label - user-specified string for the condition
%       .onset - time of the first sample in seconds in terms of the 
%                original file
%       .bad - 0 or 1 flag to allow rejection of trials.
%       .repl - for epochs that are averages - number of replications used 
%               for the average.
%       .tag  - the user can put any data here that will be attached to
%               the respective trial. This is useful e.g. to make sure the 
%               relation between regressors and data is not broken when
%               removing bad trials or merging files.
%       .events - this is a structure array describing events related to 
%                 each trial.
%
%           Subfields of .events
%
%           .type - string (e.g. 'front panel trigger')
%           .value - number or string, can be empty (e.g. 'Trig 1').
%           .time - in seconds in terms of the original file
%           .duration - in seconds
%
% .channels - This is a structure array which is a field of meeg.
%             length(channels) should equal size(.data.y, 1) and the order 
%             must correspond to the order of channels in the data.
%
%   Subfields of .channels
%
%       .label - channel label which is always a string
%       .type - a string, possible values - 'MEG', 'EEG', 'VEOG', 'HEOG', 
%               'EMG' ,'LFP' etc.
%       .units - units of the data in the channel.
%       .bad - 0 or 1 flag to mark channels as bad.
%       .X_plot2D, .Y_plot2D - positions on 2D plane (formerly in ctf). NaN
%                              for channels that should not be plotted.
%
% .sensors
%
%
%   Subfields of .sensors (optional)
%       .meg - struct with sensor positions for MEG (subfields: .pnt .ori .tra .label)
%       .eeg - struct with sensor positions for MEG (subfields: .pnt .tra .label)
%
% .fiducials - headshape and fiducials for coregistration with sMRI
%      
%   Subfiels of .fiducials (optional)
%       .pnt - headshape points
%       .fid.pnt - fiducial points
%       .fid.label - fiducial labels
%
% .transform - additional information for transfomed (most commonly time-frequency) data
%    Subfields of .transform 
%        .ID - 'time', 'TF', or 'TFphase'
%        .frequencies (optional)
%
% .history - structure array describing commands that modified the file.
%
%   Subfields of .history:
%
%       .function - string, the function name
%       .arguments - cell array, the function arguments
%       .time - when function call was made
%
% .other - structure used to store other information bits, not fitting the
%          object structure at the moment,
%       for example:
%       .inv - structure array corresponding to the forw/inv problem in MEEG.
%       .val - index of the 'inv' solution currently used.
%
% .condlist - cell array of unique condition labels defining the proper
%        condition order
%
% .montage - structure used to store info on on-line montage used
%       .M contains transformation matrix of the montage and names of 
%           original and new channels (+ new channels definition)
%       .Mind indicates which montage to use
%__________________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: meeg.m 6525 2015-08-20 10:03:16Z vladimir $

switch nargin
    case 0
        D    = struct('Nsamples', 0);
    case 1
        D = varargin{1};
    case 2
        error('Illegal number of arguments');
    case 3
        D    = struct('Nsamples', varargin{2}, ...
            'channels', struct('bad', num2cell(zeros(1, varargin{1}))), ...
            'trials', struct('bad', num2cell(zeros(1, varargin{3}))));
    case 4
        D    = struct('Nsamples', varargin{3}, ...
            'channels', struct('bad', num2cell(zeros(1, varargin{1}))), ...
            'trials', struct('bad', num2cell(zeros(1, varargin{4}))));
        D.transform.ID = 'TF';
        D.transform.frequencies = 1:varargin{2};
end

if ~isa(D, 'meeg')
    D = class(checkmeeg(D), 'meeg');
end



