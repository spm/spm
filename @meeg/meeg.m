function D = meeg(varargin)
% Function for creating meeg objects.
% FORMAT m = meeg(varargin)
%
% SPM8 format consists of a binary memory - mapped data file and a header mat file. 
% 
% The header file will contain an object of class meeg called D. All information other than 
% data is contained in this object and access to the data is via methods of the object. Also, 
% arbitrary data can be stored inside the object if their field names do not conflict with existing methods' names.
% 
% The following is a description of theinternal implementation of meeg. 
% 
% Fields of meeg:
%  
% .Fsample – sampling rate
% .data – substruct containing the information related to the bimary data file and memory mapping.
% 
%	Subfields of .data
%       .y – reference to the data file_array
%       .fnamedat – name of the .dat file
%       .scale – scaling coefficients for file_array
%       .datatype – type of the data for file_array
% 
% 
% .Nsamples – length of the trial (whole data if the file is continuous). 
% .timeOnset – the peri-stimulus time of the first sample in the trial
% 
% .fname, .path – strings added by spm_eeg_ldata to keep track of where a header struct was loaded from.
% 
% .trials –this describes the segments of the epoched file and is also a structure array. 
% 
%   Subfields of .trials
% 
%       .label – user-specified string for the condition
%       .onset – time of the first sample in seconds in terms of the original file
%       .bad – 0 or 1 flag to allow rejection of trials.
%       .repl – for epochs that are averages – number of replications used for the average.
%       .events – this is a structure array describing events related to each trial. 
% 
%           Subfields of .events
% 
%           .type – string (e.g. 'front panel trigger') 
%           .value - number or string, can be empty (e.g. 'Trig 1').
%           .time – in seconds in terms of the original file
%           .duration - in seconds 
% 
% .channels - This is a structure array which is a field of meeg. 
%             length(channels) should equal size(.data.y, 1) and the order must correspond to the order of channels in the data.
% 
%	Subfields of .channels
% 
%       .label – channel label which is always a string
%       .type –a string, possible values - 'MEG', 'EEG', 'VEOG', 'HEOG', 'EMG' ,'LFP' etc.
%       .units – units of the data in the channel. 
%       .bad – 0 or 1 flag to mark channels as bad.
%       .X_plot2D, .Y_plot2D – positions on 2D plane (formerly in ctf). NaN for channels that should not be plotted.
%       .tra – this is a vector which must be compatible with the .sensors subfield. It defines how the signal in the 
%              corresponding channel is derived from the sensors.
% 
% .sensors  
% 
% 
%	Subfields of .sensors
% 
%   Only sensors with positions should be included (i.e., EEG, MEG)
% 
%       .pnt   -  Mx3 matrix with the position of each sensor
%       .ori   - Mx3 matrix with the orientation of each MEG coil. 
%       .label - labels of the sensors (Mx1 cell array of strings).
%       .type – cell array of strings with sensor types ('EEG', 'MEG').
% 
%   M - number of sensors can be larger or smaller than the number of channels. 
% 
% 
% .fiducials – MRI fiducials for coregistration (3x3 matrix)
% 
% .artifacts - structure array with fields .start and .stop expressed in seconds in original file time.  
% 
% .history – structure array describing commands that modified the file.
% 
%	Subfields of .history:
% 
%       .function – string, the function name 
%       .arguments – cell array, the function arguments
%       .time - when function call was made
% 
% ______________________________________________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id $

if nargin==1
    if isstruct(varargin{1})
        [result D] = checkmeeg(varargin{1}, 'basic');
        D = class(D, 'meeg');
        return;
    elseif isa(varargin{1},'meeg'),
        D = varargin{1};
        return;
    end;
end;

D=[];
D.Nsamples = 0;
[res D] = checkmeeg(D);
D = class(D, 'meeg');
