function D = spm_eeg_ft2spm(ftdata, filename)
% Converter from Fieldtrip (http://www.ru.nl/fcdonders/fieldtrip/)
% data structures to SPM8 file format
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_ft2spm.m 2720 2009-02-09 19:50:46Z vladimir $

% If raw format
if iscell(ftdata.time)

    if length(ftdata.time)>1
        % Initial checks
        if any(diff(cellfun('length', ftdata.time))~=0)
            error('SPM can only handle data with equal trial lengths.');
        else
            times=cell2mat(ftdata.time(:));
            if any(diff(times(:, 1))~=0) || any(diff(times(:, end))~=0)
                error('SPM can only handle data the same trial borders.');
            end
        end
    end

    Ntrials=length(ftdata.trial);
    Nchannels  = size(ftdata.trial{1},1);
    Nsamples  = size(ftdata.trial{1},2);

    data = zeros(Nchannels, Nsamples, Ntrials);

    for n=1:Ntrials
        data(:,:,n) = ftdata.trial{n};
    end

    ftdata.time  = ftdata.time{1};
else % timelockanalysis format
    rptind=strmatch('rpt', tokenize(ftdata.dimord, '_'));
    if isempty(rptind)
        rptind=strmatch('subj', tokenize(ftdata.dimord, '_'));
    end

    timeind=strmatch('time', tokenize(ftdata.dimord, '_'));
    chanind=strmatch('chan', tokenize(ftdata.dimord, '_'));

    if ~isempty(rptind)
        if isfield(ftdata, 'trial')
            Ntrials = size(ftdata.trial, rptind);
            data =permute(ftdata.trial, [chanind, timeind, rptind]);
        else
            Ntrials = size(ftdata.individual, rptind);
            data =permute(ftdata.individual, [chanind, timeind, rptind]);
        end
    else
        Ntrials = 1;
        data =permute(ftdata.avg, [chanind, timeind]);
    end

    Nchannels  = size(ftdata.trial, chanind);
    Nsamples  = size(ftdata.trial, timeind);
end

%--------- Start making the header

D = [];

% sampling rate in Hz
D.Fsample = ftdata.fsample;

D.timeOnset = ftdata.time(1);

% Number of time bins in peri-stimulus time
D.Nsamples = Nsamples;

% Names of channels in order of the data
D.channels = struct('label', ftdata.label);

D.trials = repmat(struct('label', {'Undefined'}), 1, Ntrials);

[pathname, fname] = fileparts(filename);

D.path = pathname;
D.fname = [fname '.mat'];

D.data.fnamedat = [fname '.dat'];
D.data.datatype = 'float32-le';

if Ntrials == 1
    datafile = file_array(fullfile(D.path, D.data.fnamedat), [Nchannels Nsamples], D.data.datatype);
    % physically initialise file
    datafile(end,end) = 0;
    datafile(:, :) = data;
else
    datafile = file_array(fullfile(D.path, D.data.fnamedat), [Nchannels Nsamples Ntrials], D.data.datatype);
    % physically initialise file
    datafile(end,end) = 0;

    for i = 1:Ntrials
        datafile(:, :, i) = data(:, :, i);
    end
end

if isfield(ftdata, 'hdr') &&  isfield(ftdata.hdr, 'grad')
    D.sensors.meg = forwinv_convert_units(ftdata.hdr.grad, 'mm');
end

D = meeg(D);

% Set channel types to default
S1 = [];
S1.task = 'defaulttype';
S1.D = D;
S1.updatehistory = 0;
D = spm_eeg_prep(S1);

if Ntrials == 1
    D = type(D, 'continuous');
else
    D = type(D, 'single');
end

save(D);


function [tok] = tokenize(str, sep, rep)

% TOKENIZE cuts a string into pieces, returning a cell array
%
% Use as
%   t = tokenize(str, sep)
%   t = tokenize(str, sep, rep)
% where str is a string and sep is the separator at which you want
% to cut it into pieces.
%
% Using the optional boolean flag rep you can specify whether repeated
% seperator characters should be squeezed together (e.g. multiple
% spaces between two words). The default is rep=1, i.e. repeated
% seperators are treated as one.

% Copyright (C) 2003-2006, Robert Oostenveld


tok = {};
f = find(str==sep);
f = [0, f, length(str)+1];
for i=1:(length(f)-1)
    tok{i} = str((f(i)+1):(f(i+1)-1));
end

if nargin<3 || rep
    % remove empty cells, which occur if the separator is repeated (e.g. multiple spaces)
    tok(find(cellfun('isempty', tok)))=[];
end