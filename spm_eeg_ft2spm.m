function D = spm_eeg_ft2spm(ftdata, filename)
% Converter from Fieldtrip (http://www.ru.nl/fcdonders/fieldtrip/)
% data structures to SPM file format
%_______________________________________________________________________
% Copyright (C) 2008-2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_ft2spm.m 5438 2013-04-24 10:38:47Z vladimir $

isTF = 0;

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
else
    Nchannels  = numel(ftdata.label);
    Nsamples   = length(ftdata.time);
    
    rptind=strmatch('rpt', tokenize(ftdata.dimord, '_'));
    if isempty(rptind)
        rptind=strmatch('subj', tokenize(ftdata.dimord, '_'));
    end

    timeind=strmatch('time', tokenize(ftdata.dimord, '_'));
    chanind=strmatch('chan', tokenize(ftdata.dimord, '_'));

    if any(ismember({'trial', 'individual', 'avg'}, fieldnames(ftdata) )) % timelockanalysis
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
    elseif isfield(ftdata, 'powspctrm')
        isTF = 1;
        Nfrequencies = numel(ftdata.freq);
        freqind = strmatch('freq', tokenize(ftdata.dimord, '_'));
        if ~isempty(rptind)
            Ntrials = size(ftdata.powspctrm, rptind);
            data = permute(ftdata.powspctrm, [chanind, freqind, timeind, rptind]);
        else
            Ntrials = 1;
            data = permute(ftdata.powspctrm, [chanind, freqind, timeind]);
        end
    end
end

%--------- Start making the header

D = [];

% sampling rate in Hz
if isfield(ftdata, 'fsample')
    D.Fsample = ftdata.fsample;
else
    D.Fsample = 1./mean(diff(ftdata.time));
end

D.timeOnset = ftdata.time(1);

% Number of time bins in peri-stimulus time
D.Nsamples = Nsamples;

% Names of channels in order of the data
D.channels = struct('label', ftdata.label);

D.trials = repmat(struct('label', {'Undefined'}), 1, Ntrials);

[pathname, fname] = fileparts(filename);

D.path = pathname;
D.fname = [fname '.mat'];

fnamedat = [fname '.dat'];

if ~isTF
    if Ntrials == 1
        datafile = file_array(fullfile(D.path, fnamedat), [Nchannels Nsamples], 'float32-le');
        datafile(:, :) = data;
    else
        datafile = file_array(fullfile(D.path, fnamedat), [Nchannels Nsamples Ntrials], 'float32-le');
        datafile(:, :, :) = data;
    end
else
    if Ntrials == 1
        datafile = file_array(fullfile(D.path, fnamedat), [Nchannels Nfrequencies Nsamples], 'float32-le');
        datafile(:, :, :) = data;
    else
        datafile = file_array(fullfile(D.path, fnamedat), [Nchannels Nfrequencies Nsamples Ntrials], 'float32-le');
        datafile(:, :, :, :) = data;
    end    
    D.transform.ID = 'TF';
    D.transform.frequencies = ftdata.freq;
end

D.data = datafile;

D = meeg(D);

if  isfield(ftdata, 'hdr')
    % Uses fileio function to get the information about channel types stored in
    % the original header. This is now mainly useful for Neuromag support but might
    % have other functions in the future.
    origchantypes = ft_chantype(ftdata.hdr);
    [sel1, sel2] = spm_match_str(D.chanlabels, ftdata.hdr.label);
    origchantypes = origchantypes(sel2);
    if length(strmatch('unknown', origchantypes, 'exact')) ~= numel(origchantypes)
        D.origchantypes = struct([]);
        D.origchantypes(1).label = ftdata.hdr.label(sel2);
        D.origchantypes(1).type = origchantypes;
    end
end

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

if  isfield(ftdata, 'hdr') && isfield(ftdata.hdr, 'grad')
    D = sensors(D, 'MEG', ft_convert_units(ftdata.hdr.grad, 'mm'));
    
    S = [];
    S.task = 'project3D';
    S.modality = 'MEG';
    S.updatehistory = 0;
    S.D = D;

    D = spm_eeg_prep(S);
end

D = check(D);

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