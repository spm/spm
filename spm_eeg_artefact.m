function D = spm_eeg_artefact(S)
% performs simple artefact detection and correction on epoched EEG data
% FORMAT D = spm_eeg_artefact(S)
%
% S		            - optional struct with the following fields:
% thresholds        - struct with the following (optional) fields:
%   Check_Threshold - switch (0/1) whether to use thresholding
%   threshold       - vector of channelwise thresholds [\mu Volt]
%   External_list   - switch (0/1) whether user has information about
%                       clean or artefactual trials
%   out_list        - indices of trials that are artefactual (SPM will not
%                     examine these trials for artefacts)
%   in_list         - indices of trials that are clean (SPM will not
%                     examine these trials for artefacts)
%
% Output:
% D			- EEG data struct (also written to files)
%_______________________________________________________________________
%
% spm_eeg_artefact implements a simple artefact detection. Note that SPM
% will not remove any data, but just flag, whether it is an artefact or
% not for subsequent processing stages. This information is stored in the
% 0/1-vector D.events.reject will. Alternatively, one can also provide
% flags (e.g. from another program) in an out_list and in_list, see above.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel, Rik Henson & James Kilner
% $Id$


[Finter,Fgraph,CmdLine] = spm('FnUIsetup', 'EEG artefact setup',0);

try
    D = S.D;
catch
    D = spm_get(1, '.mat', 'Select EEG mat file');
end

P = spm_str_manip(D, 'H');

try
    D = spm_eeg_ldata(D);
catch
    error(sprintf('Trouble reading file %s', D));
end

if ~isfield(D, 'thresholds')
    D.thresholds = [];
end

try
    D.thresholds.External_list = S.thresholds.External_list;
catch
    D.thresholds.External_list = spm_input('Read own artefact list?','+1','yes|no',[1 0]);
end

MustDoWork = 1; % flag to indicate whether user already specified full artefact list

if D.thresholds.External_list
    try
        D.thresholds.out_list = S.thresholds.out_list;
    catch
        D.thresholds.out_list = ...
            spm_input('List artefactual trials (0 for none)', '+1', 'w', '', inf);
    end

    if D.thresholds.out_list == 0
        D.thresholds.out_list = [];
    end

    try
        D.thresholds.in_list = S.thresholds.in_list;
    catch
        D.thresholds.in_list = ...
            spm_input('List clean trials (0 for none)', '+1', 'w', '', inf);
    end

    if D.thresholds.in_list == 0
        D.thresholds.in_list = [];
    end

    if any([D.thresholds.out_list; D.thresholds.in_list] < 1 | [D.thresholds.out_list; D.thresholds.in_list] > D.Nevents)
        error('Trial numbers cannot be smaller than 1 or greater than %d.', D.Nevents);
    end

    % check the lists
    tmp = intersect(D.thresholds.out_list, D.thresholds.in_list);
    if ~isempty(tmp)
        error('These trials were listed as both artefactual and clean: %s', mat2str(tmp));
    end

    % Check whether user has specified all trials
    Iuser = [D.thresholds.out_list D.thresholds.in_list];
    if length(Iuser) == D.Nevents
        MustDoWork = 0;
    end
end

try
    Weighted = S.artefact.weighted;
catch
    Weighted = spm_input('robust average?','+1','yes|no',[1 0]);
end

if Weighted ==1
    try
        wtrials = S.artefact.wtrials;
    catch
        wtrials = spm_input('robust across trials?','+1','yes|no',[1 0]);
    end
end

spm('Clear',Finter, Fgraph);

[Finter,Fgraph,CmdLine] = spm('FnUIsetup', 'EEG artefact setup',0);

if MustDoWork
    try
        Check_Threshold = S.thresholds.Check_Threshold;
    catch
        Check_Threshold = spm_input('Threshold channels?','+1','yes|no',[1 0]);
    end

    if Check_Threshold
        try
            D.thresholds.threshold = S.thresholds.threshold;
        catch
            str = 'threshold[s]';
            Ypos = -1;

            while 1
                if Ypos == -1
                    [D.thresholds.threshold, Ypos] = spm_input(str, '+1', 'r', [], [1 Inf]);
                else
                    D.thresholds.threshold = spm_input(str, Ypos, 'r', [], [1 Inf]);
                end
                if length(D.thresholds.threshold) == 1
                    D.thresholds.threshold = D.thresholds.threshold * ones(1, D.Nchannels);
                end

                if length(D.thresholds.threshold) == D.Nchannels, break, end
                str = sprintf('enter a scalar or [%d] vector', D.Nchannels);
            end
        end
    else
        D.thresholds.threshold = kron(ones(1, D.Nchannels), Inf);
    end

end % MustDoWork

spm('Pointer', 'Watch'); drawnow

% matrix used for detecting bad channels
Mbad = zeros(D.Nchannels, D.Nevents);
% flag channels that were already marked as bad
if isfield(D.channels, 'Bad')
    Mbad(D.channels.Bad, :) = 1;
end

D.events.reject = zeros(1, D.Nevents);

% cell vectors of channel-wise indices for thresholded events
D.channels.thresholded = cell(1, D.Nchannels);

if MustDoWork

    Tchannel = D.thresholds.threshold;

    spm_progress_bar('Init', D.Nevents, '1st pass - Events thresholded'); drawnow;
    if D.Nevents > 100, Ibar = floor(linspace(1, D.Nevents,100));
    else, Ibar = [1:D.Nevents]; end

    % first flag bad channels based on thresholding
    for i = 1:D.Nevents

        d = squeeze(D.data(:, :, i));

        % indices of channels that are above threshold and not marked as
        % bad
        Id = find(max(abs(d')) > Tchannel & ~Mbad(:, i)');
        Mbad(Id, i) = 1;

        if ismember(i, Ibar)
            spm_progress_bar('Set', i);
            drawnow;
        end

    end

    spm_progress_bar('Clear');

    % flag channels as bad if 20% of events above threshold
    s = sum(Mbad, 2)/D.Nevents;
    ind = find(s > 0.2);

    % remove EOGs from bad channels
    ind = setdiff(ind, [D.channels.heog D.channels.veog]);

    Mbad = zeros(D.Nchannels, D.Nevents);
    Mbad(ind, :) = 1;

    % report on command line
    if isempty(ind)
        disp(sprintf('There isn''t a bad channel.'));
    else
        disp(['Bad channels: ', sprintf('%s ', D.channels.name{ind})])
        D.channels.Bad = ind;

    end

    if Weighted == 1
        % weighted averaging

        allWf=zeros(D.events.Ntypes, D.Nchannels, D.Nevents*D.Nsamples);
        tloops = [1:D.Nchannels];
        tloops(ind) = [];

        for i = 1:D.events.Ntypes

            trials = find(D.events.code == D.events.types(i));

            for j = tloops %loop across electrodes
                tdata = squeeze(D.data(j, :, trials));

                if wtrials == 1
                    [B, bc] = spm_eeg_robust_average(tdata);
                else
                    [B, bc] = spm_eeg_robust_averaget(tdata);
                end

                allWf(i, j, [1:length(trials)]*D.Nsamples) = bc';
            end
        end

        D.weights = allWf;

    else

        % 2nd round of thresholding, but excluding bad channels
        index = [];

        spm_progress_bar('Init', D.Nevents, '2nd pass - Events thresholded'); drawnow;
        if D.Nevents > 100, Ibar = floor(linspace(1, D.Nevents,100));
        else, Ibar = [1:D.Nevents]; end

        for i = 1:D.Nevents

            d = squeeze(D.data(:, :, i));

            % indices of channels that are above threshold
            Id = find(max(abs(d')) > Tchannel & ~Mbad(:, i)');
            Mbad(Id, i) = 1;

            if any(Id)
                % reject
                index = [index i];

                % stow away event indices for which good channels were
                % above threshold
                for j = Id
                    D.channels.thresholded{j} = [D.channels.thresholded{j} i];
                end
            end
            if ismember(i, Ibar)
                spm_progress_bar('Set', i);
                drawnow;
            end

        end
        spm_progress_bar('Clear');

    end

    disp(sprintf('%d rejected trials: %s', length(index), mat2str(index)))

    D.events.reject(index) = 1;

end % MustDoWork

% User-specified lists override any artefact classification
if D.thresholds.External_list
    D.events.reject(D.thresholds.out_list) = 1;
    D.events.reject(D.thresholds.in_list)  = 0;
end

% Save the data
copyfile(fullfile(D.path, D.fnamedat), fullfile(D.path, ['a' D.fnamedat]));
D.fnamedat = ['a' D.fnamedat];

spm_progress_bar('Clear');

D.data = [];
D.fname = ['a' D.fname];

if str2num(version('-release'))>=14 
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end

spm('Pointer', 'Arrow');
