function D = spm_eeg_artefact5(S)
% Simple artefact detection, optionally with robust averaging
% FORMAT D = spm_eeg_artefact5(S)
%
% S                     - input structure (optional)
% (optional) fields of S:
%   S.D                 - MEEG object or filename of M/EEG mat-file with
%                         continuous data
%   S.artefact with entries (all optional):
%    External_list      - (0: no/1: yes) flag, whether to flag trials as
%                         artefactual or clean
%    out_list           - vector of artefactual trial indices (if
%                         External_list yes)
%    in_list            - vector of clean trial indices (if
%                         External_list yes)
%    Weighted           - (0: no/1: yes) flag, whether to use robust
%                         averaging
%    Weightingfunction  - parameter used for robust averaging
%    Smoothing          - parameter used for robust averaging
%    Check_Threshold    - (0: no/1: yes) flag, whether to threshold trials
%    channels_threshold - vector of indices which channels to threshold (if
%                         Check_Threshold yes)
%    threshold          - vector of thresholds, with which the absolute
%                         data values are compared against (if Check_Threshold yes). 
%                         Vector length must be either number of selected
%                         channels, or a single number applied to all channels.
%
% Output:
% D                     - MEEG object (also written on disk)
%__________________________________________________________________________
%
% spm_eeg_artefact is an artefact detection routine. The user can specify 
% clean or artefactual trials as vector indices. These trials are not 
% checked by SPM. 
% The function uses simple thresholding to detect artefactual trials and
% bad channels. Optionally, 'robust averaging' can be used to estimate how
% much weight each trial should have in the average to compute the evoked
% response. 
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel, Rik Henson & James Kilner
% $Id: spm_eeg_artefact5.m 3258 2009-07-08 17:46:54Z vladimir $

SVNrev = '$Rev: 3258 $';

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

%-User specified list of artefacted and clean trials
%--------------------------------------------------------------------------
try
    artefact.External_list   = S.artefact.External_list;
catch
    artefact.External_list   = ...
        spm_input('Read own artefact list?','+1','yes|no',[1 0]);
    S.artefact.External_list = artefact.External_list;
end

MustDoWork = 1; % indicate whether user already specified full artefact list

if artefact.External_list
    try
        artefact.out_list    = S.artefact.out_list;
    catch
        artefact.out_list    = ...
            spm_input('List artefactual trials (0 for none)', '+1', 'w', '', Inf);
        S.artefact.out_list  = artefact.out_list;
    end

    if artefact.out_list == 0
        artefact.out_list    = [];
    end

    try
        artefact.in_list     = S.artefact.in_list;
    catch
        artefact.in_list = ...
            spm_input('List clean trials (0 for none)', '+1', 'w', '', Inf);
        S.artefact.in_list   = artefact.in_list;
    end

    if artefact.in_list == 0
        artefact.in_list     = [];
    end

    if any([artefact.out_list; artefact.in_list] < 1 | ...
           [artefact.out_list; artefact.in_list] > D.ntrials)
        error('Trial numbers cannot be smaller than 1 or greater than %d.', D.ntrials);
    end

    % check the lists
    tmp = intersect(artefact.out_list, artefact.in_list);
    if ~isempty(tmp)
        error('These trials were listed as both artefactual and clean: %s', mat2str(tmp));
    end

    % Check whether user has specified all trials
    Iuser = [artefact.out_list; artefact.in_list];
    if length(Iuser) == D.ntrials
        MustDoWork = 0;
    end
end

%-Robust averaging?
%--------------------------------------------------------------------------
try
    artefact.Weighted   = S.artefact.Weighted;
catch
    artefact.Weighted   = spm_input('Robust average?','+1','yes|no',[1 0]);
    S.artefact.Weighted = artefact.Weighted;
end

if artefact.Weighted == 1
    try
        artefact.Weightingfunction   = S.artefact.Weightingfunction;
    catch
        artefact.Weightingfunction   = ...
            spm_input('Offset weighting function by', '+1', 'r', '3', 1);
        S.artefact.Weightingfunction = artefact.Weightingfunction;
    end
    try
        artefact.Smoothing   = S.artefact.Smoothing;
    catch
        artefact.Smoothing   = ...
            spm_input('FWHM for residual smoothing (ms)', '+1', 'r', '20', 1);
        S.artefact.Smoothing = artefact.Smoothing;
    end
    artefact.Smoothing       = round(artefact.Smoothing / 1000 * D.fsample);
end

%-Thresholding details
%--------------------------------------------------------------------------
if MustDoWork
    try
        artefact.Check_Threshold   = S.artefact.Check_Threshold;
    catch
        artefact.Check_Threshold   = ...
            spm_input('Threshold channels?', '+1', 'yes|no', [1 0]);
        S.artefact.Check_Threshold = artefact.Check_Threshold;
    end

    if artefact.Check_Threshold
        try
            artefact.channels_threshold   = S.artefact.channels_threshold;
        catch
            artefact.channels_threshold   = ...
                spm_input('Select channels', '+1', 'i', num2str(sort([D.meegchannels D.eogchannels])));
            S.artefact.channels_threshold = artefact.channels_threshold;
        end

        try
            artefact.threshold     = S.artefact.threshold;
            if length(artefact.threshold) == 1
                artefact.threshold = repmat(artefact.threshold, 1, ...
                                      length(artefact.channels_threshold));
            end
        catch
            str  = 'threshold[s]';
            Ypos = -1;
            while 1
                if Ypos == -1
                    [artefact.threshold, Ypos] = spm_input(str, '+1', 'r', [], [1 Inf]);
                else
                    artefact.threshold         = spm_input(str, Ypos, 'r', [], [1 Inf]);
                end
                if length(artefact.threshold) == 1
                    artefact.threshold = repmat(artefact.threshold, 1, ...
                                      length(artefact.channels_threshold));
                end
                if length(artefact.threshold) == length(artefact.channels_threshold), break, end
                str = sprintf('Enter a scalar or [%d] vector', length(artefact.channels_threshold));
            end
            S.artefact.threshold    = artefact.threshold;
        end
    else
        artefact.channels_threshold = 1: nchannels(D);
        artefact.threshold          = Inf(1, length(artefact.channels_threshold));
    end
end % MustDoWork

%-Artefact detection
%--------------------------------------------------------------------------
if MustDoWork

    % matrix used for detecting bad channels
    Mbad = zeros(D.nchannels, D.ntrials);
    % flag channels that were already marked as bad
    Mbad(D.badchannels, :) = 1;

    % cell vectors of channel-wise indices for thresholded trials
    thresholded = cell(1, length(artefact.channels_threshold));

    Tchannel = artefact.threshold;

    %-First flag bad channels based on thresholding
    %----------------------------------------------------------------------
    spm_progress_bar('Init', D.ntrials, '1st pass - Trials thresholded');
    if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials, 100));
    else Ibar = [1:D.ntrials]; end
    
    for i = 1:D.ntrials

        d = squeeze(D(artefact.channels_threshold, :, i));

        % indices of channels that are above threshold and not marked as bad
        Id = find(max(abs(d')) > Tchannel & ~Mbad(artefact.channels_threshold, i)');
        Mbad(intersect(artefact.channels_threshold(Id), D.meegchannels), i) = 1;

        if ismember(i, Ibar), spm_progress_bar('Set', i); end

    end

    spm_progress_bar('Clear');

    %-Flag channels as bad if 20% of trials above threshold
    %----------------------------------------------------------------------
    ind  = find(mean(Mbad, 2) > 0.2);

    Mbad = zeros(D.nchannels, D.ntrials);
    Mbad(ind, :) = 1;

    %-Report on command line and set badchannels
    %----------------------------------------------------------------------
    if isempty(ind)
        fprintf('There isn''t a bad channel.\n');                       %-#
        D = badchannels(D, [D.meegchannels], zeros(length(D.meegchannels), 1));
    else
        lbl = D.chanlabels(ind);
        if ~iscell(lbl), lbl = {lbl}; end
        disp(['Bad channels: ', sprintf('%s ', lbl{:})]);               %-#
        D = badchannels(D, ind, ones(length(ind), 1));
    end

    %-Weighted averaging
    %----------------------------------------------------------------------
    if artefact.Weighted == 1
        
        cl = condlist(D);

        allWf = zeros(D.nchannels, D.ntrials * D.nsamples);
        tloops = setdiff(artefact.channels_threshold, ind);

        for i = 1:D.nconditions
            
            nbars = D.nconditions * length(tloops);
            spm_progress_bar('Init', nbars, '2nd pass - robust averaging');
            if nbars > 100, Ibar = floor(linspace(1, nbars,100));
            else Ibar = [1:nbars]; end

            trials = pickconditions(D, deblank(cl{i}), 0);

            for j = tloops %loop across electrodes
                if ismember((i-1)*length(tloops)+j, Ibar)
                    spm_progress_bar('Set', (i-1)*length(tloops)+j);
                end
                tempdata = max(abs(squeeze(D(j, :, trials))));
                itrials  = trials;

                itrials(tempdata>Tchannel(j)) = '';
                tdata    = squeeze(D(j, :, itrials));
                [B, bc]  = spm_eeg_robust_averaget(tdata, ...
                           artefact.Weightingfunction, artefact.Smoothing);
                bc  = bc(:);
                ins = 0;

                for n = itrials'
                    ins = ins + 1;
                    allWf(j, (n-1)*D.nsamples+1 : n*D.nsamples) = ...
                                  bc((ins-1)*D.nsamples+1:ins*D.nsamples)';
                end
            end
        end

        spm_progress_bar('Clear');

        artefact.weights = allWf;
        D.artefact = artefact;

    else % if artefact.Weighted == 1

        %-2nd round of thresholding, but excluding bad channels
        %------------------------------------------------------------------
        index = [];

        spm_progress_bar('Init', D.ntrials, '2nd pass - Trials thresholded');
        if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
        else Ibar = [1:D.ntrials]; end

        for i = 1:D.ntrials

            d = squeeze(D(artefact.channels_threshold, :, i));

            % indices of channels that are above threshold
            Id = find(max(abs(d')) > Tchannel & ~Mbad(artefact.channels_threshold, i)');
            Mbad(artefact.channels_threshold(Id), i) = 1;

            if any(Id), index = [index i]; end % reject

            % stow away event indices for which good channels were
            % above threshold
            for j = Id
                thresholded{j} = [thresholded{j} i];
            end
            
            if ismember(i, Ibar), spm_progress_bar('Set', i); end
        
        end

        spm_progress_bar('Clear');
        
        fprintf('%d rejected trials: %s\n', length(index), mat2str(index)); %-#

        if ~isempty(index)
            D = reject(D, index, 1);
        end
    end

    D.thresholded = thresholded;

end % MustDoWork

%-User-specified lists override any artefact classification
%--------------------------------------------------------------------------
if artefact.External_list
    if ~isempty(artefact.out_list)
        D = reject(D, artefact.out_list, 1);
    end
    if ~isempty(artefact.in_list)
        D = reject(D, artefact.in_list, 0);
    end
end

%-Save new dataset
%--------------------------------------------------------------------------
D = D.history('spm_eeg_artefact', S);
copyfile(fullfile(D.path, D.fnamedat), fullfile(D.path, ['a' D.fnamedat]));
D = fnamedat(D, ['a' D.fnamedat]);
D = fname(D, ['a' D.fname]);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG artefact detection: done'); spm('Pointer','Arrow');
