function D = spm_eeg_artefact(S)
% performs simple artefact detection and correction on epoched EEG data
% FORMAT D = spm_eeg_artefact(S)
%
% S		            - optional struct with the following fields:
% thresholds        - struct with the following (optional) fields:
%   deblock         - number of time points
%   saturation      - 12 or 16 (bits)
%   Check_Threshold - switch (0/1) whether to use thresholding
%   threshold       - vector of channelwise thresholds [\mu Volt]
%   Check_Drift     - switch (0/1) whether to apply drift detection
%   drift           - (Nchannels x 2)-matrix of drift thresholds
%   Correct blinks  - switch (0/1) whether to correct for blinks
%   veog_peaks      - threshold needed for blink correction [\mu Volt]
%   veog_width      - threshold needed for blink correction [time points]
%   blink_sample    - number of time points around a detected blink to
%                       which correction is applied
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
% spm_eeg_artefact implements a simple artefact detection. Optionally, one
% can attempt to use a blink correction. Note that SPM will not remove any
% data, but just indicate, whether it is an artefact or not. The 0/1-vector
% D.events.reject will store that information.
% Alternatively, one can also provide these information as an out_list and
% in_list, see above.
%_______________________________________________________________________
% @(#)spm_eeg_artefact.m	1.1 Stefan Kiebel, Rik Henson 04/06/28

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG artefact setup',0);

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

MustWork = 1;
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
        MustWork = 0;
    end
end

% Does user want to correct for blinks?
spm('Clear', Finter);
try
    Correct_Blinks = S.thresholds.Correct_Blinks;
catch
    Correct_Blinks = spm_input('Correct Blinks?','+1','yes|no',[1 0]);
end

if Correct_Blinks
    try
        D.thresholds.veog_peak = S.thresholds.veog_peak;
    catch
        D.thresholds.veog_peak =...
            spm_input('peak threshold for Veog (uV)', '+1', 'i', '', 1);
    end
    
    try
        D.thresholds.veog_width = S.thresholds.veog_width;
    catch
        
        D.thresholds.veog_width = ...
            spm_input('width threshold for Veog (ms)', '+1', 'i', '', 1);
    end
    
    try
        D.thresholds.blink_sample = S.thresholds.blink_sample;
    catch
        
        D.thresholds.blink_sample = ...
            spm_input('blink sampling points', '+1', 'i', '', 1);
    end
end

spm('Clear',Finter, Fgraph);

if MustWork
    try
        D.thresholds.deblock = S.thresholds.deblock;
    catch
        D.thresholds.deblock = ...
            spm_input('Deblock min sampling points', '+1', 'i', '', 1);
    end
    
    if D.thresholds.deblock == 0
        Check_deblock = 0;
    else
        Check_deblock = 1;
    end
    
    try
        D.thresholds.saturation = S.thresholds.saturation;
    catch
        D.thresholds.saturation = ...
            spm_input('Saturation in bits (eg 12 bit)', '+1', 'i', '', 1);
    end
    
        if D.thresholds.saturation == 0
        Check_saturation = 0;
    else
        Check_saturation = 1;
    end

    
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
    
    try
        Check_Drift = S.thresholds.Check_Drift;
    catch
        Check_Drift = spm_input('Check for drift?','+1','yes|no',[1 0]);
    end
    
    if Check_Drift
        try
            D.thresholds.drift = S.thresholds.drift;
        catch
            for c = 1:D.Nchannels
                D.thresholds.drift(:,c) = spm_input(sprintf('%s: drift value and corrcoeff', D.channels.name{c}), '+1', 'r', '54.9 0.8', 2);
            end
        end
    else
        D.thresholds.drift = kron(ones(1, D.Nchannels), [Inf 0]');
    end
end % MustWork

spm('Pointer', 'Watch');

if MustWork
    
    if Check_saturation == 1
        % In digitisation points plus safety buffer
        Tsat = 2^(D.thresholds.saturation-1) - D.thresholds.saturation;
        
        % Needs information about base gain, which assumed in D.channels.calib
        D.thresholds.saturation = Tsat*D.channels.calib;
        
        % Take minimum of saturation and user-specified thresholds
        Tchannel = min([D.thresholds.saturation; D.thresholds.threshold]);
    else
        Tchannel = D.thresholds.threshold;
    end
    
    D.events.reject = zeros(1, D.Nevents);
    
    % cell vectors of channel-wise indices for thresholded events
    D.channels.thresholded = cell(1, D.Nchannels);
    D.channels.blocked = cell(1, D.Nchannels);
    
    index = [];
    for i = 1:D.Nevents
        
        d = squeeze(D.data(:, :, i));
        
        % indices of channels that are above one of the threshold(s)
        Id = find(max(abs(d')) > Tchannel);
        
        % user-specified thresholding and saturation
        if any(Id)
            index = [index i];
            for j = Id
                D.channels.thresholded{j} = [D.channels.thresholded{j} i];
            end
        end
        
        if Check_deblock
            d = diff(d, 1, 2)==0;
            d = conv2([1], ones(1, D.thresholds.deblock-1), d);
            
            Id = find(any(d >= D.thresholds.deblock-1, 2))';
            if any(Id)
                index = [index i];
                for j = Id
                    D.channels.blocked{j} = [D.channels.blocked{j} i];
                end 
            end
        end
    end
    
    disp(sprintf('Blocked, above threshold or saturated trials (N=%d): %s', length(index), mat2str(index)))
    
    D.events.reject(index) = 1;
    valid_events = find(~D.events.reject);
    
    %---------------------------------------------------------
    % Drift detection
    % Rik's bit
    %---------------------------------------------------------
    index = [];
    if Check_Drift
        
        D.channels.drift = cell(1, D.Nchannels);
        
        % use a linear model 
        X  = detrend(1:D.Nsamples,0)'/D.Nsamples;
        pX = pinv(X);
        
        for i = valid_events
            
            d = detrend(squeeze(D.data(:, :, i))', 'constant');
            
            b = pX * d;
            
            % first threshold: Slope > threshold?
            Ic = find(abs(b) > D.thresholds.drift(1,:));
            
            chan_drift = [];
            for c = Ic
                cc = corrcoef(X, d(:,c));
                % second threshold: Correlation coefficient > threshold?
                if abs(cc(1,2)) > D.thresholds.drift(2,c)
                    chan_drift = [chan_drift c];
                end
            end
            
            if ~isempty(chan_drift)
                index = [index i];
                for c = chan_drift
                    D.channels.drift{c} = [D.channels.drift{c} i];
                end
            end		
        end
        
        disp(sprintf('Drift trials (N=%d): %s',length(index),mat2str(index)))
        
        D.events.reject(index) = 1;
        valid_events=find(~D.events.reject);
    end
end % MustWork

%%%%%%%%%%%%%%%%%%%%%%%
% MODIFIED BY RIK (to handle multiple blinks per epoch, for continuous data?)
%---------------------
% Veog blink detection
%---------------------
%
% COULD INSERT OPTION HERE TO SMOOTH BEFORE DETECTING BLINKS

index = [];
Iblink = [];

if Correct_Blinks
    
    Dveog = squeeze(D.data(D.channels.veog, :, :));
    
    % if there is only one epoch
    if size(Dveog,1) == 1,  Dveog = Dveog(:); end
    
    D.events.blinks = zeros(1,D.Nevents);
    
    blink_waves=[]; index=[]; mpoint=[];
    blink_max_width = D.Radc * D.thresholds.veog_width / 1000;	% in sample points
    
    for n = 1:length(valid_events)
        i = valid_events(n);
        
        Iblink = find(Dveog(:,i) > D.thresholds.veog_peak/2);	% possible blinks
        
        if(~isempty(Iblink))
            
            Sblink = [0 find(diff(Iblink)>1)' length(Iblink)];	% blink onsets
            Dblink = diff(Sblink);					% blink durations
            
            Fblink = find(Dblink > (2 * D.thresholds.blink_sample));
            
            for b = Fblink
                
                blink = Dveog(Iblink((Sblink(b)+1):Sblink(b+1)),i)';
                
                [Mblink,Pblink] = max(blink);			% ASSUMES ONE PEAK
                
                if (Mblink > D.thresholds.veog_peak & ...
                        length(find(blink > Mblink/2)) < blink_max_width)	% FWHM
                    
                    sample = [(Pblink(1)-D.thresholds.blink_sample) : ...
                            (Pblink(1)+D.thresholds.blink_sample)];	% takes first max 
                    
                    if(sample(1)>0 & sample(end)<length(blink))
                        blink_waves = [blink_waves; blink(sample)];
                        index = [index i];
                        mpoint = [mpoint Pblink(1)+Iblink(Sblink(b)+1)-1];
                    end
                end
            end
        end
    end
    disp(sprintf('Blink trials: %s',mat2str(index)))
    
    D.events.blinks(unique(index))=1;
    
    %---------------------
    % Veog blink analysis
    %---------------------
    % !!!Cannot replicate ERPSLAB for some reason!!!
    
    % indices of channels without EOGs
    in_no_veog = setdiff([1:D.Nchannels], D.channels.veog);
    
    Ywave = zeros(length(in_no_veog), 2*D.thresholds.blink_sample+1);
    
    if(~isempty(index))
        for n = 1:length(index)		
            Ywave = Ywave + squeeze(D.data(in_no_veog, ...
                mpoint(n) - D.thresholds.blink_sample : ...
                mpoint(n) + D.thresholds.blink_sample, index(n)));
        end
        Ywave = Ywave'/length(index);
        
        VEOGwave = mean(blink_waves)';
        
        X = detrend(VEOGwave,0);
        Y = detrend(Ywave,0);
        
        pX = pinv(X);
        b = pX*Y;
        
        cc = (X(:,1)'*Y)./sqrt(X(:,1)'*X(:,1)*diag(Y'*Y)');
        figure
        plot(Ywave)
        disp('      Cor       Wgt')
        disp([cc' b(1,:)'])
    end
    
end	% End correct blinks flag

% User-specified lists override any artefact classification
if D.thresholds.External_list
    D.events.reject(D.thresholds.out_list) = 1;
    D.events.reject(D.thresholds.in_list)  = 0;
end

% Save the data
D.fnamedat = ['a' D.fnamedat];

fpd = fopen(fullfile(P, D.fnamedat), 'w');

d = zeros(D.Nchannels, D.Nsamples);
k = 1;
D.scale.dim = [1 3];
D.scale.values = zeros(D.Nchannels, D.Nevents);

for i = 1:D.Nevents
    d = squeeze(D.data(:, :, i));
    
    % Veog blink correction
    if Correct_Blinks & ~isempty(Iblink)
        d(in_no_veog, :) = d(in_no_veog, :) - b(1,:)' * Dveog(:,i)';
    end
    
    D.scale.values(:, i) = max(abs(d'))./32767;
    d = int16(d./repmat(D.scale.values(:, i), 1, D.Nsamples));
    fwrite(fpd, d, 'int16');
end
fclose(fpd);

D.datatype = 'int16';

D.data = [];
D.fname = ['a' D.fname];
if str2num(version('-release'))>=14
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end
spm('Pointer', 'Arrow');
