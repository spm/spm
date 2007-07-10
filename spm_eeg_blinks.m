function D = spm_eeg_blinks(S);
% simple linear regression method for blink correction using VEOG 
% FORMAT D = spm_eeg_blinks(S)
% 
% S       - struct (optional)
% (optional) fields of S:
% D	  - filename of EEG-data file or EEG data struct
% blinks  - struct with the following fields:
%		veog_peak - minimum amplitude for a blink (eg 200uV)
%		veog_width - maximum width for a blink (eg 200ms)
% 		blink_sample - number of sample points for blink (eg 50)
%		filter_para - if want to smooth before detecting blinks
%		blink_wgts - if have already calculated weights
%		blink_method - blink shape or blink area
%
% D	   - EEG data struct (also written to files)
%_______________________________________________________________________
% 
% spm_eeg_blinks	simple method for blink correction
%
% Three basic stages:	1. Detect blinks (using peak,width,filter)
%			2. Calculate correction weights
%			3. Subtract weighted VEOG from EEG channels
%______________________________________________________________________
%
% Created Rik Henson, 2003
% Updated for SPM5 by Doris Eckstein, 2006
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

[Finter,Fgraph,CmdLine] = spm('FnUIsetup', 'EEG filter setup',0);

try
    D = S.D;
catch
    D = spm_select(1, '\.mat$', 'Select EEG mat file');
end

P = spm_str_manip(D, 'H');

try
    D = spm_eeg_ldata(D);
catch    
    error(sprintf('Trouble reading file %s', D));
end

if ~isfield(S,'blinks')
    S.blinks=[];
end

% indices of channels without VEOG
in_no_veog = setdiff([1:D.Nchannels], D.channels.veog);


if isfield(D.channels, 'blink_wgts')
    disp('Found these blink weights...')
    for c = 1:length(in_no_veog)
        disp(sprintf('%s\t%3.2f', ...
            D.channels.name{in_no_veog(c)},D.channels.blink_wgts(c)))
    end
    try 
	Detect_Blinks = S.blinks.Detect_Blinks;
    catch
    	Detect_Blinks = 0;
    end
else
    Detect_Blinks = 1;
end


if Detect_Blinks == 1

 Pos=1;
 if ~isfield(S.blinks, 'veog_peak')
    S.blinks.veog_peak =...
        spm_input('Min VEOG peak for blink (e.g, uV)', Pos, 'i', '200', 1);
    Pos = Pos+1;
 end
 D.blinks.veog_peak = S.blinks.veog_peak;

 if ~isfield(S.blinks, 'veog_width')
    S.blinks.veog_width = ...
        spm_input('Max VEOG width for blink (ms)', Pos, 'i', '200', 1);
    Pos = Pos+1;
 end
 D.blinks.veog_width = S.blinks.veog_width;

 if ~isfield(S.blinks, 'blink_sample')
    suggest = S.blinks.veog_width*D.Radc/(2*1000);
    S.blinks.blink_sample = ...
        spm_input('blink sampling points', Pos, 'i', suggest, 1);
    Pos = Pos+1;
 end
 D.blinks.blink_sample = S.blinks.blink_sample;

 if ~isfield(S.blinks, 'filter_para')
    S.blinks.filter_para =...
        spm_input('Lowpass filter cutoff (0 none)', Pos, 'r', '20', 1);
    Pos = Pos+1;
 end

 if S.blinks.filter_para > 0,
    D.blinks.filter = round( ((0.31*D.Radc / S.blinks.filter_para) - 1 )/2 );
else
    D.blinks.filter = 0;
 end

 if ~isfield(S.blinks, 'blink_method')
    Ctype = {'Area','Shape'};
    S.blinks.blink_method =...
        spm_input('correction method','+1','m',Ctype);
    Pos = Pos+1;
 end
 D.blinks.blink_method = S.blinks.blink_method;

end

try 
    Correct_Blinks = S.blinks.Correct_Blinks;
catch
    Correct_Blinks = 0;
end

try 				% EXCLUDE an epoch (eg for continuous data)
    Epoch = S.blinks.epoch;	% in ms, Eg [-100 600]
    Epoch = round(Epoch*D.Radc/1000);
catch
    Epoch = [];
end



% If wanted to check for artefacts before detecting blinks....
%
%Pos = 1;
%if ~isfield(S.blinks, 'deblock')
%	S.blinks.deblock = ...
%	 	spm_input_ui('Deblock min sampling points', Pos, 'i', '', 1);
%	Pos = Pos + 1;
%end
%
%if ~isfield(S.blinks, 'saturation')
%	S.blinks.saturation = ...
%	 	spm_input_ui('Saturation in bits (eg 12 bit)', Pos, 'i', '', 1);
%	Pos = Pos+1;
%end
%
% In digitisation points plus safety buffer
%sat_threshold = 2^(S.blinks.saturation-1) - S.blinks.saturation;
%
% Needs information about base gain, which assumed in D.channels.calib
%S.blinks.saturation = sat_threshold*D.channels.calib;
%
% Take minimum of saturation and user-specified blinks
%thr_channel = min([S.blinks.saturation; S.blinks.threshold]);


spm('Pointer', 'Watch');


%------------------------
% 1. Veog blink detection
%------------------------

if Detect_Blinks

Iblink=[];

Dveog = squeeze(D.data(D.channels.veog,:,:));

if(size(D.data,3)==1) 		% only one epoch
    Dveog = Dveog(:);
end

if D.blinks.filter > 0			% '2-step' digital filter kernel
    LF = fix(D.blinks.filter);
    F = ones(1, 2*D.blinks.filter+1)./(2*D.blinks.filter+1);
    F = conv(F,F);

    for i = 1:D.Nevents
        tmp = conv(Dveog(:,i), F');
        Dveog(:,i) = tmp([1:D.Nsamples] + floor(length(F)/2));
    end
end

blink_max_width = D.blinks.veog_width * D.Radc/1000;	% in sample points

if(isfield(D.events,'reject'))
    valid_events = find(~D.events.reject);
else
    valid_events = 1:D.Nevents;
end

D.events.blinks = zeros(1,D.Nevents);

blink_waves=[]; blink_samples=[]; pst_samples=[]; pst_times=[];
index=[]; rindex=[]; rmblink=[];

for n = 1:length(valid_events)
    i = valid_events(n);

    Iblink = find(Dveog(:,i) > D.blinks.veog_peak/2);	% possible blinks

    if(~isempty(Iblink))

        dI = find(diff(Iblink)>1);
        Sblink = [Iblink(1); Iblink(dI+1)];			% possible onsets
        Eblink = [Iblink(dI); Iblink(end)];			% possible durations

	if(~isempty(Epoch))
          for b = 1:length(Sblink)				% cannot vectorise?
	    rels = Sblink(b)-D.events.time;
	    rele = Eblink(b)-D.events.time;
	    if(max(rels(find(rels<0))) > Epoch(1) | min(rels(find(rels>0))) < Epoch(2) | max(rele(find(rele<0))) > Epoch(1))
		rmblink = [rmblink b];
            end
          end
  	end
	Sblink(rmblink)=[];
	Eblink(rmblink)=[];

        Dblink = Eblink - Sblink;				% possible durations

        for b = 1:length(Sblink)				% cannot vectorise?

            blink = Dveog((Sblink(b) : (Sblink(b)+Dblink(b)-1) ),i)';

            [Mblink,Pblink] = max(blink);

            if (Mblink > D.blinks.veog_peak)

                if (length(find(blink > Mblink/2)) < blink_max_width)	% max FWHM

                    bsample = [(Pblink(1)-D.blinks.blink_sample) : ...
                               (Pblink(1)+D.blinks.blink_sample)];

                    if(bsample(1)>0 & bsample(end)<length(blink))

                        tsample = bsample + Sblink(b);

                        if(D.Nevents==1)	% non-epoched
                            pst = find(D.events.time-tsample(1)>0);
                            if(~isempty(pst))
                                pst_samples = [pst_samples D.events.time(pst(1))-tsample(1)];
			        pst_times = [pst_times 1000*pst_samples(end)/D.Radc];
                            end
                        else
                            pst_samples = [pst_samples tsample(1)];
			    pst_times = [pst_times (1000*tsample(1)-D.events.start-1)/D.Radc];
                        end

                        blink_waves = [blink_waves; blink(bsample)];
                        blink_samples = [blink_samples; tsample];
                        index = [index i];
                    end
                else			% reject: VEOG exceeds peak, but too long
                    rindex = [rindex i];
                end
            end
        end
    end
end
disp(sprintf('Valid Blinks (N=%d): %s',length(index),mat2str(index)))
disp(sprintf('Invalid "Blinks" (N=%d): %s',length(rindex),mat2str(rindex)))

D.events.blinks = index;
%D.events.reject(unique(rindex))=1;


%-----------------------
% 2. Calculate weights
%------------------------

% Keep HEOG for correlation purposes (do not necessarily correct)

    if isempty(index)
      warning('No blinks detected!') 
      return;
    end

    VEOGwave = mean(blink_waves)';

    y = zeros(length(in_no_veog), 2*D.blinks.blink_sample+1);
    for n = 1:length(index)
        d = squeeze(D.data(in_no_veog, :, index(n)));
        y = y + d(:, blink_samples(n,:));
    end
    y = y'/length(index);

    switch D.blinks.blink_method

        case 1				% area (equivalent to mean)
	    X = VEOGwave - VEOGwave(1);
	    y = y - kron(ones(size(y,1),1),y(1,:));

%            b = mean(y)/mean(VEOGwave);
%            cc = b;				

        case 2				% shape (linear regression)
            X = detrend(VEOGwave,0);
            y = detrend(y,0);
    end

            pX = pinv(X);
            b  = pX*y;
            Y  = X*b;
            cc = sign(b').*diag(y'*Y)./sqrt(diag(y'*y).*diag(Y'*Y));


    [Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG blink detection');
    spm('Clear',Finter, Fgraph);
    subplot(1,3,1)
    plot(y)
    title('Blink shapes')
    ylabel('Amplitude')
    xlabel('Sample Window')
    subplot(1,3,2)
    plot(b(1,:)',cc,'x')
    title('Weight-Correlation')
    ylabel('Correlation')
    xlabel('Weight')
    subplot(1,3,3)
    hist(pst_times)
    title('PST Histogram')
    ylabel('No. Blinks')
    xlabel('PST')
    drawnow

    disp('      Cor       Wgt')
    for c = 1:length(in_no_veog)
        disp(sprintf('%s\t%3.2f\t%3.2f', ...
            D.channels.name{in_no_veog(c)},cc(c),b(1,c)))
    end
    D.channels.blink_wgts = b(1,:)';
end


%-------------------------
% 3. Veog blink correction
%-------------------------

if(Correct_Blinks)

    disp('Using these weights...')
    for c = 1:length(in_no_veog)
        disp(sprintf('%s\t%3.2f', ...
            D.channels.name{in_no_veog(c)},D.channels.blink_wgts(c)))
    end

    D.fnamedat = ['b' spm_str_manip(D.fnamedat, 't')];

    fpd = fopen(D.fnamedat, 'w');
    d = zeros(D.Nchannels, D.Nsamples);
    D.scale = zeros(D.Nchannels, 1, D.Nevents);

    if(size(D.data,3)==1) 					% only one epoch

        d = squeeze(D.data(:, :,1));

        d(in_no_veog, :) = d(in_no_veog, :) - ...
            D.channels.blink_wgts * d(D.channels.veog,:); %

        D.scale = spm_eeg_write(fpd, d, 2, 'int16');

    else
        for i = 1:D.Nevents

            d = squeeze(D.data(:, :, i));

            d(in_no_veog, :) = d(in_no_veog, :) - ...
                D.channels.blink_wgts * d(D.channels.veog,:);

            D.scale(:,1,i) = spm_eeg_write(fpd, d, 2, 'int16');
        end
    end

    fclose(fpd);

    D.data = [];
    D.fname = ['b' spm_str_manip(D.fname, 't')];

else
    disp('No correction performed')
end

save(D.fname, 'D');

spm('Pointer', 'Arrow');





