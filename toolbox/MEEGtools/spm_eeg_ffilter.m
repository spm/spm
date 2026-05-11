function Dnew = spm_eeg_ffilter(S)
% Filter M/EEG data (optimised for long datasets)
% FORMAT D = spm_eeg_filter(S)
%
% S           - input structure
%  Fields of S:
%   S.D       - MEEG object or filename of M/EEG mat-file
%
%   S.band    - filterband [low|high|bandpass|stop]
%   S.freq    - cutoff frequency(-ies) [Hz]
%
%  Optional fields:
%   S.type    - filter type [default: 'butterworth']
%                 'butterworth': Butterworth IIR filter
%                 'fir':         FIR filter (using MATLAB fir1 function)
%   S.order   - filter order [default: 5 for Butterworth]
%   S.dir     - filter direction [default: 'twopass']
%                 'onepass':         forward filter only
%                 'onepass-reverse': reverse filter only, i.e. backward in time
%                 'twopass':         zero-phase forward and reverse filter
%   S.prefix  - prefix for the output file [default: 'f']
%   S.chunkSize - size data segment (MB) to be filtered (default 200)
% Currently only set up for single epoch data
% D           - MEEG object (also written to disk)
%__________________________________________________________________________

% Tim Tierney
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if ~isfield(S, 'dir'),    S.dir    = 'twopass';     end
if ~isfield(S, 'chunkSize'),    S.chunkSize    = 200;     end
if ~isfield(S, 'prefix'), S.prefix = 'f';           end
if ~isfield(S, 'order'),  S.order=5; end

if ~isfield(S, 'band')
    error('A frequency band must be supplied')
end

%-Get MEEG object
%--------------------------------------------------------------------------
D = spm_eeg_load(S.D);

Ntrials=D.ntrials,

%-Check band
%--------------------------------------------------------------------------
switch lower(S.band)

    case {'low','high'}
        if numel(S.freq)~=1
            error('Cutoff frequency should be a single number.');
        end

        if S.freq < 0 || S.freq > D.fsample/2
            error('Cutoff must be > 0 & < half sample rate.');
        end

    case {'bandpass','stop'}
        if S.freq(1) < 0 || S.freq(2) > D.fsample/2 || S.freq(1) > S.freq(2)
            error('Incorrect frequency band specification.');
        end

    otherwise
        error('Incorrect filter band.')
end

fprintf('%-40s: %30s\n',...
    ['Filter ' S.band ' (' 'butterworth' ', ' S.dir ')'],...
    ['[' num2str(S.freq) '] Hz']);                                      %-#

%-Filter
%==========================================================================

%-Generate new meeg object with new filenames
Dnew = copy(D, [S.prefix fname(D)]);

%-Determine channels for filtering
Fchannels = D.indchantype('Filtered');

if isempty(Fchannels)
    warning('No channels suitable for filterning found. Please check your channel type specification.');
end


%- stability check
%--------------------------------------------------------------------------

unstable =1;
while (unstable)
    try
        [B, A] = butter(S.order,S.freq/(D.fsample/2),S.band);
    catch
        [B, A] = spm_biquad(S.order,S.freq,D.fsample,S.band);
    end
    unstable = any(abs(roots(A))>=1);
    if(unstable)
        ft_warning('instability detected - reducing the %dth order filter to an %dth order filter', S.order, S.order-1);
        S.order=S.order-1;
        if(S.order<1)
            error('Cannot stabilise filter');
        end
    end
end


%- Work out memory chunk size
%--------------------------------------------------------------------------
chunkSamples= round(S.chunkSize/(8*size(D,1))*1e6); % in MB
%chunkSamples= round(S.chunkSize) % /(8*size(D,1))*1e6)
begs=1:chunkSamples:size(D,2);
ends = (begs+chunkSamples-1);
if(ends(end)>size(D,2))
    ends(end)= size(D,2);
end

%-Forward direcion
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Filtering Channels (Forward)',spm('time'));
switch S.dir
    case {'twopass','onepass'}
        for tr=1:Ntrials,
            zf = zeros(max(length(A),length(B))-1,length(Fchannels));
            for i =1:length(begs)
                %display(['completed chunk ' num2str(i) ' of ' num2str(length(begs)) 'trial' num2str(tr)]);
                input = squeeze(D(:,begs(i):ends(i),tr))';
                output = input;
                [output(:,Fchannels), zf] = filter(B,A,input(:,Fchannels),zf);
                Dnew(:,begs(i):ends(i),tr)= output';
            end
        end;% for trials
        Dnew.save();
    otherwise
        fprintf('No forward filtering')
end; % switch

%-backward  direcion
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Filtering Channels (backward)',spm('time'));
switch S.dir
    case {'twopass','onepass-reverse'}
        zf = zeros(max(length(A),length(B))-1,length(Fchannels));
        ends = flip(ends);
        begs = flip(begs);
        for tr=1:Ntrials,
            for i =1:length(ends)
                
                input = squeeze(Dnew(:,begs(i):ends(i),tr))';
                input = flip(input);
                output = input;
                [output(:,Fchannels), zf] = filter(B,A,input(:,Fchannels),zf);
                Dnew(:,begs(i):ends(i),tr)= flip(output)';
            end
        end; % for tr

    case 'onepass'

    otherwise
        error('supplied direction not supported');
end


Dnew.save();
%
%-Cleanup
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Completed',spm('time'));

end

function [B,A]= spm_biquad(n, f,fs, band)
% Biquad filter coefficients
% FORMAT [B,A] = spm_biquad(n, f,fs, band)
%   n       - order of filter
%   Wn      - frequencies to filter[]
%   band    - filterband [low|high|bandpass|stop]
%
% [B,A]     - IIR filter coefficients
%__________________________________________________________________________

% Tim Tierney
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

Wn= f/fs;

fc = prod(Wn)^(1/length(Wn));
K = tan(pi*fc);
ncascades = ceil(n/2);
b0=[];
b1=[];
b2=[];
a1=[];
a2=[];

switch band
    case 'low'
        for i = 0:(ncascades-1)
            Q= 1/(2*cos(pi/(ncascades*4)+i*pi/(ncascades*2)));
            norm = 1/(1+K/Q +K*K);
            b0 =[b0 K*K*norm];
            b1 =[b1 2*b0(i+1)];
            b2 =[b2 b0(i+1)];
            a1 =[a1 2*(K*K-1)*norm];
            a2 =[a2 (1-K/Q+K*K)*norm];
        end
    case 'high'
        for i = 0:(ncascades-1)
            Q= 1/(2*cos(pi/(ncascades*4)+i*pi/(ncascades*2)));
            norm = 1 / (1 + K / Q + K * K);
            b0 =[b0 1*norm];
            b1 =[b1 -2*b0(i+1)];
            b2 = [b2 b0(i+1)];
            a1 = [a1 2*(K*K-1)*norm];
            a2 = [a2 (1- K/Q+K*K)*norm];
        end
    case 'bandpass'
        K = tan(pi*Wn(1));
        for i = 0:(ncascades-1)
            Q= 1/(2*cos(pi/(ncascades*4)+i*pi/(ncascades*2)));
            norm = 1 / (1 + K / Q + K * K);
            b0 =[b0 1*norm];
            b1 =[b1 -2*b0(i+1)];
            b2 = [b2 b0(i+1)];
            a1 = [a1 2*(K*K-1)*norm];
            a2 = [a2 (1- K/Q+K*K)*norm];
        end
        K = tan(pi*Wn(2));
        for i = 0:(ncascades-1)
            Q= 1/(2*cos(pi/(ncascades*4)+(i)*pi/(ncascades*2)));
            norm = 1/(1+K/Q +K*K);
            b0 =[b0 K*K*norm];
            b1 =[b1 2*b0(i+1+ncascades)];
            b2 =[b2 b0(i+1+ncascades)];
            a1 =[a1 2*(K*K-1)*norm];
            a2 =[a2 (1-K/Q+K*K)*norm];
        end
    case 'stop'
        ncascades = ncascades*2;
        Qeff = sqrt(prod(f))/(f(2)-f(1));
        Q = Qeff*sqrt(2^(ncascades)-1);
        for i = 0:(ncascades-1)
            norm = 1 / (1 + K / Q + K * K);
            b0 =[b0 (1+K*K)*norm];
            b1 =[b1 2*(K*K-1)*norm];
            b2 = [b2 b0(i+1)];
            a1 = [a1 b1(i+1)];
            a2 = [a2 (1- K/Q+K*K)*norm];
        end
end
b = [b0' b1' b2'];
a = [ones(length(b0),1) a1' a2'];

B= b(1,:);
for i = 1:(length(b0)-1)
    B= conv(B,b(i+1,:));
end

A= a(1,:);
for i = 1:(length(b0)-1)
    A= conv(A,a(i+1,:));
end


end



