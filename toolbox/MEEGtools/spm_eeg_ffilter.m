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
%
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
    [B, A] = butter(S.order,S.freq/(S.D.fsample/2),S.band);
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
chunkSamples= round(S.chunkSize/(8*size(S.D,1))*1e6);
begs=1:chunkSamples:size(S.D,2);
ends = (begs+chunkSamples-1);
if(ends(end)>size(S.D,2))
    ends(end)= size(S.D,2);
end

%-Forward direcion
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Filtering Channels (Forward)',spm('time'));
zf = zeros(max(length(A),length(B))-1,length(Fchannels));
for i =1:length(begs)
    %display(['completed chunk ' num2str(i) ' of ' num2str(length(begs))]);
    input = squeeze(S.D(:,begs(i):ends(i),:))';
    output = input;
    [output(:,Fchannels), zf] = filter(B,A,input(:,Fchannels),zf);
    Dnew(:,begs(i):ends(i),:)= output';
end
Dnew.save();

%-backward  direcion
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Filtering Channels (backward)',spm('time'));
switch S.dir
    case 'twopass'
        zf = zeros(max(length(A),length(B))-1,length(Fchannels));
        ends = flip(ends);
        begs = flip(begs);
        for i =1:length(ends)
     %       display(['completed chunk ' num2str(i) ' of ' num2str(length(begs))]);
            input = squeeze(Dnew(:,begs(i):ends(i),:))';
            input = flip(input);
            output = input;
            [output(:,Fchannels), zf] = filter(B,A,input(:,Fchannels),zf);
            Dnew(:,begs(i):ends(i),:)= flip(output)';
        end
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




