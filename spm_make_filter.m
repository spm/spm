function [K] = spm_make_filter(k,RT,filterHF,filterLF,norm)
% [K] = spm_make_filter(k,RT,filterHF,filterLF,norm)
%
% General routine to construct band pass filter convolution matrix;
% First, a high pass filter is designed (specified in filterLF),
% then a low pass (filterHF). These are combined in temporal domain.
% See spm_fmri_spm_ui for filter structure.
% norm = 'norm' for normalisation such that the signal variance
% is unchanged by the filter.
%
% k 		: number of scans
% RT		: Repetition time in seconds
% filterHF	: cut High Frequencies
% filterLF	: cut Low Frequencies
% norm 		: normalization


%------------------------------------------------------------------
% Low frequencies %
%------------------------------------------------------------------
%cLFmenu   = {'Cut LF with defaults',...
%	     'None',...
%	     'Specify'};
%

sLF=char(deblank(filterLF.Menu(filterLF.Choice)));

switch sLF

	case 'None'
		KLF=eye(k);
	
	case {'Cut LF with defaults','Specify'}
		freqLF=2*RT/filterLF.Param;
		NLF=16;	
		filterLF_t=fir1(NLF,freqLF,'high');
		%- Warning : fir1 gives appropriate results only when 
		%- the cutting frequency is not too close to the 0 
		%- frequency (in number of frequency samples)
		KLF=[eye(NLF) zeros(NLF,k-NLF);...
		toeplitz([filterLF_t(NLF+1) zeros(1,k-NLF-1)],...
		[fliplr(filterLF_t) zeros(1,k-NLF-1)])];

	otherwise
		warning('choice pb !');

end
	
%------------------------------------------------------------------
% High frequencies %
%------------------------------------------------------------------
%cHFmenu   = {'Cut HF with default hrf',...
%	     'Cut HF with hrf derivative',...
%	     'None',...
%	     'Specify with gaussian filter',...
%	     'Specify with cutting frequency'};


sHF=char(deblank(filterHF.Menu(filterHF.Choice)));

switch sHF

	case 'None'
		KHF=eye(k);

	case 'Specify with cutting frequency'
		freqHF=2*RT/filterHF.Param;
		NHF=16;	
		filterHF_t=fir1(NHF,freqHF);
		%- Warning : fir1 gives appropriate results only when 
		%- the cutting frequency is not too close to the sampling 
		%- frequency (in number of frequency samples)
		KHF=[eye(NHF) zeros(NHF,k-NHF);...
		toeplitz([filterHF_t(NHF+1) zeros(1,k-NHF-1)],...
		[fliplr(filterHF_t) zeros(1,k-NHF-1)])];


	case 'Cut HF with default hrf'  	
		[h,p]=spm_hrf(RT);p(7)=RT*k;[hrf,p]=spm_hrf(RT,p);
		filterHF_t=[zeros(k-1,1);hrf];
		KHF=toeplitz(filterHF_t(k:2*k-1),filterHF_t(k:-1:1));

	case 'Cut HF with hrf derivative'

		[h,p]=spm_hrf(RT);p(7)=RT*k;[hrf,p]=spm_hrf(RT,p);
		filterHF_t=[zeros(k-1,1);gradient(hrf)];
		KHF=toeplitz(filterHF_t(k:2*k-1),filterHF_t(k:-1:1));
	
	case 'Specify with gaussian filter'
		t=0:k-1;
		filterHF_t=exp(-t.^2/(2*(filterHF.Param)^2));
		KHF=toeplitz(filterHF_t);

	otherwise
		warning('choice pb !');

end


%------------------------------------------------------------------
% Filter %
%------------------------------------------------------------------

K= KLF*KHF;

if nargin == 5

	switch norm

	case 'norm' 
		D=sum(K'.^2);
		K=inv(diag(sqrt(D)))*K;
	otherwise 
		error('unknown normalisation');
		
	end % switch norm
end

% function fir1 taken from matlab signal toolbox
% Hamming function has been copied in.

function [b,a] = fir1(N,Wn,varargin)
%FIR1   FIR filter design using the window method.
%   B = FIR1(N,Wn) designs an N'th order lowpass FIR digital filter
%   and returns the filter coefficients in length N+1 vector B.
%   The cut-off frequency Wn must be between 0 < Wn < 1.0, with 1.0 
%   corresponding to half the sample rate.  The filter B is real and
%   has linear phase, i.e., even symmetric coefficients obeying B(k) =
%   B(N+2-k), k = 1,2,...,N+1.
%
%   If Wn is a two-element vector, Wn = [W1 W2], FIR1 returns an
%   order N bandpass filter with passband  W1 < W < W2.
%   B = FIR1(N,Wn,'high') designs a highpass filter.
%   B = FIR1(N,Wn,'stop') is a bandstop filter if Wn = [W1 W2].
%
%   If Wn is a multi-element vector, 
%          Wn = [W1 W2 W3 W4 W5 ... WN],
%   FIR1 returns an order N multiband filter with bands
%    0 < W < W1, W1 < W < W2, ..., WN < W < 1.
%   B = FIR1(N,Wn,'DC-1') makes the first band a passband.
%   B = FIR1(N,Wn,'DC-0') makes the first band a stopband.
%
%   For filters with a passband near Fs/2, e.g., highpass
%   and bandstop filters, N must be even.
%   
%   By default FIR1 uses a Hamming window.  Other available windows,
%   including Boxcar, Hanning, Bartlett, Blackman, Kaiser and Chebwin
%   can be specified with an optional trailing argument.  For example,
%   B = FIR1(N,Wn,kaiser(N+1,4)) uses a Kaiser window with beta=4.
%   B = FIR1(N,Wn,'high',chebwin(N+1,R)) uses a Chebyshev window.
%
%   By default, the filter is scaled so the center of the first pass band 
%   has magnitude exactly one after windowing. Use a trailing 'noscale' 
%   argument to prevent this scaling, e.g. B = FIR1(N,Wn,'noscale'), 
%   B = FIR1(N,Wn,'high','noscale'), B = FIR1(N,Wn,wind,'noscale').
%
%   See also KAISERORD, FIRCLS1, FIR2, FIRLS, FIRCLS, CREMEZ,
%            REMEZ, FREQZ, FILTER.

%   FIR1 is an M-file implementation of program 5.2 in the IEEE
%   Programs for Digital Signal Processing tape. 

%   Author(s): L. Shure
%              L. Shure, 4-5-90, revised
%              T. Krauss, 3-5-96, revised
%   Copyright (c) 1988-98 by The MathWorks, Inc.
%   $Revision: 2.1 $  $Date: 1999-03-05 16:28:27 $

%   Reference(s):
%     [1] "Programs for Digital Signal Processing", IEEE Press
%         John Wiley & Sons, 1979, pg. 5.2-1.

nargchk(2,5,nargin);

% Up to 3 optional input arguments, always in this order:
%   1 - Filter type flag, can be 'high', 'stop', '', 'DC-0', or 'DC-1'
%   2 - Window vector
%   3 - 'noscale' flag

% default optional parameter values:
Ftype = '';
Wind = [];
Scale = 'scale';

switch length(varargin)
case 1
    if isstr(varargin{1})&(length(varargin{1})>0)
        s = upper(varargin{1});
        switch upper(s)
        case {'SCALE','NOSCALE'}
            Scale = s;
        otherwise
            Ftype = s;
        end
    else
        Wind = varargin{1};
    end
case 2
    if isstr(varargin{1})
        Ftype = varargin{1};
    else
        Wind = varargin{1};
    end
    if isstr(varargin{2})
        Scale = varargin{2};
    else
        Wind = varargin{2};
    end
case 3
    Ftype = varargin{1};
    Wind = varargin{2};
    Scale = varargin{3};
end

switch upper(Scale)
case 'NOSCALE'
    SCALING = 0;
case 'SCALE'
    SCALING = 1;
otherwise
    error('Scaling option must be ''noscale'' or ''scale''.')
end

if isempty(N) | ~isnumeric(N) | ~isreal(N) | N~=round(N) | N<=0
    error('N must be a real, positive integer.')
end

Ftype = upper(Ftype);
if ~strncmp(Ftype,'HIGH',1) & ~strncmp(Ftype,'STOP',1) & ...
   ~strncmp(Ftype,'DC-0',4) & ~strncmp(Ftype,'DC-1',4) & ...
   ~isempty(Ftype)
    error('Filter type must be ''high'',''stop'',''DC-0'', or ''DC-1''.')
end

nw = length(Wind);

nbands = length(Wn) + 1;
if (nbands > 2) & isempty(Ftype)
    Ftype = 'DC-0';  % make sure default 3 band filter is bandpass
end
First_Band = isempty(findstr('DC-0',Ftype)) & isempty(findstr('HIGH',Ftype));
mags = rem( First_Band + (0:nbands-1), 2);

L = N + 1;
odd = rem(L, 2);
if (mags(nbands) & ~odd)
      disp('For highpass and bandstop filters, order must be even.')
      disp('Order is being increased by 1.')
      N = N + 1;  L = L + 1;
      odd = 1;
end
if nw ~= 0 & nw ~= L
   error('The window length must be the same as the filter length.')
end
if nw == 0   % replace the following with the default window of your choice.
   Wind = .54 - .46*cos(2*pi*(0:L-1)'/(L-1));
end
%
% to use Kaiser window, beta must be supplied
% att = 60; % dB of attenuation desired in sidelobe
% beta = 0.1102*(att-8.7);
% wind = kaiser(L,beta);

if  any( Wn<0 | Wn>1 )
   error('Frequencies must fall in range between 0 and 1.')
end
if  any(diff(Wn)<0)
   error('Frequencies must be increasing')
end

Wn = Wn(:)';
ff = [0,Wn(1:nbands-1); Wn(1:nbands-1),1];
mags = [mags(:)'; mags(:)'];
hh = firls(L-1,ff(:),mags(:));

b = hh.*Wind(:)'; 
a = 1;

if SCALING
    if First_Band
        b = b / sum(b);  % unity gain at DC
    else
        if ff(4)==1
            % unity gain at Fs/2
            f0 = 1;
        else
            % unity gain at center of first passband
            f0 = mean(ff(3:4));
        end
        b = b / abs( exp(-j*2*pi*(0:L-1)*(f0/2))*(b.') );
    end
end

