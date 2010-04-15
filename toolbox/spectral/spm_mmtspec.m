function [p, f, t] = spm_mmtspec (x,Fs,WinLength,nOverlap,NW)
% Moving multitaper based spectrogram
% FORMAT [p, f, t] = spm_mmtspec (x,Fs,WinLength,nOverlap,NW)
%
% x         input time series
% Fs        sampling frequency 
% WinLength length of moving window (default is 256)
% nOverlap  overlap between successive windows (default is 250)
% NW        time bandwidth parameter (e.g. 3 or 4), default 3
%
% p         p(f, t) is estimate of power at freq f and time t
% 
% Plot spectrogram using imagesc(t,f,p); axis xy
%___________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Partha Mitra and Ken Harris
% $Id: spm_mmtspec.m 3821 2010-04-15 14:26:34Z will $

if (nargin<3 | isempty(WinLength)) WinLength = 256;  end;
if (nargin<4 | isempty(nOverlap)) nOverlap = 250; end;
if (nargin<5 | isempty(NW)) NW = 3; end;

nFFT=WinLength;

WinLength=round(WinLength);
nOverlap=round(nOverlap);

Detrend = ''; 
nTapers = 2*NW -1; 

% Now do some computations that are common to all spectrogram functions

winstep = WinLength - nOverlap;

nChannels = size(x, 2);
nSamples = size(x,1);

% check for column vector input
if nSamples == 1 
	x = x';
	nSamples = size(x,1);
	nChannels = 1;
end;

if nSamples < WinLength
    disp('Error in spm_mmtspec: win length must be less than number of samples');
    return;
end

% calculate number of FFTChunks per channel
nFFTChunks = round(((nSamples-WinLength)/winstep));
% turn this into time, using the sample frequency
t = winstep*(0:(nFFTChunks-1))'/Fs;

% set up f and t arrays
if ~any(any(imag(x)))    % x purely real
	if rem(nFFT,2),    % nfft odd
		select = [1:(nFFT+1)/2];
	else
		select = [1:nFFT/2+1];
	end
	nFreqBins = length(select);
else
	select = 1:nFFT;
end
f = (select - 1)'*Fs/nFFT;


% allocate memory now to avoid nasty surprises later
y=complex(zeros(nFreqBins, nFFTChunks, nChannels, nChannels)); % output array
Periodogram = complex(zeros(nFreqBins, nTapers, nChannels, nFFTChunks)); % intermediate FFTs
Temp1 = complex(zeros(nFFT, nTapers, nFFTChunks));
Temp2 = complex(zeros(nFFT, nTapers, nFFTChunks));
Temp3 = complex(zeros(nFFT, nTapers, nFFTChunks));
eJ = complex(zeros(nFFT, nFFTChunks));

% calculate Slepian sequences.  
Tapers=spm_dpss(WinLength,NW);
Tapers=Tapers(:,1:nTapers);

% Vectorized alogrithm for computing tapered periodogram with FFT 
TaperingArray = repmat(Tapers, [1 1 nChannels]);
for j=1:nFFTChunks
	Segment = x((j-1)*winstep+[1:WinLength], :);
	if (~isempty(Detrend))
		Segment = detrend(Segment, Detrend);
	end;
	SegmentsArray = permute(repmat(Segment, [1 1 nTapers]), [1 3 2]);
	TaperedSegments = TaperingArray .* SegmentsArray;
						
	fftOut = fft(TaperedSegments,nFFT);
	Periodogram(:,:,:,j) = fftOut(select,:,:); %fft(TaperedSegments,nFFT);
end	
	
% Now make cross-products of them to fill cross-spectrum matrix
for Ch1 = 1:nChannels
	for Ch2 = Ch1:nChannels % don't compute cross-spectra twice
		Temp1 = reshape(Periodogram(:,:,Ch1,:), [nFreqBins,nTapers,nFFTChunks]);
		Temp2 = reshape(Periodogram(:,:,Ch2,:), [nFreqBins,nTapers,nFFTChunks]);
		Temp2 = conj(Temp2);
		Temp3 = Temp1 .* Temp2;
		eJ=sum(Temp3, 2);
        p(:,:, Ch1, Ch2)= eJ/nTapers;
		
		% for off-diagonal elements copy into bottom half of matrix
		if (Ch1 ~= Ch2)
            p(:,:, Ch2, Ch1) = conj(eJ) / nTapers;
		end
	end
end



