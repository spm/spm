function [xY] = spm_voice_ff(Y,FS)
% decomposition at fundamental frequency
% FORMAT [xY] = spm_voice_ff(Y,FS)
%
% file  - .wav file
%
% xY.Y  - timeseries
% xY.i  - timing
% xY.y  - unwrapped timeseries
% xY.FS - sampling frequency
% xY.ni - number of samples
% xY.ti - last sample (seconds)
% xY.ci - coefficients
%
% This routine decomposes a timeseries into a temporal basis set at the
% fundamental frequency
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_ff.m 7512 2019-01-05 21:15:16Z karl $

% find Sigma points (maxima of Hilbert transform)
%==========================================================================

% cross-correlation function of Hilbert transform
%--------------------------------------------------------------------------
H     = abs(hilbert(Y));
ccf   = fftshift(xcorr(H,FS/32,'coeff'));
k     = 1:numel(ccf);
[d,i] = max(ccf'.*(k > FS/128 & k < FS/64));

% find successive maxima
%--------------------------------------------------------------------------
dI    = i;
I     = 1;
for i = 1:256
    di    = (1:(2*dI))';                      % search range
    s     = (dI/8)^2;                         % window s.d.
    g     = exp(-(di - dI).^2/(2*s));         % Gaussian search window
    try
        [d,j] = max((H(di + I(i))).*g);       % local maximum
    catch
        break                                 % break if no maxima
    end
    dI       = min(max(j,FS/128),FS/64);      % interval
    I(i + 1) = I(i) + round(dI);              % next interval
end

% great sampling disease
%--------------------------------------------------------------------------
I     = round((I(1:end - 1) + I(2:end))/2);
s0    = 16*FS/1000;                           % smoothing (milliseconds)
sY    = spm_conv(abs(Y),s0);                  % smoothed power
i     = find(sY/max(sY) > 1/16,1,'last');     % power offset

J     = I(I < i);                             % indices of extrema
nJ    = length(J);                            % number of sample points
tJ    = J(end);                               % time length
J0    = linspace(1,tJ,nJ);                    % uniform increments
D     = spm_dctmtx(nJ,8);                     % basis set
cJ    = (J0 - J)*D;                           % fluctuations
tJ    = tJ/FS;                                % time in seconds


% unwrap fundamental segments
%--------------------------------------------------------------------------
N     = 32;
M     = 8;
i     = J;
ni    = numel(i) - 1;
y     = zeros(N,ni);
for j = 1:ni
    ii     = i(j):(i(j + 1) - 1);
    nii    = numel(ii);
    H      = spm_hanning(nii);
    [d,k]  = max(H.*Y(ii));
    k      = k - numel(ii)/2;
    jj     = max(1,ii + round(k));
    D      = spm_dctmtx(length(jj),N*4);
    D      = D*kron(speye(N,N),[1 0 -1 0]');
    Yj     = H.*(Y(jj) + flipud(Y(jj)));
    y(:,j) = D'*Yj;
end

% temporal smoothing (within epoch)
%--------------------------------------------------------------------------
D     = spm_dctmtx(ni,M);
y     = y*D;
y     = y/sum(y(:).^2);

% output structure
%--------------------------------------------------------------------------
xY.Y  = Y;                                   % timeseries
xY.y  = y;                                   % unwrapped timeseries
xY.i  = i;                                   % timing
xY.FS = FS;                                  % sampling frequency
xY.ni = nJ;                                  % number of samples
xY.ti = tJ;                                  % last sample
xY.ci = cJ;                                  % coefficients


return

% graphics
%--------------------------------------------------------------------------
subplot(2,1,1)
j    = (1:5000);
b    = spm_zeros(Y);
b(I) = max(Y);
pst  = 1000*j/FS;
plot(pst,Y(j),pst,b(j),':')






