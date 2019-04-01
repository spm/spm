function [xY] = spm_voice_ff(Y,FS)
% decomposition at fundamental frequency
% FORMAT [xY] = spm_voice_ff(Y,FS)
%
% Y     - timeseries
% FS    - sampling frequency
%
% expects
% VOX.F0 - fundamental frequency (glottal pulse rate)
% VOX.F1 - format frequency
%
% output structure
%--------------------------------------------------------------------------
% xY.Y  - timeseries
% xY.Q  - parameters - lexical
% xY.P  - parameters - prosidy
% 
% xY.P.amp - amplitude
% xY.P.ff0 - fundamental frequency (Hz)
% xY.P.ff1 - format frequency (Hz)
% xY.P.dur - duration (seconds)
% xY.P.tim - timbre
% xY.P.inf - inflection
%
% This routine transforms a timeseries using a series of discrete cosine
% transforms and segmentations into a set of lexical and prosody
% parameters. Effectively, this is a rather complicated sequence of
% straightforward operations that constitute a parameterised nonlinear
% mapping from a parameter space to a timeseries corresponding to a spoken
% word. In brief, the transform involves identifying the interval
% containing the words spectral energy and dividing it up into a sequence
% of  fundamental segments (at the fundamental frequency or glottal pulse
% rate). The spectral content (or form of transient) for each segment is
% characterised in terms of the cross covariance function whose length is
% determined by the fundamental format frequency. The resulting matrix  its
% parameterised with even functions based upon a discrete cosine transform.
% Because the basis functions are even (i.e. symmetrical) the resulting
% coefficients nonnegative. In turn, this allows a log transform  and
% subsequent normalisation, by a scaling (timbre) parameter. The normalised
% log format coefficients are finally parameterised using two discrete
% cosine transforms over time, within and between segments, respectively.
% This provides a sufficiently rich parameterisation to generate reasonably
% realistic timeseries. The  fluctuations in the fundamental frequency
% between segments are parameterised with another discrete cosine
% transform. This has two key parameters that model inflection. Please see
% the annotated code below for further details.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_ff.m 7562 2019-04-01 09:49:28Z karl $

% defaults
%--------------------------------------------------------------------------
global VOX
try, Nu = VOX.Nu; catch, Nu  = 32;    end    % DCT order   (formants)
try, Nv = VOX.Nv; catch, Nv  = 8;     end    % DCT order   (interval)
try, Tu = VOX.Tu; catch, Tu  = 4;     end    % log scaling (formants)
try, Tv = VOX.Tv; catch, Tv  = 1;     end    % log scaling (interval)
try, F0 = VOX.F0; catch, F0  = 96;    end    % fundamental frequency
try, F1 = VOX.F1; catch, F1  = 32;    end    % formant frequency
try, F2 = VOX.F2; catch, F2  = 1024;  end    % minimum formant

% parameterise fundamental frequency modulations
%==========================================================================
I     = spm_voice_frequency(Y,FS,F0)/FS;     % get intervals
nI    = length(I);                           % number of intervals
D     = spm_dctmtx(nI,3);                    % inflection basis set
p     = sqrt(nI)\I*D;                        % fluctuations around mean
dI    = sqrt(nI)*p*D';                       % parameterised intervals
I     = round([1 FS*cumsum(dI)]);            % starting from one


% unwrap fundamental segments using cross covariance functions (ccf) and
% apply DCT to formant frequency modulation
%--------------------------------------------------------------------------
Ny    = numel(Y);                            % number of samples
Ni    = round(8192/F1);                      % basis set to 8000 Hz
nj    = round(FS/F1);                        % formant interval length
ni    = numel(find((I + nj) < Ny));          % number of intervals
D     = spm_dctmtx(2*nj + 1,Ni*4);           % discrete cosine transform
D     = D*kron(speye(Ni,Ni),[1 0 -1 0]');    % retain even functions
Q     = zeros(Ni,ni);
for j = 1:ni
    ii     = I(j) + (0:nj);
    ccf    = xcov(Y(ii));
    Q(:,j) = D'*ccf;
end

% remove noise
%--------------------------------------------------------------------------
% Q = bsxfun(@minus,Q,Q(:,1));

% log transform (nonnegative) coefficients  and parameterise with a pair of
% discrete cosine transforms over formant frequnecies and intervals
% respectively to create formant parameters (Q)
%--------------------------------------------------------------------------
Q  = log(abs(Q) + eps);                      % log transform
Q  = Q - mean(Q(:));                         % detrend
B  = 1 - exp(-(1:Ni)*F1/F2);                 % auditory range
Q  = bsxfun(@times,Q,B(:));                  % balance
S  = std(Q(:));                              % timbre
Q  = Q/std(Q(:));                            % normalise

nu = exp(-Tu*(0:(Ni - 1))/Ni);               % log spacing
nv = exp(-Tv*(0:(ni - 1))/ni);               % log spacing
nu = nu - min(nu); nu = Ni - (Ni - 1)*nu/max(nu);
nv = nv - min(nv); nv = ni - (ni - 1)*nv/max(nv);
U  = spm_dctmtx(Ni,Nu,nu);                   % DCT over formants
V  = spm_dctmtx(ni,Nv,nv);                   % DCT over intervals
Q  = U\Q/V';                                 % coeficients

% assemble prosody parameters
%--------------------------------------------------------------------------
P.amp = log(max(Y));                         % amplitude (a.u.)
P.lat = log(1/32);                           % latency (sec)
P.ff1 = log(F1);                             % format frequency (Hz)
P.dur = log(Ny/FS);                          % duration (seconds)
P.tim = log(S);                              % timbre
P.inf = p;                                   % inflection

% output structure
%--------------------------------------------------------------------------
xY.Y  = Y;                                   % timeseries
xY.Q  = Q;                                   % parameters - lexical
xY.P  = P;                                   % parameters - prosidy
xY.i  = [1,Ny];                              % range – indices      

return

% uncomment 'return' for graphics
%==========================================================================
subplot(2,2,1), plot((1:Ny)/FS,Y), axis square, spm_axis tight
xlabel('time (sec)'), ylabel('amplitude')
title('Timeseries','FontSize',16)

subplot(2,2,2), imagesc((1:ni)/F0,1000*[-nj,nj]/FS,D*exp(S*U*Q*V'))
axis square, xlabel('time (seconds)'), ylabel('time (ms)')
title('transients'), set(gca,'YLim',[-16 16])

subplot(4,2,6), imagesc((1:ni)/F0,(1:Ni)*F1,U*Q*V')
xlabel('time (seconds)'), ylabel('Formants (Hz)')
title('Spectral decomposition','FontSize',16)

subplot(4,2,8), imagesc((1:ni)/F0,(1:Ni)*F1,U*Q*V')
xlabel('time (seconds)'), ylabel('Formants (Hz)')
title('Spectral decomposition','FontSize',16)

subplot(4,2,5), imagesc(Q)
xlabel('coefficients'), ylabel('coefficients')
title('Parameters','FontSize',16)

subplot(4,2,7), plot(1./dI), axis square, spm_axis tight
xlabel('time (intervals)'), ylabel('fundamental frequency')
title('Inflection','FontSize',16), drawnow

return

% auxiliary code
%==========================================================================

% snap-to grid
%--------------------------------------------------------------------------
Qi   = std(Q,[],2);
Qj   = std(Q,[],1);
i    = spm_voice_warp(Qi,6);
j    = spm_voice_warp(Qj,3);
Q    = Q(i,j);



