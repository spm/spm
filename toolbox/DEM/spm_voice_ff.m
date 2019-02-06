function [xY] = spm_voice_ff(Y,FS,F0)
% decomposition at fundamental frequency
% FORMAT [xY] = spm_voice_ff(Y,FS,F0)
%
% Y    - timeseries
% FS   - sampling frequency
% F0   - fundamental frequency (glottal pulse rate)
%
% output structure
%--------------------------------------------------------------------------
% xY.Y  -   timeseries
% xY.Q  -   parameters - lexical
% xY.P  -   parameters - prosidy
% 
% xY.P.amp -   amplitude
% xY.P.dur -   duration (seconds)
% xY.P.ff0 -   fundamental frequency (Hz)
% xY.P.ff1 -   format frequency (Hz)
% xY.P.tim -   timbre
% xY.P.inf -   inflection
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
% $Id: spm_voice_ff.m 7528 2019-02-06 19:19:49Z karl $


% defaults
%--------------------------------------------------------------------------
if nargin < 3, F0 = 100; end

% find periods of spectral energy, extract and normalise
%--------------------------------------------------------------------------
i     = spm_voice_onset(Y,FS);                % interval containing word
Y     = Y(i);
A     = max(abs(Y));
Y     = Y/A;
ny    = length(i);

% parameterise fundamental frequency modulations
%==========================================================================
[I,DJ] = spm_voice_frequency(Y,FS,F0);        % get intervals

DI    = mean(I);                              % average interval
I     = I/DI;                                 % normalised intervals
nI    = length(I);                            % number of intervals
D     = spm_dctmtx(nI,3);                     % inflection basis set
p     = I*D;                                  % fluctuations
dI    = p*D'*DI;                              % parameterised intervals
I     = round([1 cumsum(dI)]);                % starting from one

% unwrap fundamental segments using cross covariance functions (ccf)
%--------------------------------------------------------------------------
N     = 64;                                   % order of basis set
ni    = numel(I) - 4;                         % number of intervals
nj    = round(DJ);                            % interval length
D     = spm_dctmtx(2*nj + 1,N*4);             % discrete cosine transform
D     = D*kron(speye(N,N),[1 0 -1 0]');       % retain even functions
Q     = zeros(N,ni);
for j = 1:ni
    ii     = I(j) + (0:2*nj);
    ccf    = xcov(Y(ii),nj);
    Q(:,j) = D'*ccf;
end

% log transform (nonnegative) coefficients  and parameterise with a pair
% of discrete cosine transforms over time, within and between intervals
% respectively to create formant parameters (Q)
%--------------------------------------------------------------------------
Q    = abs(Q/mean(Q(:)));
Q    = log(Q + exp(-8));
Q    = Q - mean(Q(:));
U    = spm_dctmtx(N,16);
V    = spm_dctmtx(ni,8);
Q    = U'*Q*V;
Q    = Q - mean(Q(:));
S    = std(Q(:));
Q    = Q/S;
 
% assemble prosody parameters
%--------------------------------------------------------------------------
P.amp = log(A);                              % amplitude
P.dur = log(ny/FS);                          % duration (seconds)
P.ff0 = log(FS/DI);                          % fundamental frequency (Hz)
P.ff1 = log(FS/DJ);                          % format frequency (Hz)
P.tim = log(S);                              % timbre
P.inf = p/p(1);                              % inflection

% output structure
%--------------------------------------------------------------------------
xY.Y = Y;                                    % timeseries
xY.Q = Q;                                    % parameters - lexical
xY.P = P;                                    % parameters - prosidy


return

% uncomment 'return' for graphics
%==========================================================================
pst = (1:ny)/FS;
subplot(2,2,1), plot(pst,Y), axis square, spm_axis tight
xlabel('time (sec)'), ylabel('amplitude')
title('Timeseries','FontSize',16)

subplot(2,2,2), imagesc(D*exp(S*U*Q*V')), axis square
xlabel('time (seconds)'), ylabel('time (bins)'), title('transients')

subplot(4,2,6), imagesc(U*Q*V')
xlabel('time (intervals)'), ylabel('time (bins)')
title('Spectral decomposition','FontSize',16)

subplot(4,2,8), imagesc(Q)
xlabel('coefficients'), ylabel('coefficients')
title('Parameters','FontSize',16)

subplot(2,2,3), plot(FS./dI), axis square, spm_axis tight
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



