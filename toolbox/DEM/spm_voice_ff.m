function [xY] = spm_voice_ff(Y,FS)
% Decomposition at fundamental frequency
% FORMAT [xY] = spm_voice_ff(Y,FS)
%
% Y      - timeseries
% FS     - sampling frequency
%
% requires the following in the global VOX structure:
% VOX.F0 - fundamental frequency (glottal pulse rate)
% VOX.F1 - format frequency
%
% output structure
%--------------------------------------------------------------------------
% xY.Y   - timeseries
% xY.W   - parameters - lexical
% xY.P   - parameters - prosody
% 
% xY.P.amp - log amplitude
% xY.P.dur - log duration (sec)
% xY.P.lat - log latency  (sec)
% xY.P.tim - log timbre (a.u.)
% xY.P.pch - log pitch  (a.u.)
% xY.P.inf - inflection (a.u.)
%
% xY.R.F0  - fundamental frequency (Hz)
% xY.R.F1  - format frequency (Hz
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
% Because the basis functions are even (i.e., symmetrical) the resulting
% coefficients are nonnegative. In turn, this allows a log transform  and
% subsequent normalisation, by a scaling (timbre) parameter. The normalised
% log format coefficients are finally parameterised using two discrete
% cosine transforms over time, within and between segments, respectively.
% This provides a sufficiently rich parameterisation to generate reasonably
% realistic timeseries. The  fluctuations in the fundamental frequency
% between segments are parameterised with another discrete cosine
% transform. This has two key parameters that model inflection. Please see
% the annotated code below for further details.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_ff.m 7750 2019-12-05 17:54:29Z spm $


% defaults
%--------------------------------------------------------------------------
global VOX
try, F0 = VOX.F0; catch, F0  = 110;        end    % fundamental frequency
try, F1 = VOX.F1; catch, F1  = 26 + F0/16; end    % formant frequency

% sampling frequency
%--------------------------------------------------------------------------
if nargin < 2, FS = spm_voice_FS; end


% parameterise fundamental frequency modulations
%==========================================================================
I     = spm_voice_frequency(Y,FS,F0)/FS;     % get intervals (seconds)
nI    = length(I);                           % number of intervals
D     = spm_dctmtx(nI,3);                    % inflection basis set
I     = I*F0;                                % normalise to fundamental
S     = sqrt(nI)\I*D;                        % fluctuations around mean
pI    = sqrt(nI)*S*D';                       % parameterised intervals
I     = fix([1 FS*cumsum(pI/F0)]);           % starting from one

% unwrap fundamental segments using cross covariance functions (ccf) and
% apply DCT to formant frequency modulation
%--------------------------------------------------------------------------
Ny    = numel(Y);                            % number of samples
Ni    = 256;                                 % number of format bins
nj    = round(FS/F1);                        % formant interval length
ni    = numel(find((I + nj) < Ny));          % number of intervals
D     = spm_dctmtx(2*nj + 1,Ni*4);           % discrete cosine transform
D     = D*kron(speye(Ni,Ni),[1 0 -1 0]');    % retain even functions
Q     = zeros(Ni,ni);
for j = 1:ni
    ii     = I(j) + (0:nj);
    ccf    = xcov(Y(ii));
    Q(:,j) = abs(D'*ccf);
end

% log transform (nonnegative) coefficients  and parameterise with a pair of
% discrete cosine transforms over formant frequencies and intervals
% respectively to create formant parameters (Q)
%--------------------------------------------------------------------------
Q     = bsxfun(@times,Q,spm_hanning(ni)');   % window
Q     = log(Q + eps)/2;                      % log transform
Q     = Q - mean(Q(:));                      % detrend
T     = std(Q(:));                           % timbre
W     = spm_voice_iQ(Q);                     % DCT transform

% assemble prosody parameters
%--------------------------------------------------------------------------
P.amp = log(max(Y));                         % amplitude (a.u.)
P.lat = log(1/32);                           % latency (sec)
P.dur = log(Ny/FS);                          % duration (seconds)
P.tim = log(T);                              % log timbre
P.pch = zeros(4,1);                          % log pitch
P.inf = S;                                   % inflection

% assemble speaker parameters
%--------------------------------------------------------------------------
R.F0  = log(F0);                             % fundamental frequency (Hz)
R.F1  = log(F1);                             % format frequency (Hz)

% output structure
%--------------------------------------------------------------------------
xY.Y  = Y;                                   % timeseries
xY.W  = W;                                   % parameters - lexical
xY.P  = P;                                   % parameters - prosody
xY.R  = R;                                   % parameters - speaker
xY.i  = [1,Ny];                              % range - indices      


% uncomment 'return' for graphics
%==========================================================================
if ~ isfield(VOX,'analysis'), return, end
if VOX.analysis
    
    % figure
    %----------------------------------------------------------------------
    spm_figure('GetWin','Voice (analysis)'); clf;
    
    subplot(2,2,1), plot((1:Ny)/FS,Y), axis square, spm_axis tight
    xlabel('time (sec)'), ylabel('amplitude')
    title('Timeseries','FontSize',16)
    
    subplot(2,2,2), imagesc((1:ni)/F0,1000*[-nj,nj]/FS,D*exp(Q))
    axis square, xlabel('time (seconds)'), ylabel('time (ms)')
    title('transients'), set(gca,'YLim',[-8 8])
    
    subplot(4,2,5), imagesc(W), axis square
    xlabel('coefficients'), ylabel('coefficients')
    title('Parameters','FontSize',16)
    
    subplot(4,2,6), imagesc((1:ni)/F0,(1:Ni)*F1,exp(Q))
    xlabel('time (seconds)'), ylabel('Formants (Hz)')
    title('Spectral decomposition','FontSize',16)
    
    subplot(4,2,7), plot(F0./pI), axis square, spm_axis tight
    xlabel('time (intervals)'), ylabel('fundamental frequency')
    title('Inflection','FontSize',16)
    
    subplot(4,2,8), imagesc((1:ni)/F0,(1:Ni)*F1,Q/T)
    xlabel('time (seconds)'), ylabel('Formants (Hz)')
    title('Spectral (log) energy','FontSize',16), drawnow
    
end

return

% auxiliary code (not used)
%==========================================================================

% snap-to grid
%--------------------------------------------------------------------------
Qi   = std(Q,[],2);
Qj   = std(Q,[],1);
i    = spm_voice_warp(Qi,6);
j    = spm_voice_warp(Qj,3);
Q    = Q(i,j);
