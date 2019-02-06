function [Y] = spm_voice_iff(xY,FS,DT)
% inverse decomposition at fundamental frequency
% FORMAT [Y] = spm_voice_iff(xY,FS,DT)
%
% xY    -  cell array of word structures
% xY.Q  -   parameters - lexical
% xY.P  -   parameters - prosidy
% xY.P.amp  -   amplitude
% xY.P.dur  -   duration (seconds)
% xY.P.ff0  -   fundamental frequency (Hz)
% xY.P.ff1  -   format frequency (Hz)
% xY.P.tim  -   timbre
% xY.P.inf  -   inflection
%
% FS    - sampling frequency
% DT    - latency or pause (seconds) [optional]
%
% Y     - reconstructed timeseries
%
% This routine recomposes a timeseries from temporal basis sets at the
% fundamental frequency.  in other words, it applies the reverse sequence
% of inverse transforms implemented by spm_voice_ff.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_iff.m 7528 2019-02-06 19:19:49Z karl $

% check for structure arrays
%--------------------------------------------------------------------------
if numel(xY) > 1
    for s = 1:numel(xY)
        try
            L = DT(s);
        catch
            L = 1/4;
        end
        
        % recompose and play
        %------------------------------------------------------------------
        spm_voice_iff(xY(s),FS,L);
        
    end
    return
end

% reconstitute sample points
%--------------------------------------------------------------------------
A  = exp(xY.P.amp);                         % amplitude (a.u.)
T  = exp(xY.P.dur);                         % duration (seconds)
F0 = exp(xY.P.ff0);                         % fundamental frequency (Hz)
F1 = exp(xY.P.ff1);                         % formant frequency (Hz)
S  = exp(xY.P.tim);                         % timbre
P  = xY.P.inf;                              % inflection

% reconstitute intervals
%--------------------------------------------------------------------------
DI = round(FS/F0);                          % average fundamental interval
nI = round(T*F0);                           % number of intervals
D  = spm_dctmtx(nI - 1,numel(P));           % basis set for inflection
dI = D*P(:);                                % fluctuations
dI = DI*dI/mean(dI);                        % scaled fluctuations
I  = round([1; cumsum(dI)]);                % cumulative intervals

% reconstitute format coefficients
%--------------------------------------------------------------------------
N  = 64;                                    % order of DCT
ni = numel(I) - 4;                          % number of intervals
nj = round(FS/F1);                          % interval length

Q  = xY.Q;
Q  = Q - mean(Q(:));
Q  = S*Q/std(Q(:));
V  = spm_dctmtx(ni,size(Q,2));
U  = spm_dctmtx(N, size(Q,1));
Q  = exp(U*Q*V');

% reconstitute timeseries
%--------------------------------------------------------------------------
jj = 0:(2*nj);
D  = spm_dctmtx(numel(jj),N*4);
D  = D*kron(speye(N,N),[1 0 -1 0]');
Y  = zeros(I(end) + nj,1);
for j = 1:ni
    ii    = I(j)  + jj;
    Y(ii) = Y(ii) + D*Q(:,j);
end

% scale amplitude
%--------------------------------------------------------------------------
Y  = A*Y/max(abs(Y));


% play timeseries if requested
%--------------------------------------------------------------------------
global voice_options
if ~ voice_options.mute
    if nargin < 3, DT = 0; end
    pause(DT),sound(Y,FS);
end

% graphics  if requested
%--------------------------------------------------------------------------
if voice_options.graphics
    
    % peristimulus time (seconds) and plot
    %--------------------------------------------------------------------------
    pst = (1:numel(Y))/FS;
    subplot(2,2,1), plot(pst,Y), axis square, spm_axis tight
    xlabel('time (sec)'), ylabel('amplitude')
    title('Timeseries','FontSize',16)
    
    subplot(2,2,2), imagesc(D*Q), axis square
    xlabel('time (intervals)'), ylabel('time (bins)')
    title('Spectral decomposition','FontSize',16)
    
    subplot(4,2,6), imagesc(log(Q + exp(-8)))
    ylabel('frequency'), title('Formants','FontSize',16)
    
    subplot(4,2,8), imagesc(xY.Q)
    xlabel('coefficients'), ylabel('coefficients')
    title('Parameters','FontSize',16)
    
    subplot(2,2,3), plot(FS./dI), axis square, spm_axis tight
    xlabel('time (intervals)'), ylabel('fundamental frequency')
    title('Inflection','FontSize',16), drawnow
    
end





