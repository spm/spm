function [Y] = spm_voice_iff(xY,FS,GRAPHICS)
% inverse decomposition at fundamental frequency
% FORMAT [Y] = spm_voice_iff(xY,FS,GRAPHICS)
%
% xY.y  - unwrapped timeseries
% xY.FS - sampling frequency
% xY.ni - number of samples
% xY.ti - last sample (seconds)
% xY.ci - coefficients
%
% Y     - reconstructed timeseries
%
% This routine recomposes a timeseries from a temporal basis set at the
% fundamental frequency
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_iff.m 7512 2019-01-05 21:15:16Z karl $

% check for structure arrays
%--------------------------------------------------------------------------
if numel(xY) > 1
    for w = 1:numel(xY)
        Y = spm_voice_iff(xY(w),FS);
    end
    return
end

% reconstitute sample points
%--------------------------------------------------------------------------
J     = linspace(1,xY.ti*FS,xY.ni);
D     = spm_dctmtx(xY.ni, size(xY.ci,2));
dJ    = xY.ci*pinv(D);
i     = max(1,round(J - dJ));

% reconstitute timeseries
%--------------------------------------------------------------------------
[N,M] = size(xY.y);
Y     = zeros(i(end),1);
D     = spm_dctmtx(length(i) - 1,M);
y     = xY.y*pinv(D);
for j = 1:(numel(i) - 1)
    ii    = i(j):(i(j + 1) - 1);
    D     = spm_dctmtx(length(ii),N*4);
    D     = D*kron(speye(N,N),[1 0 -1 0]');
    Y(ii) = D*y(:,j);
end

% play timeseries
%--------------------------------------------------------------------------
wavplay(Y/max(Y),FS), drawnow

if nargin < 3, return, end

% graphics
%--------------------------------------------------------------------------
pst = (1:numel(Y))/FS;                        % peristimulus time (seconds)
subplot(2,2,1)
plot(pst,Y), axis square
xlabel('time (sec)'), ylabel('amplitude'), title('timeseries')

subplot(2,2,2)
imagesc(pst,1:size(xY.y,1),D*y), axis square
xlabel('time (seconds)'), ylabel('time (bins)'), title('transients')

subplot(4,2,6)
imagesc(y)
ylabel('coefficients'), title('coefficients')

subplot(4,2,8)
imagesc(xY.y)
xlabel('coefficients'), ylabel('coefficients'), title('coefficients')

subplot(2,2,3)
plot(dJ), axis square
xlabel('time (windows)'), ylabel('compression'), title('window length')






