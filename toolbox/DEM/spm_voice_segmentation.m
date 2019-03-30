function spm_voice_segmentation(wfile,SEG)
% Retrieves likelihoods from an audio device or file
% FORMAT spm_voice_segmentation(wfile,I0,str)
%
% wfile  - .wav file or audiorecorder object
% IO     - pointers to words
% str    - word string
%
% This routine plots the timeseries after segmentation and word recognition
% as implemented by spm_voice_read
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_segmentation.m 7561 2019-03-30 10:39:07Z karl $

%% get  parameters from VOX
%==========================================================================
global VOX
VOX.onsets = 0;

% get source (recorder) and FS
%--------------------------------------------------------------------------
if isa(wfile,'audiorecorder')

    FS    = get(wfile,'SampleRate');
    read  = @getaudiodata;
    
else
    
    % sound file (read)
    %----------------------------------------------------------------------
    try
        xI     = audioinfo(wfile);
        FS     = xI.SampleRate;
        read   = @audioread;
    catch
        [Y,FS] = wavread(wfile,[1 1]);
        read   = @wavread;
    end
end


%% read file and plot timeseries (and envelope)
%==========================================================================
spm_figure('GetWin','Segmentation'); clf;

Y   = read(wfile);
Y   = Y(round(1:(SEG(end).IT + FS/2)));
G   = spm_conv(abs(Y),FS*VOX.C);
G   = G - min(G);
pst = (1:numel(Y))/FS;

subplot(3,1,1)
plot(pst,Y,'c'), spm_axis tight
xlabel('time (sec)'), ylabel('amplitude')
title('Timeseries','FontSize',16),hold on

subplot(3,1,2)
plot(pst,G,'k',pst,spm_zeros(pst) + VOX.U,':r'), spm_axis tight
xlabel('time (sec)'), ylabel('power')
title('Envelope','FontSize',16)

subplot(6,1,5), imagesc(full(spm_cat({SEG.P})))
set(gca,'YTickLabel',{VOX.PRO.str})
xlabel('word'), ylabel('prodisy')
title('Prodisy','FontSize',16)

subplot(12,1,12), imagesc(full(spm_cat({SEG.R})))
set(gca,'YTickLabel',{VOX.WHO.str})
xlabel('word'), ylabel('frequency')
title('Identity','FontSize',16)

% scan through words
%--------------------------------------------------------------------------
for w = 1:numel(SEG)
    
    % colour 
    %----------------------------------------------------------------------
    col = spm_softmax(randn(3,1));
    
    % retrieve epoch 
    %----------------------------------------------------------------------
    i    = round(SEG(w).I0:SEG(w).IT);
    j    = round(SEG(w).I0 + FS/4); 
 
    % plot and label
    %----------------------------------------------------------------------
    subplot(3,1,1), plot(i/FS,Y(i),'Color',col)
    
    subplot(3,1,2), text(j/FS,G(j) + VOX.U,SEG(w).str,'Color',col,'FontWeight','bold')

end


