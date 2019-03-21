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
% $Id: spm_voice_segmentation.m 7551 2019-03-21 15:10:05Z karl $

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
Y   = Y(round(1:(FS + SEG(end).I0)));
G   = spm_voice_filter(Y,FS);
G   = spm_conv(G,FS/6);
pst = (1:numel(Y))/FS;

subplot(3,1,1)
plot(pst,Y,':k'), spm_axis tight
xlabel('time (sec)'), ylabel('amplitude')
title('Timeseries','FontSize',16),hold on

subplot(3,1,2)
plot(pst,G,'k',pst,spm_zeros(pst) + 1/128,':r'), spm_axis tight
xlabel('time (sec)'), ylabel('power')
title('Envelope','FontSize',16)

subplot(3,1,3), imagesc(full(spm_cat({SEG.P})))
set(gca,'YTickLabel',{VOX.PRO.str})
xlabel('word'), ylabel('prodisy')
title('Prodisy','FontSize',16)


% scan through words
%--------------------------------------------------------------------------
for w = 1:numel(SEG)
    
    % colour 
    %----------------------------------------------------------------------
    col = spm_softmax(randn(3,1));
    
    % retrieve epoch 
    %----------------------------------------------------------------------
    i      = round(SEG(w).I0:(SEG(w).IT + SEG(w).I0));
 
    % plot and label
    %----------------------------------------------------------------------
    subplot(3,1,1), plot(i/FS,Y(i),'Color',col)
    subplot(3,1,2), text(SEG(w).I0/FS,G(SEG(w).I0) + 1/256,SEG(w).str,'Color',col,'FontWeight','bold')

end


