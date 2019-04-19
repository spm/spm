function [E,PST] = spm_voice_segmentation(wfile,SEG)
% Retrieves likelihoods from an audio device or file
% FORMAT [EEG,PST] = spm_voice_segmentation(wfile,SEG)
%
% wfile  - .wav file or audiorecorder object
%
% This routine plots the timeseries after segmentation and word recognition
% as implemented by spm_voice_read
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_segmentation.m 7574 2019-04-19 20:38:15Z karl $

%% get  parameters from VOX
%==========================================================================
global VOX

% get source (recorder) and FS
%--------------------------------------------------------------------------
[FS,read] = spm_voice_FS(wfile);


%% read file and plot timeseries (and envelope)
%==========================================================================
Y   = read(wfile);
j   = fix((1:(SEG(end).IT) + FS/2));
Y   = Y(j(j < numel(Y)));
G   = spm_voice_check(Y,FS,VOX.C);
pst = (1:numel(Y))/FS;

subplot(4,1,1)
plot(pst,Y,'c'), spm_axis tight
xlabel('time (sec)'), ylabel('amplitude')
title('Acoustic signal','FontSize',16),hold on, box off

subplot(4,1,2)
plot(pst,G,'k',pst,spm_zeros(pst) + VOX.U,':r'), spm_axis tight
xlabel('time (sec)'), ylabel('power')
title('Spectral envelope','FontSize',16), box off

subplot(4,1,3), imagesc(full(spm_cat({SEG.P})))
set(gca,'YTickLabel',{VOX.PRO.str})
xlabel('word'), ylabel('prodisy')
title('Prodisy','FontSize',16), box off

% scan through words
%--------------------------------------------------------------------------
rng('default')
for w = 1:numel(SEG)
    
    % colour 
    %----------------------------------------------------------------------
    col = spm_softmax(randn(3,1));
    
    % retrieve epoch 
    %----------------------------------------------------------------------
    i    = round(SEG(w).I0:SEG(w).IT);
    j    = round(SEG(w).I0); 
 
    % plot and label
    %----------------------------------------------------------------------
    subplot(4,1,1), plot(i/FS,Y(i),'Color',col)
    subplot(4,1,2), text(j/FS,G(j) + VOX.U,SEG(w).str,'Color',col,'FontWeight','bold')

end

% return

%% simulated EEG (i.e. prediction error) responses - discrete updating
%==========================================================================
ni    = 128;                                  % number of iterations
di    = FS/2;                                 % time bins (16 ms.)
nw    = numel(VOX.LEX);                       % number of words
E     = zeros(numel(Y),nw);                   % length of time series
Q     = zeros(numel(Y),nw);                   % length of time series
D     = spm_dctmtx(di,ni);
D     = D*spm_dctmtx(ni,ni)';
for w = 1:numel(SEG)
    
        % gradient descent free energy
        %==================================================================
        L   = log(SEG(w).L{1} + exp(-16));    % posterior
        v   = log(SEG(w).p    + exp(-16));    % prior
        
        % evidence accumulation
        %------------------------------------------------------------------
        s   = spm_softmax(v);
        for i = 1:ni
            v      = L - log(s);
            v      = v - mean(v)*(1 - 1/8);
            s      = spm_softmax(log(s) + 4*v/ni);
            e(:,i) = v;
            q(:,i) = s;
        end
        
        % place in timeseries
        %------------------------------------------------------------------
        j      = (1:di) + SEG(w).IT;
        E(j,:) = E(j,:) + D*e';      
        Q(j,:) = D*q';
end

%% filter to simulate EEG (between one and 16 Hz)
%--------------------------------------------------------------------------
d = 32;                                       % decimation for plotting
Q = Q(1:d:end,:);
E = E(1:d:end,:);
E = E - spm_conv(E,FS/d,0);
E = spm_conv(E,FS/d/16,0);

% plot click
%--------------------------------------------------------------------------
PST = pst(1:d:end);
subplot(4,1,3), imagesc(PST,1:nw,(1 - Q'))
set(gca,'YTick',1:nw,'YTickLabel',{VOX.LEX.word})
xlabel('time (seconds)'), ylabel('word')
title('Simulated neuronal firing','FontSize',16)

subplot(4,1,4), plot(PST,E) %,':',PST,std(E,[],2))
xlabel('time (seconds)'), ylabel('a.u.'), spm_axis tight
title('Simulated EEG','FontSize',16), box off, set(gca,'YLim',[-1/4,1])

% place EEG in VOX
%--------------------------------------------------------------------------
VOX.EEG = E;
VOX.PST = PST;




