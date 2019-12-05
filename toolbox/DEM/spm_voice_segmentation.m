function [E,  PST] = spm_voice_segmentation(wfile,SEG)
% Plot the results of a segmented sound fileor audio stream
% FORMAT [EEG,PST] = spm_voice_segmentation(wfile,SEG)
%
% wfile      - (double) timeseries, .wav file or audiorecorder object
%
% SEG(s).str - lexical class
% SEG(s).p   - prior
% SEG(s).L   - posterior
% SEG(s).P   - prosody class
% SEG(s).R   - speaker class
% SEG(s).I0  - first index
% SEG(s).IT  - final index
%
% EEG        - simulated EEG for each lexical entry
% PST        - corresponding peristimulus times for plotting
%
% This routine plots the timeseries after segmentation and word recognition
% as implemented by spm_voice_read. It also returns simulated belief
% updating in the form of local field potentials or EEG for simulation
% purposes.
%
% EEG and PST are also placed in the global VOX structure.
%
% see also: spm_voice_read.m
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_segmentation.m 7750 2019-12-05 17:54:29Z spm $


%% get  parameters from VOX
%==========================================================================
global VOX

% get source (recorder) and FS
%--------------------------------------------------------------------------
[FS,read] = spm_voice_FS(wfile);


%% read file and plot timeseries (and envelope)
%==========================================================================
Y   = read(wfile);
n   = numel(Y);
j   = fix((1:(SEG(end).IT) + FS/2));
Y   = Y(j(j < n));
G   = spm_voice_check(Y,FS,1/32);
pst = (1:numel(Y))/FS;

subplot(4,1,1)
plot(pst,Y,'c'), spm_axis tight
xlabel('time (seconds)'), ylabel('amplitude')
title('Acoustic signal','FontSize',16),  hold on, box off

subplot(4,1,2)
plot(pst,G,'k',pst,spm_zeros(pst) + VOX.U,':r'), hold on
xlabel('time (seconds)'), ylabel('power'), 
title('Spectral envelope','FontSize',16), spm_axis tight, box off

subplot(6,1,5), imagesc(full(spm_cat({SEG.P})))
str = {VOX.PRO.str};
set(gca,'YTick',1:numel(str),'YTickLabel',str)
xlabel('word'), ylabel('prodisy'), colorbar
title('Prodisy','FontSize',16), box off

subplot(8,1,8), imagesc(full(spm_cat({SEG.R})))
str = {VOX.WHO.str};
set(gca,'YTick',1:numel(str),'YTickLabel',str)
xlabel('word'), ylabel('idenity'), colorbar
title('Idenity','FontSize',16), box off

% scan through words
%--------------------------------------------------------------------------
rng('default')
M     = max(G);
ns    = numel(SEG);
for w = 1:ns
    
    % colour 
    %----------------------------------------------------------------------
    col  = spm_softmax(randn(3,1));
    
    % retrieve epoch 
    %----------------------------------------------------------------------
    i    = SEG(w).I0:SEG(w).IT;
    j    = SEG(w).I0; 
    i    = i(i < n & i > 1);
 
    % plot and label spectral segmentation
    %----------------------------------------------------------------------
    subplot(4,1,1), plot(i/FS,Y(i),'Color',col)
    subplot(4,1,2), text(j/FS,M*w/ns,SEG(w).str,'Color',col,'FontWeight','bold')
    
    % plot boundaries and peaks
    %----------------------------------------------------------------------
    i     = fix(SEG(w).I(1) + FS/2);
    subplot(4,1,2), plot(pst(i),G(i),'.', 'Color',col,'MarkerSize',24)
    plot([1,1]*pst(max(0 + 1,SEG(w).I0)),[0 M],'-', 'Color',col)
    plot([1,1]*pst(min(n - 1,SEG(w).IT)),[0 M],'-.','Color',col)
    for j = 2:numel(SEG(w).I)
        i = SEG(w).I(j);
        plot(pst(i),G(i),'.', 'Color',col,'MarkerSize',16)
    end

end

% return if just prosody is requested
%----------------------------------------------------------------------
if ~nargout, return, end

%% simulated EEG (i.e. prediction error) responses - discrete updating
%==========================================================================
ni    = 128;                                  % number of iterations
di    = FS/2;                                 % time bins (16 ms.)
nw    = numel(VOX.LEX);                       % number of words
E     = zeros(numel(Y),nw);                   % length of time series
Q     = zeros(numel(Y),nw);                   % length of time series
D     = spm_dctmtx(di,ni);
D     = D*spm_dctmtx(ni,ni)';
for w = 1:ns
    
        % gradient descent free energy
        %==================================================================
        L   = spm_softmax(SEG(w).L{1});
        L   = log(L        + exp(-16));       % posterior
        v   = log(SEG(w).p + exp(-16));       % prior
        
        % evidence accumulation
        %------------------------------------------------------------------
        s   = spm_softmax(v);
        for i = 1:ni
            v      = L - log(s);
            v      = v - mean(v);
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

% plot simulated belief updating
%--------------------------------------------------------------------------
PST = pst(1:d:end);
subplot(4,1,3), imagesc(PST,1:nw,(1 - Q'))
set(gca,'YTick',1:nw,'YTickLabel',{VOX.LEX.word})
xlabel('time (seconds)'), ylabel('word')
title('Simulated neuronal firing','FontSize',16)

subplot(4,1,4), plot(PST,E) %,':',PST,std(E,[],2))
xlabel('time (seconds)'), ylabel('a.u.'), spm_axis tight
title('Simulated EEG','FontSize',16), box off, set(gca,'YLim',[-1/2,1])

% place EEG in VOX
%--------------------------------------------------------------------------
VOX.EEG = E;
VOX.PST = PST;
