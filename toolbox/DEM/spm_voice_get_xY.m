function [xY,word,NI] = spm_voice_get_xY(PATH)
% Create word arrays from sound file exemplars
% FORMAT [xY,word,NI] = spm_voice_get_xY(PATH)
%
% PATH      -  directory containing sound files of exemplar words
%
% xY(nw,ns) -  structure array for ns samples of nw words
% word(nw)  -  cell array of word names
% NI(nw,ns) -  numeric array of number of minima
%
%  This routine uses a library of sound files, each containing 32 words
%  spoken with varying prosody. The name of the sound file labels the word
%  in question. These exemplars are then transformed (using a series of
%  discrete cosine and Hilbert transforms) into a set of parameters, which
%  summarise the lexical content and prosody. The inverse transform
%  generates  timeseries that can be played to articulate a word. The
%  transform operates on a word structure xY to create lexical and prosody
%  parameters (Q and P respectively).
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_get_xY.m 7766 2020-01-05 21:37:39Z karl $


% get corpus
%==========================================================================
rng('default')

% get the list of words and sampling frequency (FS)
%--------------------------------------------------------------------------
spm_figure('GetWin','voice'); clf
cd(PATH)
wfile    = dir('*.wav');
try
    audioread(wfile(1).name,[1,1]);
    read = @audioread;
catch
    read = @wavread;
end
[Y,FS]   = read(wfile(1).name,[1,1]);

% assemble cell array of word structures for subsequent characterisation
%==========================================================================
nw    = numel(wfile);                               % number of words
ns    = 32;                                         % number of samples
nj    = 16;                                         % number of jitters
sj    = FS/64;                                      % s.d. of jitter
for w = 1:nw
    
    % get lexicon name and create structure
    %----------------------------------------------------------------------
    wname   = wfile(w).name;                        % name of word
    [d,str] = fileparts(wname);
    word{w} = str; disp(str)
    
    % get the midpoint of words from the (maxima) of successive exemplars
    %----------------------------------------------------------------------
    G     = spm_voice_check(read(wname),FS,1/4);
    I     = find((diff(G(1:end - 1)) > 0) & (diff(G(2:end)) < 0));
    [i,j] = sort(G(I),'descend');
    I     = sort(I(j(1:ns)));
    
    % and plot
    %----------------------------------------------------------------------
    subplot(2,1,1), plot(G), hold on, plot(I,G(I),'ro'), hold off
    title(sprintf('Power peaks - %s',word{w}),'FontSize',16)
    xlabel('frequency (hertz)'), ylabel('power'), box off
    drawnow
    
    for s = 1:ns
        
        % retrieve (one second) epoch around midpoint and transform
        %------------------------------------------------------------------
        i     = round([-1/2 1/2]*FS + I(s));
        if i(1) < 1
            i = i - i(1) + 1;
        end
        Y     = read(wname,i);
        i     = spm_voice_onsets(Y,FS);
        i     = i{end};
        ni    = numel(Y);
        for j = 1:nj
            
            % jitter interval
            %--------------------------------------------------------------
            k       = (s - 1)*nj + j;
            ik      = fix((i(1) + sj*randn):(i(end) + sj*randn));
            ik      = ik(ik < ni & ik > 1);
            xY(w,k) = spm_voice_ff(Y(ik),FS);
            
            % record interval
            %--------------------------------------------------------------
            xY(w,k).i(1) = ik(1)/FS   - 1/2;
            xY(w,k).i(2) = ik(end)/FS - 1/2;
        end
        
        % apply inverse transform and play, if requested
        %------------------------------------------------------------------
        spm_voice_iff(xY(w,k));
        
        % record number of spectral minima intervals
        %------------------------------------------------------------------
        [Y,IS]  = spm_voice_get_next(Y);
        j       = fix((0:FS + FS/4) + IS);
        j       = j(j < ni & j > 1);
        NI(w,s) = numel(spm_voice_onsets(Y(j),FS));
        
    end
    
end
