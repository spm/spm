function [S] = spm_voice(path)
% Create basis set for sounds
% FORMAT [S] = DEM_birdsong(file)
%
% file  - .wav file
%
% S.U   - h x 3 basis functions (Hz)
% S.V   - 3 x n basis functions (seconds)
% S.Hz  - s x 1 frequencies (Hz)
%
% Bird Song demo: These simple loads a .wav file of a real bird-song; and
% approximates the ensuing spectrogram with in terms of three
% time-frequency modes.  These modes are saved in BirdSong.mat (U) for
% illustrating DEM_demo_sequences
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice.m 7512 2019-01-05 21:15:16Z karl $

% load bird song
%==========================================================================
% if ~nargin
%     file = fullfile(fileparts(mfilename('fullpath')),'lebi3.wav');
% end

% load script and prompt for audio file - at one word per second
%--------------------------------------------------------------------------
% tfile = 'Document.txt';
% str   = textread(tfile,'%s');
% str   = unique(str);
% for s = 1:numel(str)
%     clc,  disp([str(s) 'new file']), pause
%     for i = 1:32
%         clc,disp(str(s));     pause(0.9)
%         clc,disp('    -*- '); pause(0.1)
%     end
% end

% create lexical structures for subsequent word recognition
%==========================================================================
try
    cd(path)
end
wfile      = dir('*.wav');
try
    [Y,FS] = audioread(wfile(1).name,[1,128]);
catch
    [Y,FS] = wavread(wfile(1).name,[1,128]);
end

% assemble words
%--------------------------------------------------------------------------
s0    = 64*FS/1000;                            % smoothing (milliseconds)
for w = 1:numel(wfile)
    
    % get lexicon name and create structure
    %----------------------------------------------------------------------
    wname        = wfile(w).name;
    [PATH,NAME]  = fileparts(wname);
    LEX(w,1).str = NAME;
    
    % get data features from a succession (1 Hz) exemplars
    %----------------------------------------------------------------------
    I     = 1;
    for s = 1:32
        
        % find next initial time point
        %------------------------------------------------------------------
        try
            Y    = wavread(wname,[0 4*FS] + I(s));
        catch
            Y    = wavread(wname,[0 FS] + I(s));
        end
        Y        = spm_conv(abs(Y),s0);
        i        = find(Y/max(Y) > 1/8,1);
        I(s)     = I(s) + i;
        
        % retrieve epoch and decompose
        %------------------------------------------------------------------
        Y        = wavread(wname,[0 0.6*FS] + I(s));
        xY(w,s)  = spm_voice_ff(Y,FS);
        
        % check
        %------------------------------------------------------------------
        y        = spm_voice_iff(xY(w,s),FS,1);
        I(s + 1) = I(s) + FS*(1000 - 128)/1000;
        drawnow
        
        % pause
    end
end


% temporal frequency eigenreduction
%==========================================================================
% Q     = [];
% for w = 1:size(xY,1)
%     for s = 1:size(xY,2)
%         
%         % find next initial time point
%         %----------------------------------------------------------------------
%         Q  = [Q spm_vec(xY(w,s).y)];
%         
%     end
% end
% 
% [U,S] = spm_svd(Q);
% U     = U(:,1:128);
% for w = 1:size(xY,1)
%     for s = 1:size(xY,2)
%         
%         % find next initial time point
%         %----------------------------------------------------------------------
%         xy   = xY(w,s);
%         y    = xy.y;
%         P    = (U'*y(:));
%         y(:) = U*P;
%         xy.y = y;
%         spm_voice_iff(xy,FS,1);
%         subplot(2,1,2)
%         bar(P)
%         
%     end
% end


% clustering and creation of lexical array (LEX)
%==========================================================================
K     = 8;
LEX   = LEX(:,1);
for w = 1:size(xY,1)
    
    % save priors in lexical structure
    %----------------------------------------------------------------------
    for s = 1:32
        P.y     = xY(w,s).y;
        P.ci    = xY(w,s).ci;
        P.ni    = xY(w,s).ni;
        P.ti    = xY(w,s).ti;
        
        yY(s,:) = spm_vec(P)';
    end
    
    
    % subspace over all phonetic parameters
    %----------------------------------------------------------------------
    V   = spm_svd(yY');
    uY  = yY*V;
    
    % clustering for this word
    %----------------------------------------------------------------------
    [pk,pE,pC,pp] = spm_kmeans(uY,K,'random',1);
    
    % save lexical exemplars
    %----------------------------------------------------------------------
    for k = 1:K
        LEX(w,k).pE = spm_unvec(V*pE(k,:)',P);
        LEX(w,k).pC = pC(:,:,k);
        LEX(w,k).pP = spm_inv(pC(:,:,k),norm(pC(:,:,k))/64);
        LEX(w,1).U  = V;
    end
    
end

% articulate exemplars
%--------------------------------------------------------------------------
for w = 1:size(LEX,1)
    for k = 1:size(LEX,2)
        spm_voice_iff(LEX(w,k).pE,FS,1);
    end
end

% likelihood of training set
%--------------------------------------------------------------------------
q     = [];
p     = [];
for i = 1:size(xY,1)
    for j = 1:size(xY,2)
        q(:,end + 1) = spm_voice_likelihood(xY(i,j),LEX);
        p(i,end + 1) = 1;
    end
end

% illustrate classification accuracy
%--------------------------------------------------------------------------
subplot(2,1,2); imagesc(q)
title('classification accuracy')
xlabel('exemplars'), ylabel('words')

% save lexical array
%--------------------------------------------------------------------------
save LEX LEX



