function [S] = spm_voice_read(wfile,LEX)
% Reads and translates a wav file
% FORMAT [S] = DEM_birdsong(file)
%
% wfile  - .wav file
% LEX    - lexical dictionary array
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_read.m 7512 2019-01-05 21:15:16Z karl $


% create lexical structures for subsequent word recognition
%==========================================================================
if nargin < 1, wfile = '../test.wav'; end
if nargin < 2, load LEX,              end

% get FS
%--------------------------------------------------------------------------
try
    [Y,FS] = audioread(wfile,[1,128]);
catch
    [Y,FS] = wavread(wfile,[1,128]);
end

% get data features from a wav file
%--------------------------------------------------------------------------
s0    = 64*FS/1000;                            % smoothing (milliseconds)
I     = 1;
for s = 1:64
    
    % find next initial time point
    %----------------------------------------------------------------------
    if s < 2, Dt = 4; else, Dt = 1; end
    try
        Y = wavread(wfile,[0 FS*Dt] + I(s));
    catch
        break
    end
    
    Y     = spm_conv(abs(Y),s0);
    i     = find(Y/max(Y) > 1/8,1);
    I(s)  = I(s) + i;
    
    % retrieve epoch and decompose at fundamental frequency
    %----------------------------------------------------------------------
    Dt    = 0.8;
    Y     = wavread(wfile,[0 FS*Dt] + I(s));
    xY(s) = spm_voice_ff(Y,FS);
    
    % identify the most likely word and place in structure
    %----------------------------------------------------------------------
    [L,i,str] = spm_voice_likelihood(xY(s),LEX);
    phrase{s} = str;
    PHRASE(s) = LEX(i,1).pE;
    
    % advance sample time
    %----------------------------------------------------------------------
    Dt        = LEX(i,1).pE.ti;
    I(s + 1)  = I(s) + round(FS*Dt);
    
end

% articulate
%--------------------------------------------------------------------------
spm_voice_iff(PHRASE,FS);

% illustrate classification accuracy
%--------------------------------------------------------------------------
disp(phrase)




