function [FS,read] = spm_voice_FS(wfile)
% Gets indices or word strings from lexicon 
% FORMAT [FS,read] = spm_voice_FS(wfile)
%
% wfile  - .wav file, audio object or (double) timeseries
%
% FS     - sampling frequency
% read   - function handle: Y = read(wfile);
%
%  This auxiliary routine finds the sampling frequency and returns a
%  function handle appropriate for the sound format in question.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_FS.m 7574 2019-04-19 20:38:15Z karl $

% get timeseries from audio recorder(or from a file)
%--------------------------------------------------------------------------
global VOX

% get source (recorder) and FS
%--------------------------------------------------------------------------
if isa(wfile,'audiorecorder')
    
    FS     = get(wfile,'SampleRate');
    read   = @getaudiodata;
    VOX.FS = FS;
    
elseif isnumeric(wfile)
    
    % timeseries
    %----------------------------------------------------------------------
    try
        FS = get(VOX.audio,'SampleRate');
    catch
        FS = VOX.FS;
    end
    read   = @(Y)Y;
    
else
    
    % sound file
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






