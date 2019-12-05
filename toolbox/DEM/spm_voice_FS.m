function [FS,read] = spm_voice_FS(wfile)
% Sampling frequency and function handle for handling sound signals
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
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_FS.m 7750 2019-12-05 17:54:29Z spm $


% get timeseries from audio recorder(or from a file)
%--------------------------------------------------------------------------
global VOX

% default
%--------------------------------------------------------------------------
if ~nargin
    try
        FS = get(VOX.audio,'SampleRate');
    catch
        FS = 22050;
    end
    return
end

% get source (recorder) and FS
%--------------------------------------------------------------------------
if isa(wfile,'audiorecorder')
    
    FS     = get(wfile,'SampleRate');
    read   = @getaudiodata;

elseif isnumeric(wfile)
    
    % timeseries
    %----------------------------------------------------------------------
    try
        FS = get(VOX.audio,'SampleRate');
    catch
        try
            FS = VOX.FS;
        catch
            FS = 22050;
        end
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
        [~,FS] = wavread(wfile,[1 1]);
        read   = @wavread;
    end
    
end

% place sampling frequency in global VOX structure
%----------------------------------------------------------------------
% VOX.FS = FS;
