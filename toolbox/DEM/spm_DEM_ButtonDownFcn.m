function spm_DEM_ButtonDownFcn
% ButtonDownFcn to play a movie or sound on button press
% FORMAT spm_DEM_ButtonDownFcn
%
% Require gcbo to have appropriate UserData; see spm_DEM_movie and
% spm_DEM_play_song
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_DEM_ButtonDownFcn.m 1703 2008-05-21 13:59:23Z karl $
 
% default
%--------------------------------------------------------------------------
S = get(gcbo,'Userdata');
if isstruct(S{1})
    
    % play movie
    %----------------------------------------------------------------------
    movie(S{1},1,S{2});
    return
    
    % save avi file
    %----------------------------------------------------------------------
    [FILENAME, PATHNAME] = uiputfile('*.wav','wave file');
    NAME = fullfile(PATHNAME,FILENAME);
    movie2avi(S{1},NAME,'compression','none')
    
else
    
    % play sound
    %----------------------------------------------------------------------
    soundsc(S{1},S{2});
    %return
    
    % save wav file
    %----------------------------------------------------------------------
    [FILENAME, PATHNAME] = uiputfile('*.wav','wave file');
    NAME = fullfile(PATHNAME,FILENAME);
    wavwrite(S{1},S{2},32,NAME)
end
 

