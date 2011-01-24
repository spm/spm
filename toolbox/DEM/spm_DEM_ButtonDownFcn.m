function spm_DEM_ButtonDownFcn
% ButtonDownFcn to play (or save) a movie or sound on button press
% FORMAT spm_DEM_ButtonDownFcn
%
% Requires gcbo to have appropriate UserData; see spm_DEM_movie and
% spm_DEM_play_song
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_DEM_ButtonDownFcn.m 4170 2011-01-24 18:37:42Z karl $
 
% default
%--------------------------------------------------------------------------
S = get(gcbo,'Userdata');
if isstruct(S{1})
    
    % play movie
    %----------------------------------------------------------------------
    movie(S{1},1,S{2});
    
    if strcmp(get(gcf,'SelectionType'),'normal')
        return
    else
        % save avi file
        %------------------------------------------------------------------
        [FILENAME, PATHNAME] = uiputfile('*.avi','movie file');
        NAME = fullfile(PATHNAME,FILENAME);
        movie2avi(S{1},NAME,'compression','none')
    end
    
else
    
    % play sound
    %----------------------------------------------------------------------
    soundsc(S{1},S{2});
    
    if strcmp(get(gcf,'SelectionType'),'normal')
        return
    else
        
        % save wav file
        %------------------------------------------------------------------
        [FILENAME, PATHNAME] = uiputfile('*.wav','wave file');
        NAME = fullfile(PATHNAME,FILENAME);
        S{1} = S{1}/max(S{1}(:));
        wavwrite(S{1},S{2},16,NAME)
    end
end
 

