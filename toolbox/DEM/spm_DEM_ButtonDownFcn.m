function spm_DEM_ButtonDownFcn
% ButtonDownFcn to play (or save) a movie or sound on button press
% FORMAT spm_DEM_ButtonDownFcn
%
% Requires gcbo to have appropriate UserData; see spm_DEM_movie and
% spm_DEM_play_song
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_DEM_ButtonDownFcn.m 7111 2017-06-16 09:01:09Z guillaume $
 
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
        [filename, pathname] = uiputfile('*.avi','movie file');
        if isequal(filename,0) || isequal(pathname,0), return; end
        fname = fullfile(pathname,filename);
        if spm_check_version('matlab','7.10') > 0
            writerObj = VideoWriter(fname,'Uncompressed AVI');
            writerObj.FrameRate = 15;
            open(writerObj);
            writeVideo(writerObj,S{1});
            close(writerObj);
        else
            movie2avi(S{1},fname,'compression','none','fps',15); %#ok
        end
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
        [filename, pathname] = uiputfile('*.wav','wave file');
        if isequal(filename,0) || isequal(pathname,0), return; end
        fname = fullfile(pathname,filename);
        S{1} = S{1}/max(S{1}(:));
        if spm_check_version('matlab','8.0') > 0
            audiowrite(fname,S{1},S{2},'BitsPerSample',16);
        else
            wavwrite(S{1},S{2},16,fname); %#ok
        end
    end
end
