function [M] = spm_DEM_movie(qU,S,FPS);
% displays a movie and set ButtonDownFunction to play it
% FORMAT [M] = spm_DEM_movie(qU,S,FPS);
%
% qU   - conditional moments of states (see spm_DEM) or v
% S    - .mat file or structure containing
%        S.V   - image modes (V) 
%        S.F   - image template (format, for spm_unvec)
%
% M    - movie array
% FPS  - Frames per second (Hz)
%
% A button press on the image will play the movie. The i-th frame is simply S.V*qU.v{1}(:,i)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_DEM_movie.m 3715 2010-02-08 13:57:26Z karl $
 
% load image modes
%--------------------------------------------------------------------------
try
    FPS;
catch
    FPS = 32;
end
try
    S.F;
catch
    try
        S = load(file);
    catch
        S = load('DEM_IMG');
    end
end
try
    U = qU.v{1};
catch
    U = qU;
end
 
% get parameters
%==========================================================================
Nc    = size(S.V,1);                % number of channels
Nm    = size(S.V,2);                % number of modes
Nf    = size(U,2);                  % number of frames
 
 
% reconstitute movie
%--------------------------------------------------------------------------
sw = warning('off','all');
for i = 1:Nf
    Y    = spm_unvec(S.V(:,1:Nm)*U(1:Nm,i),S.F);
    Y    = Y - min(Y(:));
    Y    = 255*Y/max(Y(:));
    M(i) = im2frame(uint8(Y));
end
warning(sw);
 
% Graphics
%==========================================================================
set(gca,'units','pixels')
position    = get(gca,'position');
position(3) = size(S.F,2)*.8;
position(4) = size(S.F,1)*.8;
set(gca,'position',position);
set(gca,'units','normalized');
imagesc(M(1).cdata);
axis off
 
% set ButtonDownFcn
%--------------------------------------------------------------------------
h = get(gca,'Children');
set(h(1),'Userdata',{M,FPS})
set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')
