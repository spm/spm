% Utility routine to load images and create basis functions using a
% discrete cosine basis set (over a feature dimension). This is written
% specifically for the images used in this demonstration and should be
% tailored for any new images.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_get_faces.m 1703 2008-05-21 13:59:23Z karl $
 
% try to read all images in current directory
%--------------------------------------------------------------------------
clear
 
file  = dir;
IMG   = [];
for i = 1:length(file)
    try
        if strcmp(file(i).name(1:4),IMG(end).name(1:4))
            IMG(end + 1).y = imread(file(i).name);
            IMG(end).name  = file(i).name;
        end
    catch
        try
            IMG(end + 1).y = imread(file(i).name);
            IMG(end).name  = file(i).name;
        end
    end
end
 
% sort images over feature space
%--------------------------------------------------------------------------
N     = length(IMG);
for i = 1:N
    j    = [find(IMG(i).name =='_') + 1:find(IMG(i).name =='%') - 1];
    E(i) = eval(IMG(i).name(j));
end
[i j] = sort(E);
IMG   = IMG(j);
 
% show images
%--------------------------------------------------------------------------
for i = 1:N
    M(i) = im2frame(IMG(i).y);
end
movie(M,1,12);
 
% basis set; DCT over feature space v =[0,1]
%--------------------------------------------------------------------------
Nm    = 3;
dv    = 4;
F     = double(IMG(1).y(1:dv:end,1:dv:end,:));
for i = 1:N
    Y(:,i) = spm_vec(double(IMG(i).y(1:dv:end,1:dv:end,:)));
end
v     = linspace(0,1,N);
for i = 1:Nm
    U(:,i) = cos((i - 1)*pi*v');
end
V     = Y*U;
 
 
% Save DEM_IMG .mat file
%--------------------------------------------------------------------------
save DEM_IMG V F U




