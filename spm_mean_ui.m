function spm_average
% promts for a series of images and averages them
% FORMAT spm_average
%____________________________________________________________________________
%
% spm_average simply averages a set of images to produce an mean image that
% is written to "average.img" (in the same subdirectory)
%
% It is assumed that the image and voxels sizes and data format are the same. 
%
%__________________________________________________________________________
% %W% %E%

%----------------------------------------------------------------------------
P    = spm_get(Inf,'.img','Select images to be averaged');
Q    = P(1,:);
Q    = [Q([1:max(find(Q == '/'))]) 'average.img'];

[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(P(1,:));
spm_hwrite(Q,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP);

spm_mean(prod(DIM),TYPE,Q,P)
spm_get_space(Q,spm_get_space(P(1,:)));

