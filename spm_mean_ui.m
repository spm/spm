function spm_average
% promts for a series of images and averages them
% FORMAT spm_average
%_______________________________________________________________________
%
% spm_average simply averages a set of images to produce an mean image that
% is written to "average.img" (in the current directory)
%
% It is assumed that the image and voxels sizes and data format are the same. 
%
%_______________________________________________________________________
% %W% John Ashburner %E%

%-Select images
%----------------------------------------------------------------------------
P = spm_get(Inf,'.img','Select images to be averaged');
%Q = P(1,:); Q = [Q([1:max(find(Q == '/'))]) 'average.img'];
Q = 'average.img';

drawnow;

[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(deblank(P(1,:)));
V=zeros(12,size(P,1));
for i=1:size(P,1)
	V(:,i)=spm_map(P(i,:));
	V(7,i)=V(7,i);
end

%-Compute mean and write headers etc.
%-----------------------------------------------------------------------
scale = spm_mean(V,Q)/size(P,1);

spm_hwrite(Q,DIM,VOX,scale,4,0,ORIGIN,'Mean image');
spm_get_space(Q,spm_get_space(P(1,:)));

for v=V, spm_unmap(v); end;
