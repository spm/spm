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

%-Get scalefactors and compute weights
%-----------------------------------------------------------------------
scales = zeros(size(P,1),1);
for i=1:size(P,1)
	[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = ...
					spm_hread(deblank(P(i,:)));
	scales(i) = SCALE;
	if (i==1)
		DIM1 = DIM;
		TYPE1 = TYPE;
	end
	if DIM ~= DIM1 | TYPE ~= TYPE1 | OFFSET ~= 0
		error(['spm_mean can not cope with image "' deblank(P(i,:)) '".']);
	end
end

%-Compute mean and write headers etc.
%-----------------------------------------------------------------------
scale = spm_mean(prod(DIM),TYPE,Q,P,scales/length(scales));
spm_hwrite(Q,DIM,VOX,scale,TYPE,OFFSET,ORIGIN,DESCRIP);
spm_get_space(Q,spm_get_space(P(1,:)));

