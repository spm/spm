function spm_average
% promts for a series of images and averages them
% FORMAT spm_average
%_______________________________________________________________________
%
% spm_average simply averages a set of images to produce an mean image
% that is written to "average.img" (in the current directory).
% Any spatial differences described by the `mat' files are ignored.
%_______________________________________________________________________
% %W% John Ashburner %E%

%-Select images
%----------------------------------------------------------------------------
P = spm_get(Inf,'.img','Select images to be averaged');
%Q = P(1,:); Q = [Q([1:max(find(Q == '/'))]) 'average.img'];
Q = 'test_average.img';

drawnow;

V=spm_vol(P);

%-Compute mean and write headers etc.
%-----------------------------------------------------------------------
scale  = spm_mean(V,Q)/size(P,1);
VOX    = sqrt(sum(V(1).mat(1:3,1:3).^2));
ORIGIN = (V(1).mat\[0 0 0 1]')';
ORIGIN = round(ORIGIN(1:3));
spm_hwrite(Q,V(1).dim(1:3),VOX,scale,4,0,ORIGIN,'Mean image');
spm_get_space(Q,V(1).mat);
