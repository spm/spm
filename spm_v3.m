function [I] = spm_v3(X,V3);
% creates an image of estimated |activity| based on 37 channel data
% FORMAT [I] = spm_v3(X,V3);
% I	-	{n,m} Image
% X	-	data {37 channels}
% V3	-	{n*m,37} matrix of coeficients
%
% spm_v3 uses a linear operation to produce an estimate of the distribution
% of activity (in terms of the vector length) over the surface of the
% brain.  The complex transformation matrix (V3) should be a global variable.
%
% i.e.		I  = V3.X
%
%__________________________________________________________________________
% %W% %E%

%---------------------------------------------------------------------------
I = reshape(abs(V3*X(:)),17,17);
