function pinfo = spm_resss(Vi,Vo,R,flags)
% Create residual sum of squares image (ResSS)
% FORMAT Vo.pinfo = spm_resss(Vi,Vo,R,flags)
% Vi       - vector of mapped image volumes to work on (from spm_vol)
% Vo       - handle structure for mapped output image volume
% R        - residual forming matrix
% Vo.pinfo - plane info matrix of written image
%_______________________________________________________________________
%
% Residuals are computed as R*Y, where Y is the data vector read from
% images mapped as Vi. The residual sum of squares image (mapped as Vo)
% is written.
%
%                           ----------------
%
% For a simple linear model Y = X*B * E, with design matrix X,
% (unknown) parameter vector(s) B, and data matrix Y, the least squares
% estimates of B are given by b = inv(X'*X)*X'*Y. If X is rank
% deficient, then the Moore-Penrose pseudoinverse may be used to obtain
% the least squares parameter estimates with the lowest L2 norm: b =
% pinv(X)*Y.
%
% The fitted values are then y = X*b = X*inv(X'*X)*X'*Y, (or
% y=X*pinv(X)*Y). Since the fitted values y are usually known as
% "y-hat", X*inv(X'*X)*X' is known as the "hat matrix" for this model,
% denoted H.
%
% The residuals for this fit (estimates of E) are e = Y - y.
% Substituting from the above, e = (I-H)*Y, and (I-H), where I is the
% identity matrix (see eye) is the residual forming matrix, denoted R.
%
% For temporally smoothed fMRI models with convolution matrix K, R is a
% little more complicated:
%          K*Y = K*X * B + K*E
% ...a little working shows that H = KX*inv(KX'*KX)*KX'*K (or
% KX*pinv(KX)*K), where KX=K*X, such that the fitted values are y=H*Y.
% The residual forming matrix is then R=(K-H).
%
%                           ----------------
%
% This function can also be used when the b's are images. The residuals
% are then e = Y - X*b, so let Vi refer to the vector of images and
% parameter estimates ([Y;b]), and then R is ([eye(n),-X]), where n is
% the number of Y images.
%
%                           ----------------
%
% Don't forget to either apply any image scaling (grand mean or
% proportional scaling global normalisation) to the image scalefactors,
% or to combine the global scaling factors in the residual forming
% matrix.
%_______________________________________________________________________
% %W% Andrew Holmes, John Ashburner %E%


%-Argument checks
%-----------------------------------------------------------------------
if nargin<4, flags=''; end, if isempty(flags), flags='-'; end
mask = any(flags=='m');
if nargin<3, error('insufficient arguments'); end
ni = size(R,2);					%-ni = #images
if ni~=prod(size(Vi)), error('incompatible dimensions'); end
if Vo.dim(4)~=16, error('only float output images supported'), end

%-Image dimension, orientation and voxel size checks
%-----------------------------------------------------------------------
V = [Vi(:);Vo];
if any(any(diff(cat(1,V.dim),1,1),1)&[1,1,1,0])	%NB: Bombs for single image
	error('images don''t all have the same dimensions'), end
if any(any(any(diff(cat(3,V.mat),1,3),3)))
	error('images don''t all have same orientation & voxel size'), end


%=======================================================================
% - C O M P U T A T I O N
%=======================================================================
Vo.pinfo=[1,0,0]';					%-Set scale & offsets
spm_create_image(Vo);					%-Write image header

Y  = zeros([Vo.dim(1:2),ni]);				%-PlaneStack data

spm_progress_bar('Init',Vo.dim(3),mfilename,['planes / 'num2str(Vo.dim(3))]);

%-Loop over planes computing ResSS
for p=1:Vo.dim(3)
	M = spm_matrix([0 0 p]);			%-Sampling matrix

	%-Read plane data
	for j=1:ni, Y(:,:,j) = spm_slice_vol(Vi(j),M,Vi(j).dim(1:2),0); end

	%-Apply implicit zero mask for lower image types
	if mask, Y(Y(:,:,cat(1,Vi.dim)*[0;0;0;1]<16)==0)=NaN; end

	e  = R*reshape(Y,prod(Vi(1).dim(1:2)),ni)';	%-residuals as DataMtx
	ss = reshape(sum(e.^2,1),Vi(1).dim(1:2));	%-ResSS plane
	spm_write_plane(Vo,ss,p);			%-Write plane
	spm_progress_bar('Set',p);
end

pinfo = Vo.pinfo;					%-Set output argument

%-End
%-----------------------------------------------------------------------
spm_progress_bar('Clear')
