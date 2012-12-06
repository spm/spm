function Vo = spm_resss(Vi,Vo,R,flags)
% Create residual sum of squares image (ResSS)
% FORMAT Vo = spm_resss(Vi,Vo,R,flags)
% Vi          - vector of mapped image volumes to work on (from spm_vol)
% Vo          - handle structure for mapped output image volume
% R           - residual forming matrix
% flags       - 'm' for implicit zero masking
% Vo (output) - handle structure of output image volume after modifications
%               for writing
%
% Note that spm_create_vol needs to be called external to this function -
% the header is not created.
%__________________________________________________________________________
%
% Residuals are computed as R*Y, where Y is the data vector read from
% images mapped as Vi. The residual sum of squares image (mapped as Vo)
% is written.
%
%--------------------------------------------------------------------------
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
% Substituting from the above, e = (I-H)*Y, where I is the identity
% matrix (see eye). (I-H) is called the residual forming matrix,
% denoted R.
%
% Geometrically, R is a projection matrix, projecting the data into the
% subspace orthogonal to the design space.
%
%                           ----------------
%
% For temporally smoothed fMRI models with convolution matrix K, R is a
% little more complicated:
%          K*Y = K*X * B + K*E
%           KY =  KX * B +  KE
% ...a little working shows that hat matrix is H = KX*inv(KX'*KX)*KX'
% (or KX*pinv(KX)), where KX=K*X. The smoothed residuals KE (=K*E) are
% then given from the temporally smoothed data KY (=K*Y) by y=H*KY.
% Thus the residualising matrix for the temporally smoothed residuals
% from the temporally smoothed data is then (I-H).
%
% Usually the image time series is not temporally smoothed, in which
% case the hat and residualising matrices must incorporate the temporal
% smoothing. The hat matrix for the *raw* (unsmoothed) time series Y is
% H*K, and the corresponding residualising matrix is R=(K-H*K).
% In full, that's
%         R = (K - KX*inv(KX'*KX)*KX'*K)
% or      R = (K - KX*pinv(KX)*K)              when using a pseudoinverse
%
%--------------------------------------------------------------------------
%
% This function can also be used when the b's are images. The residuals
% are then e = Y - X*b, so let Vi refer to the vector of images and
% parameter estimates ([Y;b]), and then R is ([eye(n),-X]), where n is
% the number of Y images.
%
%--------------------------------------------------------------------------
%
% Don't forget to either apply any image scaling (grand mean or
% proportional scaling global normalisation) to the image scalefactors,
% or to combine the global scaling factors in the residual forming
% matrix.
%__________________________________________________________________________
% Copyright (C) 1999-2012 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes & John Ashburner
% $Id: spm_resss.m 5097 2012-12-06 16:08:16Z guillaume $

%-Argument checks
%--------------------------------------------------------------------------
if nargin<4, flags=''; end, if isempty(flags), flags='-'; end
mask = any(flags=='m');
if nargin<3, error('insufficient arguments'); end
ni   = size(R,2);                                   %-ni = #images
if ni~=numel(Vi), error('incompatible dimensions'); end
if Vo.dt(1)~=spm_type('float32')
    error('only float output images supported');
end

%-Image dimension, orientation and voxel size checks
%--------------------------------------------------------------------------
spm_check_orientations([Vi Vo]);

%==========================================================================
% - C O M P U T A T I O N
%==========================================================================
Y  = zeros([Vo.dim(1:2),ni]);                       %-PlaneStack data

im = false(ni,1);
for j=1:ni
    im(j)=~spm_type(Vi(j).dt(1),'NaNrep');          %-Images without NaNrep
end

%-Loop over planes computing ResSS
for p=1:Vo.dim(3)

    M = spm_matrix([0 0 p]);                        %-Sampling matrix

    %-Read plane data
    for j=1:ni, Y(:,:,j) = spm_slice_vol(Vi(j),M,Vi(j).dim(1:2),0); end

    %-Apply implicit zero mask for image types without a NaNrep
    if mask, Y(Y(:,:,im)==0)=NaN; end

    e  = R*reshape(Y,prod(Vi(1).dim(1:2)),ni)';     %-residuals as DataMtx
    ss = reshape(sum(e.^2,1),Vi(1).dim(1:2));       %-ResSS plane
    Vo = spm_write_plane(Vo,ss,p);                  %-Write plane
end
