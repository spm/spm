function [Linv,W]=MNestimator(L,C,SNR2,trunc)
% function [Linv,W]=MNestimator(L,C,SNR2,trunc)
% This function makes a basic minimum-morm pseudoinverse operator = MN estimator.
% - L: the forward solution, lead-field matrix
% - C: noise covariance matrix
%       This is used for combining different sensortypes, whitening the
%       data and setting the regularisation parameter. If C is empty, it is
%       assumed eye matrix --- that would be the basic Tikhonov 0
%       regularization for one sensortype.
% - SNR2: the assumed ratio of variances of signal and noise, used for
%       setting the regularisation parameter. 
% - trunc: the number of (smallest) singular values of C that are set to 
%       zero before making the whitener. For example, if the data / C has
%       been SSP-projected, "trunc" needs to be at least the number of
%       components projected away.
%
% The whitening/regularization approach of this routine follows that used in
% the MNE software, 
% http://martinos.org/mne/, 
% http://www.martinos.org/meg/manuals/MNE-manual-2.7.pdf
% --- the influence of MNE software on this function is fully acknowledged!
% 
% Copyright (c) 2011--2013 Matti Stenroos (matti.stenroos@aalto.fi)
%  -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
%                !! There is no warranty of any kind !!
% -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
% 17 Oct 2013
Nsens=size(L,1);

if nargin==3 || isempty(trunc)
    trunc=0;
end
if isempty(C)
    C=eye(Nsens);
end

W=MakeWhitener(C,trunc);

LW=W*L;
lambda2=trace(LW*LW')/(Nsens*SNR2);
temp=LW*LW'+lambda2*eye(Nsens);
Linv=LW'/temp*W;

function [W,deW]=MakeWhitener(C,trunc)
% function [W,deW]=MakeWhitener(C,trunc)
% Computes whitener C^(-1/2) and de-whitener sqrt(C)
% trunc: how many smallest eigenvalues are set to zero before inversion
% (all other regularization has to be applied to C before calling this
% routine)
%
% Copyright (c) 2011--2013 Matti Stenroos (matti.stenroos@aalto.fi)
% -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
%                !! There is no warranty of any kind !!
% -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
% version 130610
if min(size(C))==1,
    C=diag(diag(C));
end

[U,D]=eig(C);
[ds,sortind]=sort(diag(D));
if trunc
    ds(1:trunc)=0;
end
Us=U(:,sortind);
dsinv=zeros(size(ds));
dsinv((trunc+1):end)=1./ds((trunc+1):end);
W=Us*diag(sqrt(dsinv))*Us';

if nargout==2
    deW=U*diag(sqrt(diag(D)))*Us';
end

