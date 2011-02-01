function pdf = spm_mvNpdf(z,Mu,V)
% Probability Density Function (PDF) of multivariate Normal distribution
% FORMAT pdf = spm_Npdf(z,Mu,V)
%
% z  - ordinates
% Mu - mean (a d-vector)
% V  - d x d variance-covariance matrix
%__________________________________________________________________________
%
% spm_Npdf returns the Probability Density Function (PDF) for the
% multivariate Normal (Gaussian) family of distributions.
%
% The dimension of the Normal distribution is taken as the length of Mu.
% V must be a d x d variance-covariance matrix.
%
% For the univariate Normal distribution (d=1), z can be a matrix of
% arbitrary dimensions - each entry is treated seperately and the PDF
% returned as the corresponding element in a matrix of the same size.
%
% For multivarate PDFs, the ordinates must be in the columns of z, so
% z must have column dimension d. Multiple columns can be entered. 
%
%__________________________________________________________________________
% Copyright (C) 1998-2011 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_mvNpdf.m 4182 2011-02-01 12:29:09Z guillaume $


%-Condition arguments
%--------------------------------------------------------------------------
if nargin<1,   pdf=[]; return, end
if isempty(z), pdf=[]; return, end
if nargin<2,   Mu=0;           end

%-Check Mu, make a column vector, get dimension
%--------------------------------------------------------------------------
if min(size(Mu)) > 1, error('Mu must be a vector'); end
Mu = Mu(:)';
d  = length(Mu);

if nargin<3, V=eye(d); end

%-Size & range checks
%--------------------------------------------------------------------------
if any(any(V~=V')),     error('V must be symmetric'); end
if any(size(V)~=[d,d]), error('V wrong dimension');   end

%-Computation
%--------------------------------------------------------------------------
if d==1
    %-Simpler computation for univariate normal
    %----------------------------------------------------------------------
    pdf = exp(-(z - Mu).^2/(2*V))./sqrt(2*pi*V);
else
    if size(z,1) ~= d, error('z wrong dimension'), end
    z   = z - Mu(:)*ones(1,size(z,2));
    pdf = exp(-0.5*sum((sqrtm(inv(V))*z).^2))/((2*pi)^(d/2)*sqrt(det(V)));
end

return

%-Notes
%==========================================================================
%-The following line computes the PDF in one go for all the ordinates,
% The diag()'s allow for the multiplicity.
% This way is inefficient for large numbers of ordinates.
%
% pdf = ...
%  exp(-0.5*diag((x-Mu)'*inv(V)*(x-Mu))' ) / ( (2*pi)^(d/2) * sqrt(det(V)) );
