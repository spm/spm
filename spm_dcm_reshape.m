function [A,B,C,H,varargout] = spm_dcm_reshape(P,m,n,r)
% converts free parameter vector to matrices
% FORMAT:
%    for bilinear DCM:  [A B C H]   = spm_dcm_reshape(P,m,n,[r]);
%    for nonlinear DCM: [A B C H D] = spm_dcm_reshape(P,m,n,[r]); 
% FORMAT 
% P     - parameter vector
% m     - number of inputs
% n     - number of regions
% [r]   - returns relative connections {without scaling by P(1)}
%
% A     - intrinsic connections
% B     - modulatory connections
% C     - direct connections
% H     - hemodynamic parameters
% D     - nonlinear connections
% 
% NB: order of parameter vector is: 
% bilinear neural params -> hemodynamic params -> nonlinear neural params
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_reshape.m 2504 2008-11-29 15:53:11Z klaas $
 
 
% scale intrinsic connections {A}
%--------------------------------------------------------------------------
P     = full(P);
if nargin == 4
    q = 1;
else
    q = exp(P(1));
end
P(1)  = [];
 
% fill in intrinsic connections {A}
%--------------------------------------------------------------------------
j     = 1:n*n;
A     = reshape(P(j),n,n)*q;
P(j)  = [];
 
% fill in modulatory connections {B}
%--------------------------------------------------------------------------
j     = 1:n*n*m;
B     = reshape(P(j),n,n,m)*q;
P(j)  = [];
 
% fill in direct connections {C}
%--------------------------------------------------------------------------
j     = 1:n*m;
C     = reshape(P(j),n,m);
P(j)  = [];

% fill in 6 hemodynamic parameters {H}
%--------------------------------------------------------------------------
j     = 1:n*6;
H     = reshape(P(j),n,6);
P(j)  = [];

if nargout>4
    % fill in nonlinear connections {D}
    % NB: there is one D matrix for each region, and each D matrix is
    %     normalized in the same way as the A and B matrices
    %--------------------------------------------------------------------------
    j            = 1:n*n*n;
    varargout{1} = reshape(P(j),n,n,n)*q;
    P(j)         = [];
end

return
