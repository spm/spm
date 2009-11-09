function [A,B,C,H,D] = spm_dcm_reshape(P,m,n,r)
% Converts free parameter vector to matrices
% FORMAT [A,B,C,H,[D]] = spm_dcm_reshape(P,m,n,[r])
% 
% P     - parameter vector
% m     - number of inputs
% n     - number of regions
% [r]   - returns relative connections {without scaling by P(1)}
%
% A     - intrinsic connections
% B     - modulatory connections
% C     - direct connections
% H     - hemodynamic parameters
% [D]   - nonlinear connections
% 
% NB: order of parameter vector is: 
% bilinear neural params -> hemodynamic params -> nonlinear neural params
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_reshape.m 3547 2009-11-09 18:29:59Z guillaume $
 
 P    = full(P);
 i    = cumsum([1 n*n n*n*m n*m n*6 n*n*n]);
 
% scaling factor (relative connections) [1]
%--------------------------------------------------------------------------
if nargin == 4
    q = 1;
else
    q = exp(P(1));
end
 
% fill in intrinsic connections {A} [n x n]
%--------------------------------------------------------------------------
A     = reshape(P((i(1)+1):i(2)),n,n) * q;
 
% fill in modulatory connections {B} [n x n x m]
%--------------------------------------------------------------------------
B     = reshape(P((i(2)+1):i(3)),n,n,m) * q;
 
% fill in direct connections {C} [n x m]
%--------------------------------------------------------------------------
C     = reshape(P((i(3)+1):i(4)),n,m);

% fill in 6 hemodynamic parameters {H} [n x 6]
%--------------------------------------------------------------------------
H     = reshape(P((i(4)+1):i(5)),n,6);

% fill in nonlinear connections {D} [n x n x n]
%--------------------------------------------------------------------------
if nargout > 4
    % NB: there is one D matrix for each region, and each D matrix is
    %     normalized in the same way as the A and B matrices
    D = reshape(P((i(5)+1):i(6)),n,n,n) * q;
end
