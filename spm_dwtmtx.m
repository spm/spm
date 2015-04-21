function H = spm_dwtmtx(N,K,T)
% Create basis functions for Discrete (Haar) Wavelet Transform
% FORMAT H = spm_dwtmtx(N,K,T)
%
%   N - dimension
%   K - order: number of basis functions = N/K
%
%   T - option flag for thinning eccentric wavelets [default: false]
%__________________________________________________________________________
%
% spm_dwtmtx creates a matrix for the first few basis functions of a one
% dimensional Haar Discrete Wavelet transform.
%__________________________________________________________________________
% Copyright (C) 2011-2015 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dwtmtx.m 6416 2015-04-21 15:34:10Z guillaume $
 
 
% Create basis set
%==========================================================================
try, T; catch, T = false; end
 
K     = max(K,1); 
H     = ones(N,1);
I     = 1;
while ~isempty(I)
    
    % Find non-zero elements and create new bases
    %----------------------------------------------------------------------
    j      = find(H(:,I(1)) > 0);
    k      = find(H(:,I(1)) < 0);
    nj     = length(j);
    nk     = length(k);
    hj     = sparse(1:fix(nj/2),1,2,nj,1) - 1;
    hk     = sparse(1:fix(nk/2),1,2,nk,1) - 1;
    
    % Add new columns if there are enough elements
    %----------------------------------------------------------------------
    if nj > K
        H(j,end + 1) = hj;
        I(end + 1)   = size(H,2);
    end
    if nk > K
        H(k,end + 1) = hk;
        I(end + 1)   = size(H,2);
    end
 
    % This column has now been expanded
    %----------------------------------------------------------------------
    I(1) = [];
    
end
 
% Thin eccentric basis functions
%==========================================================================
if T
    
    % indices of column (bases) to retain
    %----------------------------------------------------------------------
    j = logical(H(:,1));
    for i = 1:length(j);
        k = find(H(:,i));
        l = length(k);
        if any(k < (N/2 - 2*l)) || any(k > (N/2 + 2*l))
            j(i) = 0;
        end
    end
    
    % thin
    %----------------------------------------------------------------------
    H  = H(:,j);
end
