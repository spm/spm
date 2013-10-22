function Y = spm_ho_poly(P,M,U,varargin)
% General polynomial mapping with derivatives
% FORMAT Y = spm_ho_poly(P,M,U)
%
% P    - polynomial parameters (P{i} = i-th order coefficients)
% M    - model structure
% U    - (m,n) inputs
%
% Y(i) =  P{1} + P{2}*U(:,i) + P{3}*kron(U(:,i),U(:,i)) + ...
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ho_poly.m 5709 2013-10-22 11:07:29Z guillaume $


% evaluate
%--------------------------------------------------------------------------
nu  = size(U,2);
if nargin > 3
    
    % evaluate
    %----------------------------------------------------------------------
    Y     = sparse(1,nu) + P{1};
    for i = 1:nu
        X     = 1;
        for j = 2:length(P)
            X      = kron(X,U(:,i));
            Y(1,i) = Y(1,i) + spm_vec(P{j})'*X;
        end
    end
    
else
    
    % evaluate with derivatives dFdu
    %----------------------------------------------------------------------
    for i = 1:nu
        [dFdu,F] = spm_diff('spm_ho_poly',P,M,U(:,i),'no diff',3);
        Y(:,i)   = spm_vec(F,dFdu);
    end
    
end

% place samples in rows
%----------------------------------------------------------------------
Y  = Y';
