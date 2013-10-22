function Y = spm_ho_gm(P,M,U,varargin)
% General Gaussian mixture model with derivatives
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
% $Id: spm_ho_gm.m 5709 2013-10-22 11:07:29Z guillaume $


% evaluate
%--------------------------------------------------------------------------
nu  = size(U,2);
if nargin > 3
    
    % evaluate
    %----------------------------------------------------------------------
    for i = 1:nu
        Y(1,i) = 0;
        for j = 1:length(P)
            X      = U(:,i) - P(j).m;
            S      = P(j).c'*P(j).c;
            Y(1,i) = Y(1,i) + exp(P(j).a)*exp(-X'*spm_inv(S)*X/2);
        end
    end
    
    % log transform
    %----------------------------------------------------------------------
    Y = log(Y + eps);
    
else
    
    % evaluate with derivatives dFdu
    %----------------------------------------------------------------------
    for i = 1:nu
        [dFdu,F] = spm_diff('spm_ho_gm',P,M,U(:,i),'eval',3);
        Y(:,i)   = spm_vec(F,dFdu);
        
    end
    
end

% place samples in rows
%----------------------------------------------------------------------
Y  = Y';
