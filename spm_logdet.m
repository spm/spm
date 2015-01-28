function H = spm_logdet(C)
% Compute the log of the determinant of positive (semi-)definite matrix C
% FORMAT H = spm_logdet(C)
% H = log(det(C))
%
% spm_logdet is a computationally efficient operator that can deal with
% full or sparse matrices. For non-positive definite cases, the determinant
% is considered to be the product of the positive singular values.
%__________________________________________________________________________
% Copyright (C) 2008-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston and Ged Ridgway
% $Id: spm_logdet.m 6321 2015-01-28 14:40:44Z karl $

% Note that whether sparse or full, rank deficient cases are handled in the
% same way as in spm_logdet revision 4068, using svd on a full version of C


% remove null variances
%--------------------------------------------------------------------------
i       = find(diag(C));
C       = C(i,i);
[i,j,s] = find(C);
if any(isnan(s)), H = nan; return; end

% TOL = max(size(C)) * eps(max(s)); % as in MATLAB's rank function
%--------------------------------------------------------------------------
TOL   = 1e-16;

if any(i ~= j)
    
    % assymetric matrix
    %------------------------------------------------------------------
    if norm(spm_vec(C - C'),inf) > TOL
        
        s = svd(full(C));
        
    else
        
        % non-diagonal sparse matrix
        %------------------------------------------------------------------
        if issparse(C)
            
            % Note p is unused but requesting it can make L sparser
            %--------------------------------------------------------------
            [L,nondef,p] = chol(C, 'lower', 'vector');
            if ~nondef
                
                % pos. def. with Cholesky decomp L, and det(C) = det(L)^2
                %----------------------------------------------------------
                H = 2*sum(log(full(diag(L))));
                return
                
            end
            s = svd(full(C));
            
        else
            
            % non-diagonal full matrix
            %--------------------------------------------------------------
            try
                R = chol(C);
                H = 2*sum(log(diag(R)));
                return
            catch
                s = svd(C);
            end
        end
    end
end

% if still here, singular values in s (diagonal values as a special case)
%--------------------------------------------------------------------------
H     = sum(log(s(s > TOL & s < 1/TOL)));
