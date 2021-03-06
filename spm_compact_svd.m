function U = spm_compact_svd(Y,xyz,nu)
% Local SVD with compact support for large matrices
% FORMAT U = spm_compact_svd(Y,xyz,nu)
% Y     - matrix
% xyz   - location
% nu    - number of vectors
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2010-2022 Wellcome Centre for Human Neuroimaging


% get orders
%--------------------------------------------------------------------------
ns     = size(Y,1);                % number of samples
nv     = size(Y,2);                % number of components (voxels/vertices)

% get kernel (compact vectors)
%--------------------------------------------------------------------------
nc    = max(fix(nv/nu),1);         % voxels in compact support
C     = sum(Y.^2);                 % variance of Y
U     = spalloc(nv,nu,nc*nu);
J     = 1:nv;
for i = 1:nu
    
    % find maximum variance voxel
    %----------------------------------------------------------------------
    [v,j] = max(C);
    d     = 0;
    for k = 1:size(xyz,1)
        d  = d + (xyz(k,:) - xyz(k,j)).^2;
    end
    [d,j] = sort(d);
    try
        j = j(1:nc);
    end
    
    % save principal eigenvector
    %----------------------------------------------------------------------
    k        = J(j);
    u        = spm_svd(Y(:,k)');
    U(k,i)   = u(:,1);
    
    % remove compact support voxels and start again
    %----------------------------------------------------------------------
    J(j)     = [];
    C(j)     = [];
    xyz(:,j) = [];
    
end
