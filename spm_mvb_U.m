function U = spm_mvb_U(Y,priors,X0,xyz,vox)
% Constructs patterns U for Multivariate Bayesian inversion of a linear model
% FORMAT U = spm_mvb_U(Y,priors,X0,xyz,vox)
% Y      - date filature matrix
% priors - 'null'
%        - 'sparse'
%        - 'smooth'
%        - 'singular'
%        - 'support'
%
% X0     - confounds
% xyz    - location in mm for coherent priors
% vox    - voxel size for coherent priors
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_mvb_U.m 1161 2008-02-22 12:18:05Z karl $
 
% defaults
%--------------------------------------------------------------------------
try, X0;  catch, X0  = [];   end
try, xyz; catch, xyz = [];   end
 
% get orders
%--------------------------------------------------------------------------
ns     = size(Y,1);                 % number of samples
nv     = size(Y,2);                 % number of voxels
 
% confounds
%--------------------------------------------------------------------------
if ~length(X0); X0 = sparse(ns,1); end
 
% get U: X = Y*P + X0*Q + R
%        P = U*E;           
%--------------------------------------------------------------------------
% assemble empirical priors
%==========================================================================
switch priors
 
    case 'null'
        %------------------------------------------------------------------
        U     = sparse(nv,0);
 
    case 'sparse'
        %------------------------------------------------------------------
        U     = speye(nv,nv);
 
    case 'smooth'
        %------------------------------------------------------------------
        sm    = 4;                            % smoothness fixed at 4mm std
        dlim  = 16^2;                         % spatial limit (mm)^2
        s     = sm^2;                         % Smoothness variance
        xyz   = xyz';
        nlim  = 256;                          % voxel limit
        iv    = 1:4:nv;
        nu    = length(iv);
        U     = spalloc(nv,nu,nu*nlim);       % pre-allocate memory
        l     = 1;
        fprintf('Creating smooth patterns - please wait\n')
        for i = iv
            u      =        (xyz(:,1) - xyz(i,1)).^2;
            j      =   find(u    < dlim);
            u(j)   = u(j) + (xyz(j,2) - xyz(i,2)).^2;
            j      = j(find(u(j) < dlim));
            u(j)   = u(j) + (xyz(j,3) - xyz(i,3)).^2;
            j      = j(find(u(j) < dlim));
            [q k]  = sort(u(j));
            try,k  = k(1:nlim); end
            j      = j(k);
            U(j,l) = exp(-u(j)/(2*s));
            l      = l + 1;
        end
        fprintf('Thank you\n')
 
    case 'singular'
 
        % get kernel (singular vectors)
        %------------------------------------------------------------------
        R       = speye(size(X0,1)) - X0*pinv(X0);
        [u s v] = spm_svd((R*Y),1/4);              % c.f., Kaiser criterion
        U       = v/s;
 
    case 'support'
 
        % get kernel (image vectors)
        %------------------------------------------------------------------
        if nv > ns
            R  = speye(size(X0,1)) - X0*pinv(X0);
            U  = (R*Y)';
        else
            U  = speye(nv,nv);
        end
 
    otherwise
        disp('unknown prior')
        return
end
