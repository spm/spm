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
        dlim  = 8;                           % spatial extent fixed at 2*sm
        s     = sm^2;                                 % Smoothness variance
        xyz   = xyz';
        Vvx   = prod(vox);                              % volume of a voxel
        Vlr   = 4/3*pi*(2*dlim*s)^(3/2);               % voi around a voxel
        Nlr   = round(Vlr/Vvx*.9);   % estim of # of voxel in voi, keep 90% 
        U     = spalloc(nv,nv,nv*Nlr);                % pre-allocate memory
        unit  = ones(nv,1);
        kk    = floor(nv/4);
        fprintf('\n0%%')
        for i = 1:nv
            if ~rem(i,kk), fprintf('.....%2i%%',25*i/kk); end
            u      = exp(-sum((xyz(1:i,:) - ones(i,1)*xyz(i,:)).^2,2)/(2*s));
            tmp    = zeros(nv,1); tmp(1:i) = u.*(u > exp(-dlim));
            U(:,i) = sparse(tmp);
        end
        U = U + U' - speye(nv);
        fprintf('\n')
 
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
