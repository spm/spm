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
% $Id: spm_mvb_U.m 1171 2008-02-27 14:34:38Z christophe $
 
% defaults
%--------------------------------------------------------------------------
try X0;  catch X0  = [];   end %#ok<VUNUS>
try xyz; catch xyz = [];   end %#ok<VUNUS>
 
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
        dlim  = (4*sm)^2;                     % spatial limit (mm)^2
        s     = sm^2;                         % Smoothness variance
        xyz   = xyz';
        nlim  = 256;                          % voxel limit
        Vvx   = prod(vox);                    % volume of a voxel
        Vlr   = 4/3*pi*dlim^(3/2);            % voi around a voxel
        Nlr   = round(Vlr/Vvx*.85);  % estim of # of voxel in voi, keep 85% 
        U     = spalloc(nv,nv,nv*Nlr);        % pre-allocate memory
        fprintf('Creating smooth patterns - please wait\n')
        kk    = floor(nv/4);
        fprintf('\n0%%')
        for i = 1:nv
            if ~rem(i,kk), fprintf('.....%2i%%',25*i/kk); end
            u      = (xyz(:,1) - xyz(i,1)).^2;
            j      = find(u < dlim);
            u(j)   = u(j) + (xyz(j,2) - xyz(i,2)).^2;
            j      = j((u(j) < dlim));
            u(j)   = u(j) + (xyz(j,3) - xyz(i,3)).^2;
            j      = j((u(j) < dlim));
            if length(j)>nlim
                [q k]  = sort(u(j));
                k      = k(1:nlim);
                j      = j(k);
            end
            U(j,i) = exp(-u(j)/(2*s));
        end
        fprintf('\nThank you\n')
 
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
