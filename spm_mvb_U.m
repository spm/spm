function U = spm_mvb_U(Y,priors,X0,xyz)
% Construct patterns U for Multivariate Bayesian inversion of a linear model
% FORMAT U = spm_mvb_U(Y,priors,X0,xyz)
% Y      - date feature matrix
% priors - 'null'
%        - 'sparse'
%        - 'smooth'
%        - 'singular'
%        - 'support'
%
% X0     - confounds
% xyz    - location in mm for coherent priors
%__________________________________________________________________________

% defaults
%--------------------------------------------------------------------------
try, X0;  catch, X0  = [];   end
try, xyz; catch, xyz = [];   end

% get orders
%--------------------------------------------------------------------------
ns     = size(Y,1);                 % number of samples
nv     = size(Y,2);                 % unmber of voxels

% confounds
%--------------------------------------------------------------------------
if ~length(X0); X0 = sparse(ns,1); end

% get U: X = Y*P + X0*Q + R
%        P = U*E;           
%--------------------------------------------------------------------------
% assemble emprical priors
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
        s     = 4^2;
        U     = zeros(nv,nv);
        for i = 1:nv
            U(i,:) = exp(-sum((xyz - xyz(:,i)*ones(1,nv)).^2)/(2*s));
        end

    case 'singular'

        % get kernel (singular vectors)
        %------------------------------------------------------------------
        R       = speye(size(X0,1)) - X0*pinv(X0);
        [u s v] = spm_svd((R*Y),1/4);
        U       = v/s;

    case 'support'

        % get kernel (singular vectors)
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

