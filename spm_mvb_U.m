function U = spm_mvb_U(Y,priors,X0,xyz,vox,nu)
% Construct patterns U for Multivariate Bayesian inversion of a linear model
% FORMAT U = spm_mvb_U(Y,priors,X0,xyz,vox,nu)
% Y      - data-feature matrix
% priors - 'null'      % no patterns
%        - 'compact'   % reduced (ns/3); using SVD on local compact support
%        - 'sparse'    % a pattern is a voxel
%        - 'smooth'    % patterns are local Gaussian kernels
%        - 'singular'  % patterns are global singular vectors
%        - 'support'   % the patterns are the images
%
% X0     - confounds
% xyz    - location in mm
% vox    - voxel size in mm
% nu     - number of patterns (for 'compact')
%
% U      - pattern or mode weights
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2007-2022 Wellcome Centre for Human Neuroimaging

 
% defaults
%--------------------------------------------------------------------------
try, X0;  catch, X0  = [];   end
try, xyz; catch, xyz = [];   end
 
% get orders
%--------------------------------------------------------------------------
ns     = size(Y,1);                      % number of samples
nv     = size(Y,2);                      % number of voxels
if nargin < 6
    nu = min(ceil(ns/3),nv);             % number of patterns
end
 
% confounds
%--------------------------------------------------------------------------
if isempty(X0); X0 = zeros(ns,1); end
 
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
        Vlr   = 4/3*pi*dlim^(3/2);            % VOI around a voxel
        Nlr   = round(Vlr/Vvx*.85);           % # of voxel in VOI; keep 85% 
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
                [q,k]  = sort(u(j));
                k      = k(1:nlim);
                j      = j(k);
            end
            U(j,i) = exp(-u(j)/(2*s));
        end
        fprintf('\nThank you\n')
 
    case 'singular'
 
        % get kernel (singular vectors)
        %------------------------------------------------------------------
        Y       = Y - X0*(pinv(X0)*Y);         % remove confounds
        [u,s,v] = spm_svd(Y,1/4);              % c.f., Kaiser criterion
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
 
    case 'compact'
 
        % get kernel (compact vectors)
        %------------------------------------------------------------------
        nc    = max(fix(nv/nu),1);           % voxels in compact support
        X0    = spm_svd(X0);                 % confounds
        Y     = Y - X0*(X0'*Y);              % remove confounds
        C     = sum(Y.^2);                   % variance of Y
        U     = spalloc(nv,nu,nc*nu);
        J     = 1:nv;
        for i = 1:nu
            
            % find maximum variance voxel
            %--------------------------------------------------------------
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
            %--------------------------------------------------------------
            k        = J(j);
            u        = spm_svd(Y(:,k)');
            U(k,i)   = u(:,1);
            
            % remove compact support voxels and start again
            %--------------------------------------------------------------
            J(j)     = [];
            C(j)     = [];
            xyz(:,j) = [];
 
        end
        
    otherwise
        error('Unknown prior')
end
