function [J] = spm_dcm_J(Y,U,X0,dt,R)
% VOI extraction of adjusted data and Markov Blanket decomposition
% FORMAT [J] = spm_dcm_J(Y,U,X0,dt,R)
%
% Y  - response variable
% U  - exogenous inputs
% XO - confounds
% dt - time bin
% R  - element wise retriction matrix
%
%__________________________________________________________________________
% This routine evaluates the effective connectivity of a dynamic causal
% model based upon the Jacobian (i.e., state matrix) of a stochastic
% differential equation. In other words, it approximates the coupling among
% hidden states to first order, under some simplifying assumptions.
% Starting from a linear state space model, in which the response variable
% (y) is a linear convolution (K) of some hidden states (x) subject to
% observation and system noise (r and e) respectively, we have:
%
% D*x = x*J' + e  => K*D*x = K*x*J' + K*e = D*y = y*J' + K*e + D*r - r*J'
% y   = K*x  + r  =>   D*y = K*D*x  + D*r
%
% This means we can approximate the system with a general linear model:
%
% D*y = y*J' + w:   cov(w) = h(1)*K*K' + h(2)*D*D' + h(3)*I
%
% Where, h(3)*I  = h(2)*J*J', approximately; noting that the leading
% diagonal of J will dominate (and be negative). If K is specified in terms
% of convolution kernels, then the covariance components of the linearised
% system can be expressed as:
%
% K = k(1)*K{1} + k(2)*K{2} + ...
%   => K*K' = k(1)*k(1)*K{1}*K{1}' + k(1)*k(2)*K{1}*K{2}' ...
%
% Where k(i)*k(j) replaces the [hyper]parameter h(1) above. This linearized
% system can be solved using parametric empirical Bayes (PEB) for each
% response variable, under the simplifying assumption that there are the
% same number of observations and hidden states.
%
% This allows large graphs to be inverted by considering the afferents
% (i.e., influences on) to each node sequentially. Redundant elements of
% the Jacobian (i.e., connections) are subsequently removed using Bayesian
% model reduction (BMR). The result is a sparse Jacobian that corresponds
% to the coupling among hidden states that generate observed double
% responses, to first-order.
%
% See: Frässle S, Lomakina EI, Kasper L, Manjaly ZM, Leff A, Pruessmann KP,
% Buhmann JM, Stephan KE. A generative model of whole-brain effective
% connectivity.Neuroimage. 2018 Oct 1;179:505-529.
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_J.m 7679 2019-10-24 15:54:07Z spm $

% first order Jacobian
%==========================================================================

% convolution operators
%--------------------------------------------------------------------------
[nt,nu]    = size(U);
[nt,nv]    = size(Y);
xBF.dt     = dt;
xBF.name   = 'hrf';
xBF.length = 32;
xBF        = spm_get_bf(xBF);
K{1}       = spm_convmtx(xBF.bf,nt,'square');
K{2}       = spm_convmtx([1; 0],nt,'square');
K{3}       = spm_convmtx([1;-1],nt,'square');

if nargin < 5
    R      = ones(nv,nv);
end

% covariance components
%--------------------------------------------------------------------------
for k = 1:numel(K)
    C{k} = K{k}*K{k}';
    C{k} = nt*C{k}/trace(C{k})/64;
end

% remove confounds
%--------------------------------------------------------------------------
Y       = Y - X0*pinv(X0)*Y;
Y       = Y/std(Y(:));
U       = U/std(U(:));
J       = zeros(nv + nu,nv + nu);
for v   = 1:nv
    
    % Parametric (empirical) Bayes
    %----------------------------------------------------------------------
    dY         = gradient(Y(:,v),dt);
    pE{v}      = zeros(nv + nu,1);
    pC{v}      = speye(nv + nu,nv + nu);
    
    pE{v}(v)   = -1;
    pC{v}(v,v) = 1/512;
    r          = find( R(:,v));
    s          = find(~R(:,v));
    u          = [r; (1:nu)' + nv];
    pE{v}(s)   = 0;
    pC{v}(s,s) = 0;
    
    P{1}.X     = [Y(:,r) U];
    P{1}.C     = C;
    P{2}.X     = pE{v}(u);
    P{2}.C     = pC{v}(u,u);
    PEB        = spm_PEB(dY,P,512);
    Ep{v}      = pE{v};
    Cp{v}      = pC{v};
    Ep{v}(u)   = PEB{2}.E;
    Cp{v}(u,u) = PEB{2}.C;
    
end

% Bayesian Model Reduction
%==========================================================================
OPT = 'symmetric';
switch OPT
    
    case 'none'
        
        % no reduction
        %------------------------------------------------------------------
        for v   = 1:nv
            
            % Jacobian for this voxel
            %--------------------------------------------------------------
            J(:,v)   = Ep{v};
            
        end
        
    case 'symmetric'
        
        % reciprocal coupling shrinkage priors
        %------------------------------------------------------------------
        for i = 1:nv
            for j = (i + 1):nv
                
                if R(i,j)
                    
                    % afferent connection
                    %----------------------------------------------------------
                    rCi          = pC{i};
                    rCi(j,j)     = 0;
                    [Fi,sEi,sCi] = spm_log_evidence(Ep{i},Cp{i},pE{i},pC{i},pE{i},rCi);
                    
                    % efferent connection
                    %----------------------------------------------------------
                    rCj          = pC{j};
                    rCj(i,i)     = 0;
                    [Fj,sEj,sCj] = spm_log_evidence(Ep{j},Cp{j},pE{j},pC{j},pE{j},rCj);
                    
                    
                    % accept new priors if F has increased
                    %----------------------------------------------------------
                    if (Fi + Fj) > 3
                        Ep{i}  = sEi;
                        Ep{j}  = sEj;
                        Cp{i}  = sCi;
                        Cp{j}  = sCj;
                        pC{i}  = rCi;
                        pC{j}  = rCj;
                    end
                end
                
            end
            fprintf('evaluating (%i) of %i states\n',i,nv)
        end
        
        % assemble Jacobian
        %------------------------------------------------------------------
        for v   = 1:nv
            J(:,v)   = Ep{v};
        end
        
        
    case 'all'
        
        % Jacobian for this voxel
        %------------------------------------------------------------------
        for v   = 1:nv
            
            % Bayesian Model Reduction
            %--------------------------------------------------------------
            BMR.M.pE = pE;
            BMR.M.pC = pC;
            BMR.Ep   = Ep;
            BMR.Cp   = Cp;
            BMR      = spm_dcm_bmr_all(BMR,'all','BMS');
            
            % Jacobian for this voxel
            %--------------------------------------------------------------
            J(:,v)   = BMR.Ep;
            
        end
        
    otherwise
        
end

% retain coupling between states
%--------------------------------------------------------------------------
J  = J(1:nv,1:nv);

