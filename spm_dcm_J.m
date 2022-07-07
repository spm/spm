function [J,K] = spm_dcm_J(Y,U,X0,dt,R,I,D)
% VOI extraction of adjusted data and Markov Blanket decomposition
% FORMAT [J,K] = spm_dcm_J(Y,U,X0,dt,R,I,D)
%
% Y  - response variable
% U  - exogenous inputs
% XO - confounds
% dt - time bin
% R  - distance matrix
% I  - self inhibition [default: -1]
% D  - upper bound on distance (mm) [default: 64]
%
% J(nv,nv)  - Jacobian
% K(nv,nu)  - input block
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
% See: Frassle S, Lomakina EI, Kasper L, Manjaly ZM, Leff A, Pruessmann KP,
% Buhmann JM, Stephan KE. A generative model of whole-brain effective
% connectivity.Neuroimage. 2018 Oct 1;179:505-529.
%
% GRAPHICAL OUTPUT Sparse connectivity: the figure illustrates the sparsity
% of effective connectivity using Bayesian model reduction. The left panel
% shows the log evidence for a series of models that preclude connections
% beyond a certain distance or radius. This log evidence is been normalised
% to the log evidence of the model with the least marginal likelihood. The
% middle panel shows the ensuing sparse coupling (within the upper bound of
% D mm) as an adjacency matrix, where particles have been ordered using a
% nearest neighbour scheme in voxel space. The blue dots indicate
% connections that have been removed by Bayesian model reduction. The right
% panel zooms in on the first 32 particles, to show local connections that
% were retained (red) or removed (blue).
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% preliminaries
%--------------------------------------------------------------------------
if nargin < 7, D = 32;          end             % maximum length  (mm)
if nargin < 6, I = 1;           end             % self inhibition (Hz)

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
K{4}       = spm_convmtx(gradient(xBF.bf),nt,'square');
K{5}       = spm_convmtx(gradient(gradient(xBF.bf)),nt,'square');

% covariance components
%--------------------------------------------------------------------------
for k = 1:numel(K)
    C{k} = K{k}*K{k}';
    C{k} = nt*C{k}/trace(C{k})/32;
end

% remove confounds
%--------------------------------------------------------------------------
Y       = Y - X0*pinv(X0)*Y;
Y       = Y/std(Y(:));
U       = U/std(U(:));
I       = -abs(I);
for v   = 1:nv
    
    % Parametric (empirical) Bayes
    %----------------------------------------------------------------------
    dY         = gradient(Y(:,v),dt);
    pE{v}      = zeros(nv + nu,1);
    pC{v}      = speye(nv + nu,nv + nu);
    
    pE{v}(v)   = I;
    pC{v}(v,v) = 1/512;
    d          = R(:,v) < D;
    r          = find( d);
    s          = find(~d);
    u          = [r; (1:nu)' + nv];
    pE{v}(s)   = 0;
    pC{v}(s,s) = 0;
    
    P{1}.X     = [Y(:,r) U];
    P{1}.C     = C;
    P{2}.X     = pE{v}(u);
    P{2}.C     = pC{v}(u,u);
    PEB        = spm_PEB(dY,P,256);
    Ep{v}      = pE{v};
    Cp{v}      = pC{v};
    Ep{v}(u)   = PEB{2}.E;
    Cp{v}(u,u) = PEB{2}.C;
    
    if false % (not used)
        
        % Bayesian model reduction
        %------------------------------------------------------------------
        BMR.M.pE = pE{v};
        BMR.M.pC = pC{v};
        BMR.Ep   = Ep{v};
        BMR.Cp   = Cp{v};
        BMR      = spm_dcm_bmr_all(BMR,'all','BMS');
        
        % Jacobian for this voxel
        %------------------------------------------------------------------
        pE{v}    = BMR.M.pE;
        pC{v}    = BMR.M.pC;
        Ep{v}    = BMR.Ep;
        Cp{v}    = BMR.Cp;
        
    end
    
end

% Bayesian Model Reduction
%==========================================================================
spm_figure('GetWin','Bayesian Model Reduction');clf

% % dissipation (not used)
% %--------------------------------------------------------------------------
% r     = linspace(2*I,I/2,8);                      % self-inhibition (Hz)
% for k = 1:numel(r);
%     for i = 1:nv
%         
%         % self-connections
%         %------------------------------------------------------------------
%         rE      = pE{i};
%         rC      = pC{i};
%         rE(i)   = r(k);
%         rC(i,i) = rC(i,i)/16;
%         Fi      = spm_log_evidence(Ep{i},Cp{i},pE{i},pC{i},rE,rC);
%         
%         % log evidence (free energy)
%         %------------------------------------------------------------------
%         F(k,i) = Fi;
%         
%     end
% end
% 
% % pool evidence over regions
% %--------------------------------------------------------------------------
% F  = sum(F,2);
% 
% % Graphics
% %--------------------------------------------------------------------------
% 
% subplot(2,2,1)
% bar(r,F - min(F))
% title('Dissipation','fontsize',16)
% xlabel('self-inhibition (Hz)'), ylabel('log evidence')
% box off, axis square
% 
% 
% % reduced posterior
% %--------------------------------------------------------------------------
% [F,k] = max(F);
% for i = 1:nv
%     
%     % self-connections
%     %----------------------------------------------------------------------
%     rE      = pE{i};
%     rC      = pC{i};
%     rE(i)   = r(k);
%     rC(i,i) = rC(i,i)/4;
%     [F,sE,sC] = spm_log_evidence(Ep{i},Cp{i},pE{i},pC{i},rE,rC);
%     
%     % reduced posterior
%     %----------------------------------------------------------------------
%     Ep{i}  = sE;
%     Cp{i}  = sC;
%     pC{i}  = rC;
%     
% end


% long-range connectivity
%==========================================================================
r     = linspace(8,D,8);                        % connectivity radii
for k = 1:numel(r);
    for i = 1:nv
        
        % afferent connections
        %------------------------------------------------------------------
        j       = find(R(:,i) > r(k));
        rC      = pC{i};
        rC(j,j) = 0;
        Fi      = spm_log_evidence(Ep{i},Cp{i},pE{i},pC{i},pE{i},rC);
        
        % log evidence (free energy)
        %------------------------------------------------------------------
        F(k,i) = Fi;
        
    end
end

% pool evidence over regions
%--------------------------------------------------------------------------
F  = sum(F,2);

% Graphics
%--------------------------------------------------------------------------
subplot(1,3,1)
bar(r,F - min(F))
title('Connectivity radius','fontsize',16)
xlabel('radius (mm)'), ylabel('log evidence')
box off, axis square, spm_axis tight

% reduced posterior
%--------------------------------------------------------------------------
[F,k] = max(F);
r     = r(k);
for i = 1:nv
    
    % afferent connections
    %----------------------------------------------------------------------
    j         = find(R(:,i) > r);
    rC        = pC{i};
    rC(j,j)   = 0;
    [F,sE,sC] = spm_log_evidence(Ep{i},Cp{i},pE{i},pC{i},pE{i},rC);
    
    % reduced posterior
    %----------------------------------------------------------------------
    Ep{i}  = sE;
    Cp{i}  = sC;
    pC{i}  = rC;
    
end

% Symmetry priors on reciprocal connections
%==========================================================================

% reciprocal coupling shrinkage priors
%--------------------------------------------------------------------------
for i = 1:nv
    for j = (i + 1):nv
        
        if R(i,j) < r
            
            % afferent connection
            %--------------------------------------------------------------
            rCi          = pC{i};
            rCi(j,j)     = 0;
            [Fi,sEi,sCi] = spm_log_evidence(Ep{i},Cp{i},pE{i},pC{i},pE{i},rCi);
            
            % efferent connection
            %--------------------------------------------------------------
            rCj          = pC{j};
            rCj(i,i)     = 0;
            [Fj,sEj,sCj] = spm_log_evidence(Ep{j},Cp{j},pE{j},pC{j},pE{j},rCj);
            
            % accept new priors if F has increased
            %--------------------------------------------------------------
            F       = Fi + Fj;
            FR(i,j) = F;
            FR(j,i) = F;
            if F > 0
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

% reciprocal coupling - anti-symmetry priors (not used)
%--------------------------------------------------------------------------
% m     = 1/128;
% for i = 1:nv
%     for j = (i + 1):nv
%         
%         if Ep{i}(j)
%             
%             % reciprocal connections
%             %--------------------------------------------------------------
%             PC           = [pC{i}(j,j),0;0,pC{j}(i,i)];
%             CP           = [Cp{i}(j,j),0;0,Cp{j}(i,i)];
%             EP           = [Ep{i}(j);Ep{j}(i)];
%             PE           = [pE{i}(j);pE{j}(i)];
%             RC           = PC + [0,(m - 1);(m - 1),0];
%             
%             [Fi,sEi,sCi] = spm_log_evidence(EP,CP,PE,PC,PE,RC);
%             
%             % accept new priors if F has increased
%             %--------------------------------------------------------------
%             % if Fi > 0
%                 Ep{i}(j) = sEi(1);
%                 Ep{j}(i) = sEi(2);
%             % end
%         end
%     end
%     fprintf('evaluating (%i) of %i states\n',i,nv)
% end

% assemble Jacobian
%--------------------------------------------------------------------------
J     = zeros(nv + nu,nv);
for v = 1:nv
    J(:,v) = Ep{v};
end

% retain coupling between states
%--------------------------------------------------------------------------
K     = J((1:nu) + nv,:)';
J     = J((1:nv),(1:nv))';

% condition Jacobian (not used)
%==========================================================================
for i = 1:0
    
    % remove near zero coupling
    %----------------------------------------------------------------------
    J(abs(J) < 1/1024) = 0;
    
    % condition Jacobian
    %----------------------------------------------------------------------
    [e,v] = eig(J,'vector');
    v     = v - max(real(v)) - 1/64;
    J     = real(e*diag(v)*pinv(e));
    
    % remove near zero coupling
    %----------------------------------------------------------------------
    J(abs(J) < 1/1024) = 0;

end

% graphics
%==========================================================================

% reorder connections according to Euclidean proximity
%--------------------------------------------------------------------------
k     = 1;
j     = 1;
for v = 2:nv
    r     = R(:,j);
    r(j)  = Inf;
    r(k)  = Inf;
    [r,j] = min(r);
    k(v)  = j;
end

% Fraction of reciprocal connections eliminated
%--------------------------------------------------------------------------
fr  = sum(FR(:) > 0)/sum(FR(:) ~= 0);
str = sprintf('Redundant %2.0f %s',100*fr,'%');

% Sparse connectivity image
%--------------------------------------------------------------------------
subplot(1,3,2)
spy(FR(k,k) > 0,'.r',1/32), hold on
spy(FR(k,k) < 0,'.b',1/32), hold off
title('Sparse connectivity','fontsize',16)
xlabel('parcel'), ylabel('parcel')
box off, axis square

subplot(1,3,3), m = 32;
spy(FR(k(1:m),k(1:m)) > 0,'.r'), hold on
spy(FR(k(1:m),k(1:m)) < 0,'.b'), hold off
title(str,'fontsize',16)
xlabel('parcel'), ylabel('parcel')
box off, axis square
