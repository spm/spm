function [LEX,PRO,WHO] = spm_voice_get_LEX(xY,word)
% Creates lexical, prosody and identity structures from word structures
% FORMAT [LEX,PRO,WHO] = spm_voice_get_LEX(xY,word)
%
% xY(nw,ns) -  structure array for ns samples of nw words
% word(nw)  -  cell array of word names
%
% LEX(nw,nk) -  structure array for nk variants of nw words
% PRO(np)    -  structure array for np aspects of prosody
% WHO(nq)    -  structure array for nq aspects of idenity
%
%  This routine creates a triplet of structure arrays used to infer the
%  lexical content and prosody of a word -  and the identity of the person
%  talking. It uses  exemplar word files, each containing 32 words spoken
%  with varying prosody.  each structure contains the expectations and
%  precisions of lexical and prosody parameters (Q and P respectively) -
%  and associated eigenbases. This allows the likelihood of any given word
%  (summarised in a word structure xY)  to be evaluated  under Gaussian
%  assumptions about random fluctuations in parametric space. The identity
%  and prosody likelihoods are based upon the prosody parameters, while the
%  lexical likelihood is based upon the lexical parameters.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_get_LEX.m 7535 2019-03-03 20:27:09Z karl $


% defaults
%--------------------------------------------------------------------------
global voice_options
try, E = voice_options.E; catch, E  = 64; end     % regularisation


%% assemble parameters for subsequent analysis
%==========================================================================
[nw,ns] = size(xY);
for w   = 1:nw
    for s = 1:ns
        Q{w}(:,s) = spm_vec(xY(w,s).Q);
        P{w}(:,s) = spm_vec(xY(w,s).P);
    end
end

% concatenate and illustrate distribution of prosody parameters
%==========================================================================
spm_figure('GetWin','Parameter distributions'); clf
PP    = full(spm_cat(P)');
Pstr  = {'ff0','ff1','dur','timbre','const','inf','bif'};
for i = 1:size(PP,2)
    subplot(3,3,i)
    hist(exp(PP(:,i)),32), axis square
    title(sprintf('%s mean: %.2f',Pstr{i},mean(exp(PP(:,i)))))
end


%% Eigenmodes of lexical (U) and prosody (V) parameters
%==========================================================================
spm_figure('GetWin','Eigen reduction'); clf

% lexical
%--------------------------------------------------------------------------
[U,S] = spm_svd(cov(spm_cat(Q)'));
S     = log(diag(S));
s     = find([S;1] > (max(S) - 4),1,'last');
subplot(2,2,1), bar(S)
title('Eigenvalues - lexical','FontSize',16), ylabel('log value')
hold on, plot([s,s],[0,max(S)],'r:'), hold off
xlabel(sprintf('eigenbasis (%i)',s)), axis square

% prosody (based on correlation)
%--------------------------------------------------------------------------
[V,S] = spm_svd(cov(PP));
S     = diag(S);
subplot(2,2,2), bar(S)
title('Eigenvalues - prosidy','FontSize',16)
xlabel('eigenbasis'), ylabel('eigenvalue'), axis square

subplot(2,1,2), bar(abs(V)')
title('Eigenmodes - prosidy','FontSize',16)
xlabel('prosody mode'), ylabel('amplitude'), axis square
legend(Pstr)


%% prosody {'ff0','ff1','dur','timbre','const','inf','bif'};
%==========================================================================

% prosidy ranges
%--------------------------------------------------------------------------
R(1,:)   = log([96  350]);                           % ff0
R(2,:)   = log([24   64]);                           % ff1
R(3,:)   = log([1/4 3/4]);                           % dur
R(4,:)   = log([1.8 2.8]);                           % timbre
R(5,:)   = [1 1];                                    % const
R(6,:)   = [-1/4 1/4];                               % inf
R(7,:)   = [-1/4 1/4];                               % bif

% select prosidy features and specify prior precision
%--------------------------------------------------------------------------
i        = [3,4,6,7];
ni       = numel(i);
p0       = mean(R,2);
PRO(1).P = spm_unvec(p0(:),xY(1).P);
PRO(1).i = i;

% mixture of Gaussians
%--------------------------------------------------------------------------
k     = 8;      
for p = 1:ni
    
    % prior densities
    %----------------------------------------------------------------------
    pE  = linspace(R(i(p),1),R(i(p),2),k) - p0(i(p));
    pC  = diff(R(i(p),:))/(k - 1)/4;
    pC  = pC^2;
   
    % save prior densities
    %----------------------------------------------------------------------
    PRO(p).str = Pstr{i(p)};
    PRO(p).pE  = pE(:);
    PRO(p).pC  = ones(k,1)*pC;
    PRO(p).pP  = ones(k,1)/pC;
   
end

%% identity
%==========================================================================

% select identity features and specify prior precision
%--------------------------------------------------------------------------
i        = [1 2];
ni       = numel(i);
WHO(1).i = i;

% mixture of Gaussians
%--------------------------------------------------------------------------
k     = 8;
for p = 1:ni
    
    % prior densities
    %----------------------------------------------------------------------
    pE  = linspace(R(i(p),1),R(i(p),2),k) - p0(i(p));
    pC  = diff(R(i(p),:))/(k - 1)/4;
    pC  = pC^2;
   
    % save prior densities
    %----------------------------------------------------------------------
    WHO(p).str = Pstr{i(p)};
    WHO(p).pE  = pE(:);
    WHO(p).pC  = ones(k,1)*pC;
    WHO(p).pP  = ones(k,1)/pC;
   
end


%% lexical
%==========================================================================

% evaluate word specific means and covariances in eigenmode space
%--------------------------------------------------------------------------
N     = 1;                                     % number of variants
q0    = mean(spm_cat(Q),2);                    % average word
Q0    = spm_unvec(q0,xY(1).Q);
for w = 1:nw
    
    % lexical coefficients
    %----------------------------------------------------------------------
    LEX(w,1).word   = word{w};
    Qw              = bsxfun(@minus,Q{w},q0)';
    [pL,pE,pC]      = spm_kmeans(Qw,N,'fixed-points');
    [j,k]           = sort(pL,'descend');
  
    % word variants: mean and variance
    %----------------------------------------------------------------------
    for i = k
        LEX(w,i).qE = pE(i,:)';
        LEX(w,i).qC = pC(:,:,i);
        LEX(w,i).Q  = Q0;
    end
end

% covariance normalisation; based on within word covariance
%--------------------------------------------------------------------------
QC    = 0;
for w = 1:nw
    for i = N
        QC = QC + LEX(w,i).qC;
    end
end
QC    = QC/(nw*N);
q0    = trace(QC)*eye(size(QC))/E;
for w = 1:nw
    for k = 1:N
        qC          = LEX(w,k).qC;
        qC          = trace(QC)*qC/trace(qC);
        qC          = qC + q0;
        qP          = spm_inv(qC);
        
        LEX(w,k).qC = qC;
        LEX(w,k).qP = qP;
    end
end

% find most likely exemplar
%--------------------------------------------------------------------------
for w = 1:nw
    for k = 1:N
        
        % likelihood
        %------------------------------------------------------------------
        L           = spm_voice_likelihood(xY(w,:),LEX(w,:),PRO,WHO);
        [L,i]       = max(spm_vec(L));
        LEX(w,k).rE = spm_vec(xY(w,i).Q) - spm_vec(LEX(1).Q);
        
    end
end

return
