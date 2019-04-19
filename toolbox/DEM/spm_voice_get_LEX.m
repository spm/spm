function [PP] = spm_voice_get_LEX(xY,word)
% Creates lexical, prosody and speaker structures from word structures
% FORMAT [P] = spm_voice_get_LEX(xY,word)
%
% xY(nw,ns)      -  structure array for ns samples of nw words
% word(nw)       -  cell array of word names
%
% VOX.LEX(nw,nk) -  structure array for nk variants of nw words
% VOX.PRO(np)    -  structure array for np aspects of prosody
% VOX.WHO(nq)    -  structure array for nq aspects of speaker
%
% P              -  lexical parameters for exemplar (training) words
%
%  This routine creates a triplet of structure arrays used to infer the
%  lexical content and prosody of a word -  and the identity of the person
%  talking. It uses  exemplar word files, each containing 32 words spoken
%  with varying prosody. Each structure contains the expectations and
%  precisions of lexical and prosody parameters (Q and P respectively) -
%  and associated eigenbases. This allows the likelihood of any given word
%  (summarised in a word structure xY)  to be evaluated  under Gaussian
%  assumptions about random fluctuations in parametric space. The identity
%  and prosody likelihoods are based upon the prosody parameters, while the
%  lexical likelihood is based upon the lexical parameters. These (LEX,
%  PRO, and WHO)structures are placed in the VOX structure, which is a
%  global variable. In addition, the expected value of various coefficients
%  are stored in VOX.Q and VOX.P.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_get_LEX.m 7574 2019-04-19 20:38:15Z karl $


% defaults
%--------------------------------------------------------------------------
global VOX
try, E = VOX.E; catch, E  = 64; end     % regularisation


%% assemble parameters for subsequent analysis
%==========================================================================
[nw,ns] = size(xY);
for w   = 1:nw
    for s = 1:ns
        Q{w}(:,s) = spm_vec(xY(w,s).Q);
        P{w}(:,s) = spm_vec(xY(w,s).P);
        I{w}(:,s) = spm_vec(xY(w,s).i);
    end
end

% concatenate and illustrate distribution of prosody parameters
%==========================================================================
spm_figure('GetWin','Parameter distributions'); clf
PP    = full(spm_cat(P)');
Pstr  = {'amp','lat','ff1','dur','timbre','ff0','inf','bif'};
ii    = [1 4 5 6 7 8];
for i = 1:numel(ii);
    subplot(3,3,i)
    
    if ii(i) > 5
        hist(PP(:,ii(i)),32), axis square
        title(sprintf('%s mean: %.2f',Pstr{ii(i)},mean(PP(:,ii(i)))))
    else
        hist(exp(PP(:,ii(i))),32), axis square
        title(sprintf('%s mean: %.2f',Pstr{ii(i)},mean(exp(PP(:,ii(i))))))
    end
end

% onsets and offsets
%--------------------------------------------------------------------------
i    = full(spm_cat(I));
subplot(3,2,5), hist(i(1,:),32,'Color','c'), axis square
title(sprintf('%s mean (sd): %.2f (%.3f)','onset',mean(i(1,:)),std(i(1,:))))
subplot(3,2,6), hist(i(2,:),32,'c'), axis square
title(sprintf('%s mean (sd): %.2f (%.3f)','offset',mean(i(2,:)),std(i(2,:))))


%% Eigenmodes of lexical (U) and prosody (V) parameters
%==========================================================================
spm_figure('GetWin','Eigen-reduction'); clf

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
[V,S] = spm_svd(corr(PP));
S     = diag(S);
subplot(2,2,2), bar(S)
title('Eigenvalues - prosidy','FontSize',16)
xlabel('eigenbasis'), ylabel('eigenvalue'), axis square

subplot(2,1,2), bar(abs(V)')
title('Eigenmodes - prosidy','FontSize',16)
xlabel('prosody mode'), ylabel('amplitude'), axis square
legend(Pstr)


%% prosody {'amp','lat','ff1','dur','timbre','ff0','inf','bif'}
%==========================================================================

% prosidy ranges
%--------------------------------------------------------------------------
R      = [min(PP); max(PP)]';
R(1,:) = log([1/32   1]);                          % amp
R(2,:) = log([1/32   1]);                          % lat
R(3,:) = log([24    48]);                          % ff1


% select prosidy features and specify prior precision
%--------------------------------------------------------------------------
i     = [1,2,4,5,6,7,8];
ni    = numel(i);
p0    = mean(R,2);
VOX.P = spm_unvec(p0(:),xY(1).P);
VOX.i = i;

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

%% speaker
%==========================================================================

% select speaker features and specify prior precision
%--------------------------------------------------------------------------
i     = [3];
ni    = numel(i);
VOX.j = i;

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
c0    = trace(QC)*eye(size(QC))/E;
for w = 1:nw
    for k = 1:N
        qC          = LEX(w,k).qC;
        qC          = trace(QC)*qC/trace(qC);
        qC          = qC + c0;
        qP          = spm_inv(qC);
        
        LEX(w,k).qC = qC;
        LEX(w,k).qP = qP;
    end
end


% save mean and precision in voice structure
%--------------------------------------------------------------------------
VOX.Q  = spm_unvec(q0,xY(1).Q);
VOX.qP = spm_inv(QC + c0);

% place lexical and other structures voice structure
%--------------------------------------------------------------------------
VOX.LEX = LEX;
VOX.PRO = PRO;
VOX.WHO = WHO;

% find most likely exemplar for word generation
%--------------------------------------------------------------------------
for w = 1:nw
    for k = 1:N
        
        % likelihood
        %------------------------------------------------------------------
        L               = spm_voice_likelihood(xY(w,:),w);
        [L,i]           = max(L(w,1,1,:));
        VOX.LEX(w,k).rE = spm_vec(xY(w,i).Q) - spm_vec(VOX.Q);
        
    end
end

return
