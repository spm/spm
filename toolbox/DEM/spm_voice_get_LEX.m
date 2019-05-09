function [PP] = spm_voice_get_LEX(xY,word)
% Creates lexical, prosody and speaker structures from word structures
% FORMAT [P] = spm_voice_get_LEX(xY,word)
%
% xY(nw,ns)      -  structure array for ns samples of nw words
% word(nw)       -  cell array of word names
%
% updates or completes the global structure VOX:
%
% VOX.LEX(nw,nk) -  structure array for nk variants of nw words
% VOX.PRO(np)    -  structure array for np aspects of prosody
% VOX.WHO(nq)    -  structure array for nq aspects of speaker
%
% P              -  lexical parameters for exemplar (training) words
%
%  This routine creates a triplet of structure arrays used to infer the
%  lexical content and prosody of a word - and the identity of the person
%  talking (in terms of the vocal tract, which determines F1). It uses
%  exemplar word files, each containing 32 words spoken with varying
%  prosody. Each structure contains the expectations and precisions of
%  lexical and prosody parameters (Q and P respectively) - and associated
%  eigenbases. This allows the likelihood of any given word (summarised in
%  a word structure xY)  to be evaluated  under Gaussian assumptions about
%  random fluctuations in parametric space. The identity and prosody
%  likelihoods are based upon the prosody parameters, while the lexical
%  likelihood is based upon the lexical parameters. These (LEX, PRO, and
%  WHO)structures are placed in the VOX structure, which is a global
%  variable. In addition, the expected value of various coefficients are
%  stored in VOX.Q and VOX.P.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_get_LEX.m 7589 2019-05-09 12:57:23Z karl $


% defaults
%--------------------------------------------------------------------------
global VOX
try, E       = VOX.E; catch, E  = 4; end                % regularisation
VOX.analysis = 0;
VOX.graphics = 0;
VOX.interval = 0;
VOX.mute     = 0;
VOX.onsets   = 0;


%% assemble parameters for subsequent analysis
%==========================================================================
[nw,ns] = size(xY);
for w   = 1:nw
    for s = 1:ns
        Q{w}(:,s) = spm_vec(xY(w,s).Q);
        P{w}(:,s) = spm_vec(xY(w,s).P);
        R{w}(:,s) = spm_vec(xY(w,s).R);
        I{w}(:,s) = spm_vec(xY(w,s).i);
    end
end

% concatenate and illustrate distribution of prosody parameters
%==========================================================================
spm_figure('GetWin','Parameter distributions'); clf

%       {'amp','lat','dur','tim','p0','p1','p2','p3'}
%--------------------------------------------------------------------------
Pstr  = {'amp','lat','dur','tim','p0','p1','p2','p3'};
Rstr  = {'F0','F1',};
PP    = full(spm_cat(P)');

% indices for plotting
%--------------------------------------------------------------------------
for i = 1:numel(Pstr);
    subplot(3,3,i)
    
    if i > 4
        hist(PP(:,i),32), axis square
        title(sprintf('%s mean: %.2f',Pstr{i},mean(PP(:,i))))
    else
        hist(exp(PP(:,i)),32), axis square
        title(sprintf('%s mean: %.2f',Pstr{i},mean(exp(PP(:,i)))))
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


%% prosody {'amp','lat','dur','tim','p0','p1','p2','p3'}
%==========================================================================

% prosidy ranges
%--------------------------------------------------------------------------
D      = mean(PP);
sd     = std(PP);
D      = [D - 4*sd; D + 4*sd]';                  % all prosody parameters
D(1,:) = log([1/32 1]);                          % amp
D(2,:) = log([1/32 1]);                          % lat

% select prosidy features and specify prior precision
%--------------------------------------------------------------------------
p0    = mean(D,2);
VOX.P = spm_unvec(p0(:),xY(1).P);

% mixture of Gaussians
%--------------------------------------------------------------------------
k     = 8;      
for p = 1:size(P{1},1)
    
    % prior densities
    %----------------------------------------------------------------------
    pE  = linspace(D(p,1),D(p,2),k) - p0(p);
    pC  = diff(D(p,:))/(k - 1)/4;
    pC  = pC^2;
   
    % save prior densities
    %----------------------------------------------------------------------
    PRO(p).str = Pstr{p};
    PRO(p).pE  = pE(:);
    PRO(p).pC  = ones(k,1)*pC;
   
end

%% speaker
%==========================================================================
clear D
D(1,:) = log([80 300]);                          % ff0
D(2,:) = log([25  50]);                          % ff1

% select prosidy features and specify prior precision
%--------------------------------------------------------------------------
p0    = mean(D,2);
VOX.R = spm_unvec(p0(:),xY(1).R);

% mixture of Gaussians
%--------------------------------------------------------------------------
k     = 16;
for p = 1:size(R{1},1)
    
    % prior densities
    %----------------------------------------------------------------------
    pE  = linspace(D(p,1),D(p,2),k) - p0(p);
    pC  = diff(D(p,:))/(k - 1)/4;
    pC  = pC^2;
   
    % save prior densities
    %----------------------------------------------------------------------
    WHO(p).str = Rstr{p};
    WHO(p).pE  = pE(:);
    WHO(p).pC  = ones(k,1)*pC;

end


%% lexical
%==========================================================================

% evaluate word specific means and covariances
%--------------------------------------------------------------------------
N     = 1;                                     % number of variants
np    = size(Q{1},1);                          % number of parameters
q0    = mean(spm_cat(Q),2);                    % average word
for w = 1:nw
    
    % lexical coefficients
    %----------------------------------------------------------------------
    LEX(w,1).word   = word{w};
    Qw              = bsxfun(@minus,Q{w},q0)';
    [pL,pE,pC]      = spm_kmeans(Qw,N,'fixed-points');
    [j,k]           = sort(pL,'descend');
    
    
    % mean and precision of duration
    %----------------------------------------------------------------------
    dE    = mean(P{w},2);
    dC    = cov(P{w}')/8;
    
    % word variants: mean and variance
    %----------------------------------------------------------------------
    for i = k
        LEX(w,i).qE = pE(i,:)';
        LEX(w,i).qC = pC(:,:,i);
        
        % duration priors
        %----------------------------------------------------------------------
        LEX(w,i).dE = dE;
        LEX(w,i).dC = dC;
    end
end


% precision normalisation; ensuring trace(qP*QC) = np
%--------------------------------------------------------------------------
QC    = 0;
for w = 1:nw
    for k = N
        QC = QC + LEX(w,k).qC;
    end
end
QC    = QC/(nw*N);
c0    = trace(QC)*eye(size(QC))/np/E;
for w = 1:nw
    for k = 1:N
        qC          = LEX(w,k).qC;
        qC          = qC + c0;
        qC          = qC*trace(qC\QC)/np;
        LEX(w,k).qC = qC;
    end
end


% save mean and precision in voice structure
%--------------------------------------------------------------------------
VOX.Q   = spm_unvec(q0,xY(1).Q);
VOX.qC  = QC + c0;

% place lexical and other structures voice structure
%--------------------------------------------------------------------------
VOX.LEX = LEX;
VOX.PRO = PRO;
VOX.WHO = WHO;


return
