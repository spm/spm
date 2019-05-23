function [PP] = spm_voice_get_LEX(xY,word)
% Creates lexical, prosody and speaker structures from word structures
% FORMAT [P] = spm_voice_get_LEX(xY,word)
%
% xY(nw,ns)      -  structure array for ns samples of nw words
% word(nw)       -  cell array of word names
%
% updates or completes the global structure VOX:
%
% VOX.LEX(nw)    -  structure array for nw words (lexical features)
% VOX.PRO(np)    -  structure array for np features of prosody
% VOX.WHO(nq)    -  structure array for nq features of speaker
%
% P              -  prosidy parameters for exemplar (training) words
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
%  stored in VOX.W and VOX.P.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_get_LEX.m 7597 2019-05-23 18:42:38Z karl $


% defaults
%--------------------------------------------------------------------------
global VOX
try, E       = VOX.E; catch, E = 2; end    % regularisation
VOX.analysis = 0;
VOX.graphics = 0;
VOX.interval = 0;
VOX.mute     = 0;
VOX.onsets   = 0;


%% assemble parameters for subsequent analysis
%==========================================================================
[nw,ns] = size(xY);                        % number of words and exemplars
[Nu,Nv] = size(xY(1).W);                   % fifth order of basis functions
for w   = 1:nw
    for s = 1:ns
        Q{w}(:,s) = spm_vec(xY(w,s).W);
        P{w}(:,s) = spm_vec(xY(w,s).P);
        R{w}(:,s) = spm_vec(xY(w,s).R);
        I{w}(:,s) = spm_vec(xY(w,s).i);
    end
end

% label strings
%--------------------------------------------------------------------------
Pstr  = {'amp','lat','dur','tim','Tu','Tv','Tw','p0','p1','p2'};
Rstr  = {'F0','F1',};


%% lexical
%==========================================================================

% evaluate word specific means and covariances
%--------------------------------------------------------------------------
nq    = size(Q{1},1);                      % number of lexical parameters
np    = size(P{1},1);                      % number of prosody parameters
nr    = size(R{1},1);                      % number of speaker parameters
for w = 1:nw
    
    % lexical label
    %----------------------------------------------------------------------
    LEX(w).word = word{w};
    
    % moments of lexical parameters
    %----------------------------------------------------------------------
    LEX(w).qE   = mean(Q{w},2);
    LEX(w).qC   =  cov(Q{w}');
    
    % moments of prosody parameters
    %----------------------------------------------------------------------
    LEX(w).dE   = mean(P{w},2);
    LEX(w).dC   =  cov(P{w}')/8;
    
end

% condition covariances – ensuring trace(qP*QC) = nq
%--------------------------------------------------------------------------
QC    = 0;
for w = 1:nw
    QC = QC + LEX(w).qC;
end
QC    = QC/nw;
c0    = E*trace(QC)*eye(size(QC))/nq;
for w = 1:nw
    qC        = LEX(w).qC;
    qC        = qC + c0;
    qC        = qC*trace(qC\QC)/nq;
    LEX(w).qC = qC;
end

% save mean and precision in voice structure
%--------------------------------------------------------------------------
VOX.W   = spm_zeros(xY(1).W);
VOX.qC  = QC;


%% estimate timbre parameters
%==========================================================================
fS    = @(Q,G) spm_vec(spm_voice_iQ(spm_voice_Q(Q,G)));
G     = xY(1).P.pch;                                    % expansion point
j     = find(ismember(Pstr,{'Tu','Tv','Tw'}));          % pitch indices
for w = 1:nw
    
    % setup Taylor approximation to variations in timbre
    %----------------------------------------------------------------------
    X     = spm_diff(fS,reshape(LEX(w).qE,Nu,Nv),G,2);

    % save MAP projector in LEX
    %----------------------------------------------------------------------
    LEX(w).dWdP  = X;
    
    % set up PEB
    %----------------------------------------------------------------------
    for s = 1:ns
        
        % esimate variations in pitch from expansion point
        %------------------------------------------------------------------
        qC   = LEX(w).qC;                    % covariance
        dW   = xY(w,s).W(:) - LEX(w).qE;      % residuals
        dP   = (X'*(qC\X))\X'*(qC\dW);
        
        % add to timbre parameters
        %------------------------------------------------------------------
        P{w}(j,s) = P{w}(j,s) + dP;
        
    end
    
    % moments of pitch
    %----------------------------------------------------------------------
    LEX(w).pE  = mean(P{w}(j,:),2) - G(:);
    LEX(w).pC  =  cov(P{w}(j,:)');
    
    
    % graphics
    %----------------------------------------------------------------------
    if w < 5
        spm_figure('GetWin','Basis functions');
        qX    = [LEX(w).qE(:) X];
        nx    = size(qX,2);
        for i = 1:nx
            subplot(4,nx,nx*(w - 1) + i)
            imagesc(spm_voice_Q(reshape(qX(:,i),Nu,Nv),G))
        end
    end
end


%% prosody {'amp','lat','dur','tim','Tu','Tv','Tw','p0','p1','p2'}
%==========================================================================
PP     = full(spm_cat(P)');
RR     = full(spm_cat(R)');

% prosidy ranges
%--------------------------------------------------------------------------
D      = mean(PP);                               % mean
sd     = std(PP);                                % and standard deviation
D      = [D - 3*sd; D + 3*sd]';                  % ranges of prosody
D(1,:) = log([1/32  1]);                         % amplitude
D(2,:) = log([1/32  1]);                         % latency
D(5,:) = log([3.5 4.5]);                         % pitch (Tu)

% select prosidy features and specify prior precision
%--------------------------------------------------------------------------
p0    = mean(D,2);
VOX.P = spm_unvec(p0(:),xY(1).P);

% mixture of Gaussians
%--------------------------------------------------------------------------
k     = 8;
for p = 1:np
    
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


%% identity
%==========================================================================
clear D
D(1,:) = log([80 300]);                          % ff0
D(2,:) = log([30  45]);                          % ff1

% select prosidy features and specify prior precision
%--------------------------------------------------------------------------
p0    = mean(D,2);
VOX.R = spm_unvec(p0(:),xY(1).R);

% mixture of Gaussians
%--------------------------------------------------------------------------
k     = 16;
for p = 1:nr
    
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


% place lexical and other structures voice structure
%--------------------------------------------------------------------------
VOX.LEX = LEX;
VOX.PRO = PRO;
VOX.WHO = WHO;


% illustrate distribution of lexical and prosody parameters
%==========================================================================
spm_figure('GetWin','Parameter distributions'); clf

% indices for plotting
%--------------------------------------------------------------------------
for i = 1:numel(Pstr);
    subplot(4,3,i)
    if i > 6
        hist(PP(:,i),32), axis square
        str = sprintf('%s mean: %.2f (%.2f)',Pstr{i},mean(PP(:,i)),std(PP(:,i)));
        title(str)
    else
        hist(exp(PP(:,i)),32), axis square
        str = sprintf('%s mean: %.2f (%.2f)',Pstr{i},mean(exp(PP(:,i))),std(exp(PP(:,i))));
        title(str)
    end
end

%% frequencies, onsets and offsets
%--------------------------------------------------------------------------
spm_figure('GetWin','Durations and frequencies');

for i = 1:numel(Rstr);
    subplot(4,2,i)
    hist(exp(RR(:,i)),32), axis square
    str = sprintf('%s mean: %.2f (%.2f)',Rstr{i},mean(exp(RR(:,i))),std(exp(RR(:,i))));
    title(str)
end

i    = full(spm_cat(I));
subplot(4,2,3), hist(i(1,:),32,'Color','c'), axis square
title(sprintf('%s mean (sd): %.2f (%.3f)','onset',mean(i(1,:)),std(i(1,:))))
subplot(4,2,4), hist(i(2,:),32,'c'), axis square
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



return
