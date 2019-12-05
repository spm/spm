function [PP] = spm_voice_get_LEX(xY,word,NI)
% Create lexical, prosody and speaker structures from word structures
% FORMAT [P] = spm_voice_get_LEX(xY,word)
%
% xY(nw,ns) -  structure array for ns samples of nw words
% word(nw)  -  cell array of word names
% NI(nw,ns) -  numeric array of number of minima
%
% updates or completes the global structure VOX:
%
% VOX.LEX(nw) -  structure array for nw words (lexical features)
% VOX.PRO(np) -  structure array for np features of prosody
% VOX.WHO(nq) -  structure array for nq features of speaker
%
% P           -  prosody parameters for exemplar (training) words
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
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_get_LEX.m 7750 2019-12-05 17:54:29Z spm $


% defaults
%--------------------------------------------------------------------------
global VOX
try, VOX.E; catch, VOX.E = 32; end    % regularisation (lexical)
try, VOX.G; catch, VOX.G = 8;  end    % covariance scaling (lexical)

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
Pstr  = {'amp','lat','dur','tim','Tu','Tv','Tf','Tw','p0','p1','p2'};
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
    LEX(w).dC   =   cov(P{w}');
    
end

% condition covariances - ensuring trace(qP*QC) = nq
%--------------------------------------------------------------------------
QC    = 0;
for w = 1:nw
    QC = QC + LEX(w).qC;
end
QC    = QC/nw;
c0    = trace(QC)*eye(size(QC))/nq/VOX.E;
for w = 1:nw
    qC        = LEX(w).qC;
    qC        = qC + c0;
    qC        = qC*trace(qC\QC)/nq;
    LEX(w).qC = VOX.G*qC;
end

% save mean and precision in voice structure
%--------------------------------------------------------------------------
VOX.W  = spm_zeros(xY(1).W);
VOX.qC = QC;


%% estimate prosody parameters
%==========================================================================
fS    = @(Q,G) spm_vec(spm_voice_iQ(spm_voice_Q(Q,G)));
G     = xY(1).P.pch;                                  % expansion point
j     = find(ismember(Pstr,{'Tu','Tv','Tf','Tw'}));   % pitch indices
jstr  = {'Const','W','Tu','Tv','Tf','Tw'};
for w = 1:nw
    
    % setup Taylor approximation to variations prosody
    %----------------------------------------------------------------------
    W         = LEX(w).qE/std(LEX(w).qE(:));        % expansion point
    dWdP      = spm_diff(fS,reshape(W,Nu,Nv),G,2);  % partial derivatives
    LEX(w).qE = W;                                  
    LEX(w).X  = [ones(numel(W),1) W(:) dWdP];       % design matrix (GLM)
    LEX(w).pE = [0 1 0 0 0 0]';                     % pitch expectation
    LEX(w).pC = diag([1 1 1 1 1 1]/8);              % pitch covariance
    
end

for i = 1:4
    for w = 1:nw
        
        % general linear model
        %------------------------------------------------------------------
        X    = LEX(w).X;                            % GLM
        qC   = LEX(w).qC;                           % covariance (lexical)
        pE   = LEX(w).pE;                           % expectation (pitch)
        pC   = LEX(w).pC;                           % covariance  (pitch)
        
        iqC  = spm_inv(qC);                         % precision (W)
        ipC  = spm_inv(pC);                         % precision (P)
        Cp   = (X'*iqC*X + ipC);                    % conditional precision
        
        % pitch parameters
        %------------------------------------------------------------------
        for s = 1:ns
            
            % esimate variations in prosody from expansion point
            %--------------------------------------------------------------
            dP(:,s)  = Cp\(X'*iqC*xY(w,s).W(:) + ipC*pE);
            
            % augment pitch parameters
            %--------------------------------------------------------------
            P{w}(j,s) = P{w}(j,s) + dP(3:end,s);
            
        end
        
        % expected formant coeficients and moments of pitch
        %------------------------------------------------------------------
        LEX(w).pE  = LEX(w).pE/2 + mean(dP,2)/2;
        LEX(w).pC  = LEX(w).pC/2 + cov(dP')/2;
               
        % graphics
        %------------------------------------------------------------------
        if w < 5 && i < 2
            spm_figure('GetWin','Basis functions');
            nx    = size(X,2);
            for ii = 1:nx
                subplot(4,nx,nx*(w - 1) + ii)
                imagesc(spm_voice_Q(reshape(X(:,ii),Nu,Nv),G))
                title(jstr{ii})
            end
        end
    end
end

% precompute projectors and entropy terms
%--------------------------------------------------------------------------
for w = 1:nw
    
    X    = LEX(w).X;                            % GLM
    qC   = LEX(w).qC;                           % covariance (lexical)
    pE   = LEX(w).pE;                           % expectation (pitch)
    pC   = LEX(w).pC;                           % covariance  (pitch)
    
    iqC  = spm_inv(qC);                         % precision (W)
    ipC  = spm_inv(pC);                         % precision (P)
    Cp   = (X'*iqC*X + ipC);                    % conditional precision
    
    LEX(w).pP = Cp\(ipC*pE);
    LEX(w).M  = Cp\(X'*iqC);
    LEX(w).L  = - spm_logdet(qC)/2 ...           % precision
                - spm_logdet(pC*Cp)/2;           % entropy
    
end

%% prosody {'amp','lat','dur','tim','Tu','Tv','Tf','Tw','p0','p1','p2'}
%==========================================================================
PP     = full(spm_cat(P)');
RR     = full(spm_cat(R)');

% prosody ranges
%--------------------------------------------------------------------------
D      = mean(PP);                               % mean
sd     = std(PP);                                % and standard deviation
D      = [D - 3*sd; D + 3*sd]';                  % ranges of prosody
D(1,:) = log([1/64  1]);                         % amplitude
D(2,:) = log([1/32  1]);                         % latency (sec)

% select prosody features and specify prior precision
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
D(1,:) = log([85 300]);                           % ff0
D(2,:) = log([30  45]);                           % ff1

% select prosody features and specify prior precision
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

% number of minima for each word
%--------------------------------------------------------------------------
FI    = zeros(3,nw);
for w = 1:nw
    for s = 1:size(NI,2)
        i       = NI(w,s);
        FI(i,w) = FI(i,w) + 1;
    end
end
FI     = bsxfun(@rdivide,FI,sum(FI));
VOX.FI = FI;

% Create a default prosody states for each word
%--------------------------------------------------------------------------
p0     = spm_vec(VOX.P);
pros   = zeros(np,nw);
for w  = 1:nw
    for p = 1:np
        d         = abs(LEX(w).dE(p) - (PRO(p).pE + p0(p)));
        [d,i]     = min(d);
        pros(p,w) = i;
    end
end
VOX.prosody = pros;



% illustrate distribution of lexical and prosody parameters
%==========================================================================
spm_figure('GetWin','Parameter distributions'); clf

% indices for plotting
%--------------------------------------------------------------------------
for i = 1:numel(Pstr)
    subplot(4,3,i)
    if i > 7
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

for i = 1:numel(Rstr)
    subplot(4,2,i)
    hist(exp(RR(:,i)),32), axis square
    str = sprintf('%s mean: %.2f (%.2f)',Rstr{i},mean(exp(RR(:,i))),std(exp(RR(:,i))));
    title(str)
end

i    = full(spm_cat(I));
subplot(4,2,3), hist(i(1,:),32), axis square
title(sprintf('%s mean (sd): %.2f (%.3f)','onset',mean(i(1,:)),std(i(1,:))))
subplot(4,2,4), hist(i(2,:),32), axis square
title(sprintf('%s mean (sd): %.2f (%.3f)','offset',mean(i(2,:)),std(i(2,:))))


%% Eigenmodes of lexical (U) and prosody (V) parameters
%==========================================================================
spm_figure('GetWin','Eigen-reduction'); clf

% lexical
%--------------------------------------------------------------------------
[U,S] = spm_svd(cov(spm_cat(Q)'));
S     = log(diag(S));
s     = find(S > (max(S) - 4),1,'last');
subplot(2,2,1), bar(S)
title('Eigenvalues - lexical','FontSize',16), ylabel('log value')
hold on, plot([s,s],[0,max(S)],'r:'), hold off
xlabel(sprintf('eigenbasis (%i)',s)), axis square

% prosody (based on correlation)
%--------------------------------------------------------------------------
[V,S] = spm_svd(corrcoef(PP));
S     = diag(S);
subplot(2,2,2), bar(S)
title('Eigenvalues - prosody','FontSize',16)
xlabel('eigenbasis'), ylabel('eigenvalue'), axis square

subplot(2,1,2), bar(abs(V)')
title('Eigenmodes - prosody','FontSize',16)
xlabel('prosody mode'), ylabel('amplitude'), axis square
legend(Pstr)
