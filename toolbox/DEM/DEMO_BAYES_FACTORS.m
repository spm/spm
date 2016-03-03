function DEMO_BAYES_FACTORS(pC,hE,hC,N,b)
% FORMAT DEMO_BAYES_FACTORS(pC,hE,hC,N)
% Demonstration Bayes factors and classical p-values
%--------------------------------------------------------------------------
% pC   - prior covariance             (e.g., 4)
% hE   - expectation of log precision (e.g., 1)
% hC   - covariance of log precision  (e.g., 1/8)
% N    - number of observations       (e.g., 16)
% b    - relative variance under alternate and null (e.g., 1/32)
%
% This demonstration routine uses a simple linear model to examine the
% relationship between free energy differences or log Bayes factors and
% classical F statistics. Using re-randomisation of a design matrix, it
% computes the null distribution over both statistics and plots them
% against each other.  There is a linear relationship, which allows one to
% evaluate the false-positive rate for any threshold on the Bayes factor.
% Ideally, one would like to see a positive log Bayes factor map to a
% classical threshold of p=0.05. The offset and slope of the linear
% relationship between the two statistics depends upon prior beliefs about
% the covariance of the parameters and the log precision. These can be
% changed by editing the code below (or supplying input arguments).
%__________________________________________________________________________
% Copyright (C) 2010-2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Peter Zeidman
% $Id: DEMO_BAYES_FACTORS.m 6737 2016-03-03 12:05:51Z karl $


% set up
%--------------------------------------------------------------------------
rng('default')

try, pC; catch, pC = 1;     end
try, hE; catch, hE = 0;     end
try, hC; catch, hC = 1/16;  end
try, N;  catch, N  = 32;    end
try, b;  catch, b  = 1/128; end

Y    = randn(N,1);
XX   = kron([1 1; 1 -1],ones(N/2,1));

% Model specification
%==========================================================================
M.nograph = 1;
M.noprint = 1;

M.IS = @(P,M,U) U*P;
M.pE = zeros(2,1);
M.pC = speye(2,2)*pC;
M.hE = hE;
M.hC = hC;
 

% re-randomisation
%--------------------------------------------------------------------------
Ns  = 64;
for i = 1:Ns
    X{i} = XX;
    X{i}(:,2) = XX(randperm(N),2);
end

pE    = M.pE;                    % full prior expectations
pC    = M.pC;                    % full prior covariance
rC    = pC; rC(2,2) = b;         % restricted or reduced priors
R     = M; R.pC = rC;            % reduced model  
Cr    = [0 0;0 1];               % classical contrast
for i = 1:Ns
    
    % Bayesian analysis (full comparison and model reduction)
    %----------------------------------------------------------------------
    [qE,qC,Eh,f] =  spm_nlsi_GN(M,X{i},Y);
    F(i,1)       = -spm_log_evidence(qE,qC,pE,pC,pE,rC);
    [qE,qC,Eh,r] =  spm_nlsi_GN(R,X{i},Y);
    G(i,1)       =  f - r; disp(i)
    
    % classical analysis
    %----------------------------------------------------------------------
    T(i,1)       =  spm_ancova(X{i},[],Y,Cr);
    
end

% classical threshold
%--------------------------------------------------------------------------
u   = spm_invFcdf(0.95,[1 (N - 2)]);
r   = sort(T);
r   = r(fix((1 - 0.05)*Ns));

% (linear) mapping between free energy difference and F ratio
%--------------------------------------------------------------------------
j   = abs(F) < 32;
b   = pinv([F(j) ones(size(F(j)))])*T(j);
Fq  = (-32:32)';
Tq  = [Fq, ones(size(Fq))]*b;

% show results
%==========================================================================
spm_figure('GetWin','Graphics');clf

subplot(2,2,1)
hist(F,32), hold on
xlabel('Log Bayes Factor','FontSize',16), ylabel('Frequency')
title('Null distribution','FontSize',16)
axis square

subplot(2,2,2)
hist(T,32), hold on
plot([u u],[0 Ns/4],'--r'), hold on
plot([r r],[0 Ns/4],'--b'), hold off
xlabel('Classical F-ratio','FontSize',16), ylabel('Frequency')
title('Null distribution','FontSize',16)
axis square

subplot(2,1,2)
plot(F,T,'.b','Markersize',8), hold on
plot(G,T,'.r','Markersize',8), hold on
plot(Fq,Tq,'b'), hold on
plot([3 3],[0 16],':r'), hold on
plot([0 0],[0 16],'--r'), hold on
plot([-32, 32],[r r],':b'), hold off
xlabel('free energy difference','FontSize',16), ylabel('Classical F-ratio')
title('Null distribution','FontSize',16)
axis([-8 8 0 16])
axis square


