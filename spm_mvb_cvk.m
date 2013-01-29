function [p,pc,R2] = spm_mvb_cvk(MVB,k)
% K-fold cross validation of a multivariate Bayesian model
% FORMAT [p_value,percent,R2] = spm_mvb_cvk(MVB,k)
%
% MVB - Multivariate Bayes structure
% k   - k-fold cross-validation ('0' implies a leave-one-out scheme)
%
% p   - p-value: under a null GLM
% percent: proportion correct (median threshold)
% R2  - coefficient of determination
%
% spm_mvb_cvk performs a k-fold cross-validation by trying to predict
% the target variable using training and test partitions on orthogonal 
% mixtures of data (from null space of confounds)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mvb_cvk.m 5219 2013-01-29 17:07:07Z spm $
 
 
%-partition order
%--------------------------------------------------------------------------
try
    k;
catch
    str   = 'k-fold cross-validation';
    k     = spm_input(str,'!+1','b',{'2','4','8','loo'},[2 4 8 0]);
end
 
%-Get figure handles and set title
%--------------------------------------------------------------------------
Fmvb = spm_figure('GetWin','MVB');
spm_clf(Fmvb);
 
% get MVB results
%==========================================================================
try
    MVB;
catch
    mvb  = spm_select(1,'mat','please select models',[],pwd,'MVB_*');
    MVB  = load(mvb(1,:));
    MVB  = MVB.MVB;
end
 
% check under null hypothesis
%--------------------------------------------------------------------------
% MVB.Y = randn(size(MVB.Y));
 
% whiten target and predictor (X) variables (Y) (i.e., remove correlations)
%--------------------------------------------------------------------------
K     = MVB.K;
X     = K*MVB.X;
Y     = K*MVB.Y;
X0    = K*MVB.X0;
U     = MVB.M.U;
Ni    = length(MVB.M.F) - 1;
 
 
% create orthonormal projection to remove confounds
%--------------------------------------------------------------------------
Ns    = length(X);
X0    = spm_svd(X0);
R     = speye(Ns) - X0*X0';
R     = spm_svd(R);
X     = R'*X;
Y     = R'*Y;
V     = R'*R;
 
% leave-one-out
%--------------------------------------------------------------------------
if ~k, k = length(X); end
Ns    = length(X);
qX    = zeros(Ns,1);
qE    = zeros(size(Y,2),k);
P     = zeros(size(Y,2),k);
 
% k-fold cross-validation
%==========================================================================
for i = 1:k
 
    % specify indices of training and test data
    %----------------------------------------------------------------------
    ns     = floor(Ns/k);
    test   = (1:ns) + (i - 1)*ns;
 
    % orthogonalise test and training partition
    %----------------------------------------------------------------------
    tran       = 1:Ns;
    tran(test) = [];
 
    % Training
    %======================================================================
    M        = spm_mvb(X(tran,:),Y(tran,:),[],U,[],Ni,MVB.sg);
 
    % Test
    %======================================================================
    qX(test) = qX(test) + Y(test,:)*M.qE;
    
    % record feature weights
    %----------------------------------------------------------------------
    qE(:,i)  = M.qE;
    
    % and posterior probabilities
    %----------------------------------------------------------------------
    P(:,i)   = 1 - spm_Ncdf(0,abs(M.qE),M.qC);
 
end
 
% parametric inference
%==========================================================================
 
% test correlation
%--------------------------------------------------------------------------
[T,df] = spm_ancova(X,V,qX,1);
p      = 1 - spm_Tcdf(T,df(2));
 
% percent correct (after projection)
%--------------------------------------------------------------------------
pX     = R*X;
qX     = R*qX;
T      = sign(pX - median(pX)) == sign(qX - median(qX));
pc     = 100*sum(T)/length(T);
R2     = corrcoef(pX,qX);
R2     = 100*(R2(1,2)^2);
 
% assign in base memory
%--------------------------------------------------------------------------
MVB.p_value = p;
MVB.percent = pc;
MVB.R2      = R2;
MVB.cvk     = struct('qX',qX,'qE',qE,'P',P);
 
% save results
%--------------------------------------------------------------------------
save(MVB.name,'MVB', spm_get_defaults('mat.format'))
assignin('base','MVB',MVB)

% display and plot validation
%--------------------------------------------------------------------------
spm_mvb_cvk_display(MVB)
