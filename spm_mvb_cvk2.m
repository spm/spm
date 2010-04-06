function [p,pc,R2] = spm_mvb_cvk2(MVB,k)
% k-fold cross validation of a multivariate Bayesian model
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
% mixtures of data (from null space of confounds).
% This version uses the optimised covariance model from spm_mvb.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mvb_cvk2.m 3806 2010-04-06 14:42:32Z ged $
 
 
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
 

% whiten target and predictor (X) variables (Y) (i.e., remove correlations)
%--------------------------------------------------------------------------
X     = MVB.X;
X0    = MVB.X0;
V     = MVB.V;
 
% residual forming matrix
%--------------------------------------------------------------------------
Ns    = length(X);
R     = speye(Ns) - X0*pinv(X0);
 
% leave-one-out
%--------------------------------------------------------------------------
if ~k
    k = Ns;
end
pX    = zeros(Ns,1);
qX    = zeros(Ns,1);
qE    = zeros(size(MVB.Y,2),k);
 
 
% k-fold cross-validation
%==========================================================================
for i = 1:k
    [px qx qe] = mvb_cv(MVB,i,k);
    pX         = pX + px;
    qX         = qX + qx;
    qE(:,i)    = qe;
end
 
 
% parametric inference
%==========================================================================
 
% test correlation
%--------------------------------------------------------------------------
[T df] = spm_ancova(X,V,qX,1);
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
MVB.cvk     = struct('qX',qX,'qE',qE);

% save results
%--------------------------------------------------------------------------
save(MVB.name,'MVB')
assignin('base','MVB',MVB)
 
% display and plot validation
%--------------------------------------------------------------------------
spm_mvb_cvk_display(MVB)
 
 
return
 
%==========================================================================
function [X,qX,qE] = mvb_cv(MVB,n,k)
%==========================================================================
% MVB - multivariate structure
% n   - subset
% k   - partition
 
% Unpack MVB and create test subspace
%--------------------------------------------------------------------------
V     = MVB.V;
U     = MVB.M.U;
X     = MVB.X;
Y     = MVB.Y;
X0    = MVB.X0;
h     = MVB.M.h;
Cp    = MVB.M.Cp;
 
% Specify indices of training and test data
%--------------------------------------------------------------------------
Ns    = length(X);
ns    = floor(Ns/k);
test  = [1:ns] + (n - 1)*ns;
tran  = [1:Ns];
tran(test) = [];
 
test  = full(sparse(test,test,1,Ns,Ns));
tran  = full(sparse(tran,tran,1,Ns,Ns)); 
 
% Training - add test space to confounds
%==========================================================================
R     = speye(Ns) - [X0 test]*pinv([X0 test]);
R     = spm_svd(R);
L     = R'*Y*U;
 
% get error covariance
%--------------------------------------------------------------------------
Ce    = sparse(Ns,Ns);
if isstruct(V)
    for i = 1:length(V)
        Ce = Ce + h(i)*V{i};
    end
else
    Ce = V*h(1);
end
Ce     = R'*Ce*R;
 
% MAP estimates of pattern weights from training data
%----------------------------------------------------------------------
MAP    = Cp*L'*inv(Ce + L*Cp*L');
qE     = MAP*R'*X;
 
% Test - add training space to confounds and get predicted X
%==========================================================================
R      = speye(Ns) - [X0 tran]*pinv([X0 tran]);
X      = R*X;                                              % test data
qE     = U*qE;                                             % weights
qX     = R*Y*qE;                                           % prediction
