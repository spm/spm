function [p,percent] = spm_mvb_cvk(MVB,k)
% Split-half cross validation of a multivariate Bayesian model
% FORMAT [p_value,percent] = spm_mvb_cvk(MVB,k);
%   p_value: under a null GLM
%   percent: proportion correct
%
% spm_mvb_cv performs a k-fold (def. k=2) cross-validation by trying to 
% predict the target variable using a split-in-k test sample.
%__________________________________________________________________________
 % Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_mvb_cvk.m 1131 2008-02-06 11:17:09Z spm $

%-Get figure handles and set title
%--------------------------------------------------------------------------
Fmvb = spm_figure('GetWin','MVB');
if isempty(Fmvb)
    Fmvb = spm_figure('Create','MVB','Multivariate Bayes');
    figure(Fmvb)
else
    figure(Fmvb);
    clf
end
 
% get MVB results
%==========================================================================
try
    MVB;
catch
    mvb  = spm_select(1,'mat','please select models',[],pwd,'MVB_*');
    MVB  = load(mvb(1,:));
    MVB  = MVB.MVB;
end
 
% k-fold cross validation
%==========================================================================
if nargin<2
    k     = 2;
end
pX    = 0;
qX    = 0;
for i = 1:k
    [px qx q] = mvb_cv(MVB,i,k);
    pX        = pX + px;
    qX        = qX + qx;
    Q{i}      = q;
end
 
% parametric inference
%==========================================================================
 
% ReML estimate of non-sphericity
%--------------------------------------------------------------------------
X       = [pX MVB.K*MVB.X0];
V       = spm_reml_sc(qX*qX',X,Q);
C       = sparse(1,1,1,size(X,2),1);
[T df]  = spm_ancova(X,V,qX,C);
p       = 1 - spm_Tcdf(T,df(2));
 
 
% percent correct (after smoothing)
%--------------------------------------------------------------------------
S       = inv(MVB.K);
pX      = S*pX;
qX      = S*qX;
T       = sign(pX - median(pX)) == sign(qX - median(qX));
percent = 100*sum(T)/length(T);
 
% plot
%--------------------------------------------------------------------------
subplot(2,1,1)
s       = 1:length(pX);
plot(s,pX,s,qX,'-.')
xlabel('sample')
ylabel('response (adjusted)')
title('Cross validation')
axis square
legend('true','predicted')
 
subplot(2,1,2)
plot(pX,qX,'.')
xlabel('true')
ylabel('predicted')
title(sprintf('p-value (parametric) = %.5f',p))
axis square
 
% display and assign in base memory
%--------------------------------------------------------------------------
fprintf('\np-value = %.4f; percent: %.1f\n',p,percent)
MVB.p_value = p;
MVB.percent = percent;
assignin('base','MVB',MVB)
 
return
 
%==========================================================================
function [X,qX,Q] = mvb_cv(MVB,n,k)
%==========================================================================
% MVB - multivariate structure
% n   - subset
% k   - partition
 
% unpack MVB and create test subspace
%--------------------------------------------------------------------------
V     = MVB.V;
U     = MVB.M.U;
 
% whitening matrix
%--------------------------------------------------------------------------
K     = MVB.K;
X     = K*MVB.X;
Y     = K*MVB.Y;
X0    = K*MVB.X0;
 
% specify indices of training and test data
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
R      = speye(Ns) - [X0 test]*pinv([X0 test]);
Qe     = speye(Ns);
Qp     = MVB.M.Cp;
L      = R*Y*U;
Q      = {Qe L*Qp*L'};
 
% re-estimate covariance components
%----------------------------------------------------------------------
[Cy,h] = spm_reml_sc(R*X*X'*R',[X0 test],Q,size(X,2));
 
% MAP estimates of pattern weights
%----------------------------------------------------------------------
MAP    = h(2)*Qp*L'*R'*inv(Cy)*R;
qE     = MAP*X;
 
% Test - add training space to confounds
%==========================================================================
R      = speye(Ns) - [X0 tran]*pinv([X0 tran]);
X      = R*X;                                              % test data
qX     = R*Y*U*qE;                                         % prediction
Q      = test;                                             % test indices
