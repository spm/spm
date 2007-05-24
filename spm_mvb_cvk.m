function [p,percent] = spm_mvb_cvk(MVB);
% Split-half cross validation of a multivariate Bayesian model
% FORMAT [p_value,percent] = spm_mvb_cvk(MVB);
%   p_value: under a null GLM
%   percent: proportion correct
%
% spm_mvb_cv performs a two-fold cross-validation by trying to predict
% the target variable using a split halve test sample.
%__________________________________________________________________________
 
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

% k=fol cross validation
%==========================================================================
k     = 2;
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

% ReML esimate of non-sphericity
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
 
subplot(2,1,2)
plot(pX,qX,'.')
xlabel('true')
ylabel('predicted')
title(sprintf('p-value (parametric) = %.5f',p))
axis square

% displaye and assigin in base memory
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

% unpack MVB and create test subpsace
%--------------------------------------------------------------------------
V     = MVB.V;
U     = MVB.M.U;
G     = MVB.M.G;

% whitening matrix
%--------------------------------------------------------------------------
K     = MVB.K;
X     = K*MVB.X;
Y     = K*MVB.Y;
X0    = K*MVB.X0;

Ns    = length(X);
ns    = floor(Ns/k);
test  = [1:ns] + (n - 1)*ns;
tran  = [1:Ns];
tran(test) = [];

test  = full(sparse(test,test,1,Ns,Ns));
tran  = full(sparse(tran,tran,1,Ns,Ns)); 

% Training - add test space to confounds
%==========================================================================
M     = spm_mvb(X,Y,[X0 test],U,[],8);
 
% Test - add training space to confounds
%==========================================================================
R     = speye(Ns) - [X0 tran]*pinv([X0 tran]);
X     = R*X;
qX    = R*Y*U*M.qE;
Q     = test;

