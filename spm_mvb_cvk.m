function [p,pc,R2] = spm_mvb_cvk(MVB,k);
% Split-half cross validation of a multivariate Bayesian model
% FORMAT [p_value,percent,R2] = spm_mvb_cvk(MVB)
%
% MVB - Multivariate Bays structure
% k   - k-fold cross-validation
%
% p   - p_value: under a null GLM
% percent: proportion correct (median threshold)
% R2  - coeficient of determination
%
% spm_mvb_cvk performs a k-fold cross-validation by trying to predict
% the target variable using training and test partitions on orthogonal 
% mixtures of data (from null space of confounds)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mvb_cvk.m 1912 2008-07-11 18:02:03Z guillaume $
 
 
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
 
 
% create orthonormal projection to remove confounds
%--------------------------------------------------------------------------
Ns    = length(X);
X0    = orth(X0);
R     = speye(Ns) - X0*X0';
R     = orth(R);
X     = R'*X;
Y     = R'*Y;
V     = R'*R;
 
% leave-one-out
%--------------------------------------------------------------------------
Ns    = length(X);
qX    = zeros(Ns,1);
qE    = zeros(size(Y,2),Ns);
P     = zeros(size(Y,2),Ns);
if ~k, k = Ns - 1; end
 
% k-fold cross-validation
%==========================================================================
for i = 1:k
 
    % specify indices of training and test data
    %----------------------------------------------------------------------
    ns     = floor(Ns/k);
    test   = [1:ns] + (i - 1)*ns;
 
    % orthogonalise test and training partition
    %----------------------------------------------------------------------
    tran       = [1:Ns];
    tran(test) = [];
 
    % Training
    %======================================================================
    M        = spm_mvb(X(tran,:),Y(tran,:),[],U,[],MVB.Ni,MVB.sg);
 
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
 
 
% plot validation
%--------------------------------------------------------------------------
subplot(2,2,1)
s      = 1:length(pX);
plot(s,pX,s,qX,'-.')
xlabel('sample')
ylabel('response (adjusted)')
title('cross-validation')
axis square
 
subplot(2,2,2)
plot(pX,qX,'.')
xlabel('true')
ylabel('predicted')
title(sprintf('p-value (parametric) = %.5f',p))
axis square
 
 
% plot feature wietghts
%--------------------------------------------------------------------------
subplot(2,2,3)
imagesc(corrcoef(qE))
colorbar
xlabel('biparititon (k)')
title({'correlations among';'k-fold feature weights'})
axis square
 
subplot(2,2,4)
spm_mip(prod(P,2),MVB.XYZ(1:3,:),MVB.VOX)
title({[MVB.name ' (' MVB.contrast ')'];'prod( P(|weights| > 0) )'})
axis square
 
 
% display and assign in base memory
%--------------------------------------------------------------------------
fprintf('\np-value = %.4f; classification: %.1f%s; R-squared %.1f%s\n',p,pc,'%',R2,'%')
MVB.p_value = p;
MVB.percent = pc;
MVB.R2      = R2;
MVB.cvk     = struct('qX',qX,'qE',qE,'P',P);

 
assignin('base','MVB',MVB)
