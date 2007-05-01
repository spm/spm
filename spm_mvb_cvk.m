function [p,percent] = spm_mvb_cv2(MVB);
% Split-half cross validation of a multivariate Bayesian model
% FORMAT [p_value,percent] = spm_mvb_cv2(MVB);
%   p_value: binocdf(K,N,1/2)
%   percent: proproation correct
%
% spm_mvb_cv performs a two-fold cross-validation by trying to predict
% the target variable using a split halve test sample.  To ensure the test
% and training features are independent, the target and prediction are
% decimated
%__________________________________________________________________________


% cycle over subsets
%==========================================================================
k     = 2;
X     = [];
qX    = [];
for i = 1:k
    [x qx] = mvb_cv2(MVB,i,k);
    X      = [X; x];
    qX     = [qX; qx];
end

% parameteric inference (with spaced smapling for correlations)
%==========================================================================
c        = 1 + size(MVB.X0,2);
c        = sparse(1,1,1,c,1);
i        = 1:8:length(X);
[T df]   = spm_ancova([qX(i) MVB.X0(i,:)],[],X(i),c);
p        = 1 - spm_Tcdf(T,df(2));

% percent correct
%--------------------------------------------------------------------------
T        = sign(X - median(X)) == sign(qX - median(qX));
percent  = 100*sum(T)/length(T);

% plot
%--------------------------------------------------------------------------
subplot(2,1,1)
s       = 1:length(X);
plot(s,X,s,qX,'-.')
xlabel('sample')
ylabel('response (decorrelated)')
title('Cross validation')
axis square

subplot(2,1,2)
plot(X,qX,'.')
xlabel('true')
ylabel('predicted')
title(sprintf('p-value (parametric) = %.5f',p))
axis square



function [X,qX] = mvb_cv2(MVB,n,k)
% MVB - multivariate structure
% n   - subset
% k   - partition

% Train
%==========================================================================

% unpack MVB
%--------------------------------------------------------------------------
Ns    = length(MVB.X);
ns    = floor(Ns/k);
test  = [1:ns] + (n - 1)*ns;
train = [1:Ns];
train(test) = [];

% remove confounds
%--------------------------------------------------------------------------
X0    = orth(MVB.X0);
R     = speye(size(X0,1)) - X0*inv(X0'*X0)*X0';
X     = R*MVB.X;
Y     = R*MVB.Y;

% de-correlation matrix
%--------------------------------------------------------------------------
try
    Qe = {MVB.V(train,train)};
catch
    Qe = {speye(length(train))};
end



% MAP estimate of voxel weights (qE)
%==========================================================================
U     = MVB.M.U;
G     = MVB.M.G;
H     = MVB.M.h;
ns    = length(train);
nv    = size(Y,2);
nk    = size(U,2);

% random effects (and serial correlations)
%--------------------------------------------------------------------------
Ne    = length(Qe);
Np    = length(H) - Ne;
he    = H([1:Ne]);
hp    = H([1:Np] + Ne);

% Covariances: sensor space - Ce and source space - L*Cp
%--------------------------------------------------------------------------
LCp   = sparse(ns,nk);
Ce    = sparse(ns,ns);
Cp    = sparse(nk,nk);
for i = 1:Ne
    Ce = Ce + he(i)*Qe{i};
end
for i = 1:Np
    Cp = Cp + hp(i)*sparse(diag(G(:,i)));
end
L     = Y(train,:)*U;
LCp   = L*Cp;
qE    = U'*LCp'*inv(LCp*L' + Ce)*X(train,:);


% Test
%==========================================================================

% remove confounds
%--------------------------------------------------------------------------
X     = X(test,:);
Y     = Y(test,:);
qX    = Y*qE;

return




