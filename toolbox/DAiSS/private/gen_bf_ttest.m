
function [t_stat,B,SE,F_stat] = gen_bf_ttest(X,Y,c,U)
% Run a t-test

if nargin<4,
    U=[];
end;
if isempty(U),
    U=eye(size(Y,2));
end;

X0  = X - X*c*pinv(c);  %% make sure X0 is orthogonal to X
Xred   = full(X*c); %% reduced design matrix
X0  = spm_svd(X0); %% X0 is null space i.e. everything that is happening in other columns of X

%   ==========================================================================
% remove null space of contrast
%--------------------------------------------------------------------------
Y     = Y - X0*(X0'*Y); %% eg remove DC level or drift terms from all of Y
Xred     = Xred - X0*(X0'*Xred);

P     = pinv(Xred);


[n,b] = size(Xred);
[n,m] = size(Y); %% n is number of epochs, m is number of features
b     = rank(Xred);
h     = min(b,m); %% either number of features or rank of X

Ym=mean(Y,2);
B  = pinv(Xred)*Ym;
RSS   = sum((Ym - Xred*B).^2);
MRSS  = RSS / (n-b);
SE    = sqrt(MRSS*(pinv(Xred'*Xred)));
t_stat=B./SE;
F_stat=(B./SE).^2;






end

