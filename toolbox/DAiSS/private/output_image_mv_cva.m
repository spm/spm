function [chi,BIC,cva] = output_image_mv_cva(X,Y,c,U)
%% cva code to plug into bf_output_image_mv function
%function [chi,cva,t_stat] = bf_output_image_mv_cva(X,Y,c,U)
% CVA. See Chatfield and Collins

if nargin<4,
    U=[];
end;

if isempty(U),
    U=eye(size(Y,2));
end;

X0  = X - X*c*pinv(c);  %% make sure X0 is orthogonal to X
Xred   = full(X*c); %% reduced design matrix
X0  = spm_svd(X0); %% X0 is null space i.e. everything that is happening in other columns of X

%-Canonical Variates Analysis
%   ==========================================================================
% remove null space of contrast
%--------------------------------------------------------------------------
Y     = Y - X0*(X0'*Y); %% eg remove DC level or drift terms from all of Y
Xred     = Xred - X0*(X0'*Xred);
%% Xred=X*c-X0*(X0'*X*c)

P     = pinv(Xred);


[n,b] = size(Xred);
[n,m] = size(Y); %% n is number of epochs, m is number of features
b     = rank(Xred);
h     = min(b,m); %% either number of features or rank of X
f     = n - b - size(X0,2); %% number of epochs- rank(X) - number of confounds


% generalised eigensolution for treatment and residual sum of squares
%--------------------------------------------------------------------------
T     = Xred*(P*Y); %% predticon of Y based on X (P*Y is the coefficient)

SST   = T'*T;
SSR   = Y - T; %% residuals in Y (unexplained by X)
SSR   = SSR'*SSR;

[v,d] = eig(SSR\SST); %% explained/unexplained variance

[q,r] = sort(-real(diag(d)));
r     = r(1:h);
d     = real(d(r,r));
v     = real(v(:,r));
V     = U*v;                          % canonical vectors  (data)
v     = Y*v;                          % canonical variates (data)
W     = P*v;                          % canonical vectors  (design)
w     = Xred*W;                          % canonical variates (design)
%   C     = c*W;                          % canonical contrast (design)




% inference on dimensionality - p(i) test of D >= i; Wilk's Lambda := p(1)
%--------------------------------------------------------------------------

cval  = log(diag(d) + 1);
chi=[];df=[];p=[];p05thresh=[];
for i1 = 1:h
    chi(i1) = (f - (m - b + 1)/2)*sum(cval(i1:h));
    df(i1)  = (m - i1 + 1)*(b - i1 + 1); % m is the number of features, b is the rank of the design matrix
    p(i1)   = 1 - spm_Xcdf(chi(i1),df(i1));
    p05thresh(i1) = spm_invXcdf(1-0.05,df(i1));
end


ccorr=(W'*Xred'*Y*V)/sqrt((W'*Xred'*Xred*W)*(V'*Y'*Y*V));%% canonical correlation
ccorr=real(diag(ccorr));
BIC=n*log(1-ccorr.^2)+(b+m)*log(n); %http://www.sciencedirect.com/science/article/pii/S0167947311002660
% Comparison of penalty functions for sparse canonical correlation
% analysis,  Chalise , Brooke ,  Fridley: Computational Statistics and Data
% Analysis 2012, pp245-254
%% d=p+q, where x=N*p, Y=N*q

BIC=-0.5*BIC; %% to make this compatible with Model evidence estimates (see penny 2012)
cva.ccorr=ccorr;
cva.r=r;
cva.d=d;
cva.df=df;

cva.V=V;
cva.v=v;
cva.W=W;
cva.w=w;
%    cva.C=C;
cva.p05thresh=p05thresh;


end
