function [F,df,xX,xCon,beta,V] = spm_reml_ancova(y,P,Fc)
% Classical inference for [hierarchial] linear models using ReML
% FORMAT [F,df,xX,xCon,beta,V] = spm_reml_ancova(y,P,Fc);
%
% y       - (n x 1)     response variable
% P{i}.X  - (n x m)     ith level design matrix i.e:
% P{i}.C  - {q}(n x n)  ith level contraints on the form of Cov{e{i}}
% Fc      - (m x q)     contrast matrix for the last level
%
% F     -  T or F values
% df    -  degrees of freedom
% beta  -  parameter estimates
% xX    -  design matrix structure
% xCon  -  contrast structure
%__________________________________________________________________________
%
% spm_ancova uses a General Linear Model of the form:
%
%                            y = X{1}*b{1} + e{1}
%                         b{1} = X{2}*b{2} + e{2}
%                                 ...
%
%                     b{n - 1} = X{n}*b{n} + e{n}
%
% e{n} ~ N{0,Ce{n}} 
%
% An F ratio is formed using OLS estimators of the parameters and ReML 
% estimators of the hyperparamters.
% If Fc has only one column a T statistic is returned,
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_reml_ancova.m 5219 2013-01-29 17:07:07Z spm $


% get ReML hyperparameter estimates
%--------------------------------------------------------------------------
[C,P] = spm_PEB(y,P);
L     = P{1}.h;         % ReML Hyperparemter estimates
W     = P{1}.W;         % ReML precisions = ddln(p(y|h)/dhdh

% check covariance constraints - assume i.i.d. errors conforming to X{i}
%--------------------------------------------------------------------------
for i = 1:length(P)
    if ~isfield(P{i},'C')
        [n,m] = size(P{i}.X);
        if i == 1
            P{i}.C        = {speye(n,n)};
        else
            for j = 1:m
                k         = find(P{i}.X(:,j));
                P{i}.C{j} = sparse(k,k,1,n,n);
            end
        end
    end
end

% form non-hierachical model
%--------------------------------------------------------------------------
X     = 1;
Q     = {};
for i = 1:length(P)
    
    % covariance components
    %----------------------------------------------------------------------
    for j = 1:length(P{i}.C)
        Q{end + 1} = X*P{i}.C{j}*X';
    end

    % design matrix
    %----------------------------------------------------------------------
    X   = X*P{i}.X;
end


% create design matrix structure and get pseudoinverse
%--------------------------------------------------------------------------
xX    = spm_sp('Set',X);
xX.pX = spm_sp('x-',xX);

% OLS parameter  estimates
%--------------------------------------------------------------------------
beta  = xX.pX*y;

% contrast
%--------------------------------------------------------------------------
xCon  = spm_FcUtil('Set','','F','c',Fc,xX);
h     = spm_FcUtil('Hsqr',xCon,xX);
X1o   = spm_FcUtil('X1o',xCon,xX);


% Note tr{MQ} = tr{h*xX.pX*Q{i}*xX.pX'*h'} because
% M = R0 - R = X1o*pinv(X1o) = xX.pX'*h'*h*xX.pX
%--------------------------------------------------------------------------
V     = sparse(0);
T     = sparse(1,length(Q));
for i = 1:length(Q);
    V    = V + L(i)*Q{i};
    T(i) = trace(h*xX.pX*Q{i}*xX.pX'*h');
end
V     = V*length(V)/trace(V);
TL    = T*L;


% degrees of freedom
%--------------------------------------------------------------------------
[trMV,trMVMV] = spm_SpUtil('trMV',X1o,V);
df            = [trMV^2/trMVMV (TL)^2/(T*inv(W)*T')];

if size(Fc,2) == 1

    % T statistics
    %----------------------------------------------------------------------
    F     = h*beta./sqrt(TL);
else
    % F statistics
    %----------------------------------------------------------------------
    F     = sum((h*beta).^2)./TL;

end 
