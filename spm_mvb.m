function model = spm_mvb(X,Y,X0,U,V,nG)
% Multivariate Baysian inversion of a linear model
% FORMAT model = spm_mvb(X,Y,X0,U,V,nG)
% X      - contrast or target vector
% Y      - date feature matrix
% X0     - confounds
% xBF    - kernel structure or correlation matrix for errors
% nG     - number of Greedy iterations (nG = 1 => uniform hyperpiors)
%        - if not specifed, the search will terminate when F falls
%
% returns model:
%                F: log-evidence [F(0), F(1),...]
%               qE: conditional expectation of weights
%               qC: conditional variance of weights
%                h: covariance hyperparameters
%                G: covariance Partition indices
%               nk: number of patterns
%__________________________________________________________________________
%
% model: X = Y*P + X0*Q + R
%        P = U*E;           
%   cov(E) = h1*diag(G(:,1)) + h2*diag(G(:,2)) + ...
%__________________________________________________________________________

% defaults
%--------------------------------------------------------------------------
try, V;  catch, V  = []; end
try, nG; catch, nG = 8;  end

% get orders
%--------------------------------------------------------------------------
ns     = size(Y,1);                 % number of samples
nv     = size(Y,2);                 % unmber of voxels
nk     = size(U,2);                 % unmber of parameters

% confounds
%--------------------------------------------------------------------------
if ~length(X0); X0 = zeros(ns,1);  end
if ~length(U);  U  = sparse(nv,0); end

% null model
%--------------------------------------------------------------------------
model  = spm_mvb_G(X,Y,X0,sparse(nv,0),[],V);
if ~nk
    model.G = [];
    return
end

% Initialise G and greedy search
%==========================================================================
F  = model.F;
G  = ones(nk,1);
for  i = 1:nG

    % invert
    %----------------------------------------------------------------------
    M  = spm_mvb_G(X,Y,X0,U,G,V);
    
    % record free energy
    %----------------------------------------------------------------------
    F(end + 1)  = max([F M.F])
    
    % record condotioal esimtates
    %----------------------------------------------------------------------
    if M.F > F(i) || i == 1
        model   = M;
        model.G = G;
    end

    % create new spatial support
    %----------------------------------------------------------------------
    g            = find(G(:,end));
    ng           = ceil(length(g)/2);
    [q j]        = sort(-abs(M.qE(g)));
    q            = g(j(1:ng));
    G(q,end + 1) = 1;
    
    % break if cluster is one
    %----------------------------------------------------------------------
    if ng < 2, break, end

end
model.F  = F;
model.U  = U;

% project pattern weights to feature weights
%--------------------------------------------------------------------------
model.qE = U*model.qE;
model.qC = sum((U*model.qC).*U,2);


return

%==========================================================================
function model = spm_mvb_G(X,Y,X0,U,G,V);
%==========================================================================

% null space of confounds
%--------------------------------------------------------------------------
R      = orth(speye(size(X0,1)) - X0*pinv(X0));
Y      = R'*Y;
X      = R'*X;

% residual forming matrix
%--------------------------------------------------------------------------
ns     = size(Y,1);
nv     = size(Y,2);
nk     = size(U,2);

% random effects (and serial correlations)
%--------------------------------------------------------------------------
if length(V)
    try
        K  = convmtx(V.bf(1:V.T:end,1),size(X0,1))*V.T;
        K  = K*R;
        Qe = {speye(ns,ns) K'*K};
    catch
        Qe = {R'*V*R};
    end
else
    Qe = {speye(R'*R)};
end

% assemble emprical priors
%==========================================================================
Qp    = {};
LQp   = {};
LQpL  = {};
for i = 1:size(G,2)
    Qp{i}   = sparse(diag(G(:,i)));
    LQp{i}  = Y*U*Qp{i};
    LQpL{i} = LQp{i}*U'*Y';
end

% Inverse solution
%==========================================================================

% ReML
%--------------------------------------------------------------------------
if size(U,2)
    Q  = {Qe{:} LQpL{:}};
else
    Q  = Qe;
end
[Cy,h,P,F] = spm_reml_sc(X*X',[],Q,size(X,2));


% Covariances: sensor space - Ce and source space - L*Cp
%--------------------------------------------------------------------------
Ne    = length(Qe);
Np    = length(Qp);
L     = Y*U;
LCp   = sparse(ns,size(L,2));
Ce    = sparse(ns,ns);
Cp    = sparse(nk,nk);
he    = h([1:Ne]);
hp    = h([1:Np] + Ne);
for i = 1:Ne
    Ce = Ce + he(i)*Qe{i};
end
for i = 1:Np
    Cp  = Cp  + hp(i)*Qp{i};
    LCp = LCp + hp(i)*LQp{i};
end

% MAP estimates of instantaneous sources
%==========================================================================
iC    = inv(Cy);
MAP   = LCp'*iC;
qE    = MAP*X;

% conditional covariance (leading diagonal)
% Cq    = Cp - Cp*L'*iC*L*Cp;
%--------------------------------------------------------------------------
qC    = Cp - MAP*LCp;

% assemble results
%==========================================================================
model.F  = F;
model.qE = qE;
model.qC = qC;
model.h  = h;




