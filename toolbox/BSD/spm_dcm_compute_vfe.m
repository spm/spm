function [F, L, Fp, Fy] = spm_dcm_compute_vfe(varargin)
% Usage: spm_dcm_compute_vfe(DCM)
% Usage: spm_dcm_compute_vfe(DCM, P)
% Usage: spm_dcm_compute_vfe(DCM, P, true)

DCM = varargin{1};

fast = false; 
if nargin > 1
    DCM.M.P = varargin{2}; 
    if nargin > 2
        fast = varargin{3};
    end
end


M = DCM.M; 
U = DCM.xU; 
Y = DCM.xY;
Cp = DCM.Cp; 

% check integrator or generation scheme
%--------------------------------------------------------------------------
try
    M.IS;
catch
    try
        M.IS = M.G;
    catch
        M.IS = 'spm_int';
    end
end

% check feature selection
%--------------------------------------------------------------------------
try
    M.FS;
catch
    M.FS = @(x)x;
end


% composition of feature selection and prediction (usually an integrator)
%--------------------------------------------------------------------------
try
    y  = Y.y;
catch
    y  = Y;
end

if isa(M.IS,'function_handle') && isa(M.FS,'function_handle')
    
    % components feature selection and generation functions
    %----------------------------------------------------------------------
    IS = @(P,M,U)M.FS(M.IS(P,M,U));
    y  = feval(M.FS,y);
    
else
    
    try
        
        % try FS(y,M)
        %------------------------------------------------------------------
        try
            y  = feval(M.FS,y,M);
            IS = inline([M.FS '(' M.IS '(P,M,U),M)'],'P','M','U');
            
            % try FS(y)
            %--------------------------------------------------------------
        catch
            y  = feval(M.FS,y);
            IS = inline([M.FS '(' M.IS '(P,M,U))'],'P','M','U');
        end
        
    catch
        
        % otherwise FS(y) = y
        %------------------------------------------------------------------
        try
            IS = inline([M.IS '(P,M,U)'],'P','M','U');
        catch
            IS = M.IS;
        end
    end
end

% converted to function handle
%--------------------------------------------------------------------------
IS = spm_funcheck(IS);

% parameter update eqation
%--------------------------------------------------------------------------
if isfield(M,'f'), M.f = spm_funcheck(M.f); end
if isfield(M,'g'), M.g = spm_funcheck(M.g); end
if isfield(M,'h'), M.h = spm_funcheck(M.h); end


% size of data (samples x response component x response component ...)
%--------------------------------------------------------------------------
if iscell(y)
    ns = size(y{1},1);
else
    ns = size(y,1);
end
ny     = length(spm_vec(y));           % total number of response variables
nr     = ny/ns;                        % number of response components
M.ns   = ns;                           % number of samples

% initial states
%--------------------------------------------------------------------------
try
    M.x;
catch
    if ~isfield(M,'n'), M.n = 0; end
    M.x = sparse(M.n,1);
end

% input
%--------------------------------------------------------------------------
try
    U;
catch
    U = [];
end

% initial parameters
%--------------------------------------------------------------------------
try
    spm_vec(M.P) - spm_vec(M.pE);
catch
    M.P = M.pE;
end

% time-step
%--------------------------------------------------------------------------
try
    dt = Y.dt;
catch
    dt = 1;
end

% precision components Q
%--------------------------------------------------------------------------
try
    Q = Y.Q;
    if isnumeric(Q), Q = {Q}; end
catch
    Q = spm_Ce(ns*ones(1,nr));
end
nh    = length(Q);                  % number of precision components
nq    = ny/length(Q{1});            % for compact Kronecker form of M-step


% prior moments (assume uninformative priors if not specified)
%--------------------------------------------------------------------------
pE       = M.pE;
try
    pC   = M.pC;
catch
    np   = spm_length(M.pE);
    pC   = speye(np,np)*exp(16);
end

% confounds (if specified)
%--------------------------------------------------------------------------
try
    nb   = size(Y.X0,1);            % number of bins
    nx   = ny/nb;                   % number of blocks
    dfdu = kron(speye(nx,nx),Y.X0);
catch
    dfdu = sparse(ny,0);
end
if isempty(dfdu), dfdu = sparse(ny,0); end


% hyperpriors - expectation (and initialize hyperparameters)
%--------------------------------------------------------------------------
try
    hE  = M.hE(:);
    if length(hE) ~= nh
        hE = hE + sparse(nh,1);
    end
catch
    hE  = sparse(nh,1) - log(var(spm_vec(y))) + 4;
end

try h =  -log(DCM.Ce); catch, h = hE; end

% hyperpriors - covariance
%--------------------------------------------------------------------------
try
    ihC = spm_inv(M.hC);
    if length(ihC) ~= nh
        ihC = ihC*speye(nh,nh);
    end
catch
    ihC = speye(nh,nh)*exp(4);
end



% unpack covariance
%--------------------------------------------------------------------------
if isstruct(pC)
    pC = spm_diag(spm_vec(pC));
end

% dimension reduction of parameter space
%--------------------------------------------------------------------------
V     = spm_svd(pC,0);
nu    = size(dfdu,2);                 % number of parameters (confounds)
np    = size(V,2);                    % number of parameters (effective)
ip    = (1:np)';
iu    = (1:nu)' + np;

% second-order moments (in reduced space)
%--------------------------------------------------------------------------
pC    = V'*pC*V;
Cp    = V'*Cp*V; 
uC    = speye(nu,nu)/1e-8;
ipC   = inv(spm_cat(spm_diag({pC,uC})));

% initialize conditional density
%--------------------------------------------------------------------------
Eu    = spm_pinv(dfdu)*spm_vec(y);
p     = [V'*(spm_vec(M.P) - spm_vec(M.pE)); Eu];
Ep    = spm_unvec(spm_vec(pE) + V*p(ip),pE);

if fast
    f = IS(Ep, M, U); 
else
    % gradients
    %------------------------------------------------------------------
    [dfdp,f] = spm_diff(IS,Ep,M,U,1,{V});
    dfdp     = reshape(spm_vec(dfdp),ny,np);
    
    J     = -[dfdp dfdu];
end

e     =  spm_vec(y) - spm_vec(f) - dfdu*p(iu);

% precision and conditional covariance
%------------------------------------------------------------------
iS    = sparse(0);
for i = 1:nh
    iS = iS + Q{i}*(exp(-32) + exp(h(i)));
end
iS    = kron(speye(nq),iS);

if fast
    d     = h     - hE;
    L(3) = 0; 
else
    S     = spm_inv(iS);
    Pp    = real(J'*iS*J);
    Cp    = spm_inv(Pp + ipC);
    
    % precision operators for M-Step
    %------------------------------------------------------------------
    for i = 1:nh
        P{i}   = Q{i}*exp(h(i));
        PS{i}  = P{i}*S;
        P{i}   = kron(speye(nq),P{i});
        JPJ{i} = real(J'*P{i}*J);
    end
    
    % derivatives: dLdh = dL/dh,...
    %------------------------------------------------------------------
    for i = 1:nh
        dFdh(i,1)      =   trace(PS{i})*nq/2 ...
            - real(e'*P{i}*e)/2 ...
            - spm_trace(Cp,JPJ{i})/2;
        for j = i:nh
            dFdhh(i,j) = - spm_trace(PS{i},PS{j})*nq/2;
            dFdhh(j,i) =   dFdhh(i,j);
        end
    end
    
    % add second order terms; noting diS/dh(i)h(i) = diS/dh(i) = P{i}
    %------------------------------------------------------------------
    dFdhh = dFdhh + diag(dFdh);
    
    % add hyperpriors
    %------------------------------------------------------------------
    d     = h     - hE;
    dFdh  = dFdh  - ihC*d;
    dFdhh = dFdhh - ihC;
    Ch    = spm_inv(real(-dFdhh));

    L(3) = spm_logdet(ihC*Ch)/2; 
end

A = spm_logdet(iS)*nq/2 ;
B =  - real(e'*iS*e)/2;
C =  - ny*log(8*atan(1))/2;
D = spm_logdet(ipC*Cp)/2;
E =  - p'*ipC*p/2; 
G = L(3);
H = - d'*ihC*d/2;
L(1) = A+B+C;
L(2) = D+E;
L(3) = G+H;
F    = sum(L);

spm_figure('GetWin', 'Free-energy insights');
subplot(3,2,1)
bar([B A C; E D 0; H G 0; L(1) L(2) L(3); F 0 0], 'stacked'); 
xticklabels({'Log likelihood', 'Parameters', 'Hyperparam.', 'Contrib.', 'VFE'})


subplot(3,2,2); 
Fy = [- real(e'*iS*e)/2; spm_logdet(iS)*nq/2; - ny*log(8*atan(1))/2]; 
bar(Fy)
title('Log likelihood'); 
xticklabels({'Accuracy', 'Precision', 'Normalisation'});

Fp = zeros(size(V, 1), 2); 
for i = 1:size(V, 1)
    Vp = V(i, :); 
    Kp = Vp' * Vp; 
    Fp(i,1) = - p'*Kp*ipC*Kp*p/2;
    Fp(i,2) = spm_logdet(V * ipC* Kp * Cp * V');
end
subplot(3,2,[3 4]); 
[~, Io] = sort(sum(Fp,2), 1);
bar(Fp(Io,:), 'stacked');
xticks(1:size(V,1)); 
xticklabels(spm_fieldindices(M.pE, Io));
legend({'distance', 'variance'})

dFdpE = zeros(size(V, 1), 1); 
dFdpC = zeros(size(V, 1), 1); 

pC  = spm_vec(M.pC); 
ipC = 1./pC; 
ipC(pC == 0) = 0;
p = spm_vec(M.P) - spm_vec(M.pE); 
dFdpE = tanh( p .* ipC/2);
dFdpC = tanh((p.*ipC).^2/2 - pC/2);


subplot(3,2,[5,6]); 
bar([dFdpE dFdpC], 'grouped');
xticks(1:size(V,1)); 
xticklabels(spm_fieldindices(M.pE, Io));
legend({'distance', 'precision'})
ylim([-1.1 1.1]);

return