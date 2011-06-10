function [y,w,s] = spm_csd_mtf(P,M,U)
% Spectral response of a NMM (transfer function x noise spectrum)
% FORMAT [y,w,s] = spm_csd_mtf(P,M,U)
%
% P - parameters
% M - neural mass model structure
% U - trial-specific effects
%
% y - {y(N,nc,nc}} - cross-spectral density for nc channels {trials}
%                  - for N frequencies in M.Hz [default 1:64Hz]
% w - frequencies
% s - directed transfer functions (complex)
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_csd_mtf.m 4348 2011-06-10 20:50:23Z karl $


% compute log-spectral density
%==========================================================================

% frequencies of interest
%--------------------------------------------------------------------------
try
    dt = 1/(2*round(M.Hz(end)));
    N  = 1/dt;
    If = round(linspace(M.Hz(1),M.Hz(end),length(M.Hz)));
catch
    N  = 128;
    dt = 1/N;
    If = 1:N/2;
end
f    = [1:N/2]';                         % frequencies
w    = f(If);                            % frequencies selected

% number of channels and exogenous (neuronal) inputs
%--------------------------------------------------------------------------
nc   = M.l;
nu   = M.m;
nx   = size(M.x,2);
M.u  = sparse(nu,1);


% solve for fixed point (with 64 ms burn in)
%--------------------------------------------------------------------------
S    = M;
S.g  = {};
V.u  = sparse(8,M.m);
V.dt = 8/1000;
x    = spm_int_L(P,S,V);
x    = spm_unvec(x(end,:),S.x);
M.x  = x;

% spectrum of innovations (Gu) and noise (Gs and Gn)
%--------------------------------------------------------------------------
[Gu,Gs,Gn] = spm_csd_mtf_gu(P,M);

% get prior means (delays)
%--------------------------------------------------------------------------
try
    di = M.pF.D(1);                    % intrinsic delays
    de = M.pF.D(2);                    % extrinsic delays
catch
    de = 16;
    di = 1;
end


% trial-specific effects
%==========================================================================
try, X = U.X; catch, X = sparse(1,0); end
if isempty(X),       X = sparse(1,0); end


% cycle over trials
%--------------------------------------------------------------------------
for  c = 1:size(X,1)
    
    % baseline parameters
    %----------------------------------------------------------------------
    Q  = P;
    
    % trial-specific effective connectivity
    %----------------------------------------------------------------------
    for i = 1:size(X,2)
        
        
        % extrinsic connections
        %------------------------------------------------------------------
        for j = 1:length(A)
            Q.A{j} = Q.A{j} + X(c,i)*P.B{i};
        end
        
        % intrinsic connections
        %----------------------------------------------------------------------
        try
            Q.H(:,1) = Q.H(:,1) + X(c,i)*diag(P.B{i});
        catch
            Q.G(:,1) = Q.G(:,1) + X(c,i)*diag(P.B{i});
        end
        
    end
    
    
    % delays
    %----------------------------------------------------------------------
    try
        
        % evaluate delay matrix
        %------------------------------------------------------------------
        De = exp(P.D);
        Di = diag(diag(De));
        De = De - Di;
        De = De*de/1000;
        Di = Di*di/1000;
        De = kron(ones(nx,nx),De);
        Di = kron(ones(nx,nx) - speye(nx,nx),Di);
        D  = Di + De;
        
        % get delay operator
        %------------------------------------------------------------------
        D  = spm_dcm_delay(M,Q,D);
        
    catch
        D  = 1;
    end
    
    
    % augment and bi-linearise (with intrinsic delays)
    %----------------------------------------------------------------------
    [M0,M1,L] = spm_bireduce(M,Q,D);
    
    % project onto spatial modes
    %----------------------------------------------------------------------
    try,    L = M.U'*L;   end
    
    % kernels
    %----------------------------------------------------------------------
    [K0,K1]   = spm_kernels(M0,M1,L,N,dt);
    
    % Transfer functions (FFT of kernel)
    %----------------------------------------------------------------------
    S     = fft(K1);
    
    % [cross]-spectral density from neuronal innovations
    %----------------------------------------------------------------------
    G     = zeros(N/2,nc,nc);
    for i = 1:nc
        for j = 1:nc
            for k = 1:nu
                Gij      = S(:,i,k).*conj(S(:,j,k));
                Gij      = Gij((1:N/2) + 1).*Gu(:,k);
                G(:,i,j) = G(:,i,j) + Gij;
            end
        end
    end
    
    % save trial-specific frequencies of interest
    %----------------------------------------------------------------------
    y{c} = G(If,:,:);
    s{c} = S(If,:,:);
    
end

% and add channel noise
%--------------------------------------------------------------------------
for c = 1:length(y)
    
    G     = y{c};
    for i = 1:nc
        
        % channel specific noise
        %------------------------------------------------------------------
        G(:,i,i) = G(:,i,i) + Gs(If,i);
        
        % and cross-spectral density from common channel noise
        %------------------------------------------------------------------
        for j = 1:nc
            G(:,i,j) = G(:,i,j) + Gn(If);
        end
    end
    y{c} = G;
end

