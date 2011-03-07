function [y,w,s] = spm_csd_mtf(P,M,U)
% Spectral response of a NMM (transfer function x noise spectrum)
% FORMAT [y,w,s] = spm_csd_mtf(P,M,U)
%
% P - parameters
% M - neural mass model structure
% U - trial-specific effects
%
% G - {G(N,nc,nc}} - cross-spectral density for nc channels {trials}
%                  - for N frequencies in M.Hz [default 1:64Hz]
% w - frequencies
% s - directed transfer functions (complex)
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_csd_mtf.m 4232 2011-03-07 21:01:16Z karl $
 
 
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
W    = [0 1:N/2 -(N/2 - 1):-1]';         % frequencies of Fourier transform
W    = conj(2*pi*sqrt(-1)*W);            % radial frequencies
 
% number of channels and exogenous (neuronal) inputs
%--------------------------------------------------------------------------
nc   = size(spm_lx_erp(P,M),1);
nu   = M.m;
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
 
 
% trial-specific effects
%==========================================================================
try, X = U.X; catch, X = sparse(1,0); end
 
 
% cycle over trials
%--------------------------------------------------------------------------
for  c = 1:size(X,1)
 
    % baseline parameters
    %----------------------------------------------------------------------
    Q  = P;
 
    % trial-specific effective connectivity
    %----------------------------------------------------------------------
    for i = 1:size(X,2)
        Q.A{1} = Q.A{1} + X(c,i)*P.B{i};         % forward   connections
        Q.A{2} = Q.A{2} + X(c,i)*P.B{i};         % backward  connections
        Q.A{3} = Q.A{3} + X(c,i)*P.B{i};         % lateral   connections
       try
            Q.H = Q.H + X(c,i)*diag(P.B{i});     % intrinsic connections
        catch
            Q.G = Q.G + X(c,i)*diag(P.B{i});
        end
    end
 
    % get (extrinsic) delays and (intrinsic) delay operator
    %----------------------------------------------------------------------
    try
        Q.D         = diag(diag(P.D) + 32) - 32;
        [fx dfdx D] = feval(M.f,M.x,M.u,Q,M);
        d           = 16*exp(P.D)/1000;
        d           = d - diag(diag(d));
    catch
        D           = 1;                        % intrinsic delays
        d           = sparse(nc,nc);            % extrinsic delays
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
    S   = fft(K1);

 
    % [cross]-spectral density from neuronal innovations
    %----------------------------------------------------------------------
    G     = zeros(N/2,nc,nc);  
    for i = 1:nc
        for j = 1:nc
            for k = 1:nu
                Si       = S(:,i,k).*exp(W*d(i,k));
                Sj       = S(:,j,k).*exp(W*d(j,k));
                Gij      = Si.*conj(Sj);
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
sc     = norm(abs(spm_vec(y)))*4;
Gs    = Gs*sc;
Gn    = Gn*sc;
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

