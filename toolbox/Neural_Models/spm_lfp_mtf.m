function [y,w] = spm_lfp_mtf(P,M,U)
% Spectral response of a NMM (transfer function x noise spectrum)
% FORMAT [G,w] = spm_lfp_mtf(P,M,U)
%
% P - parameters
% M - neural mass model structure
% U - trial-specific effects
%
% G - {G(N,nc,nc}} - cross-spectral density for nc channels {trials}
%                  - for N frequencies in M.Hz [default 1:64Hz]
% w - frequencies
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_lfp_mtf.m 3119 2009-05-13 10:53:31Z rosalyn $


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
f    = [1:N/2]';
w    = f(If);

% exogenous (neuronal) inputs
%--------------------------------------------------------------------------
M.u  = sparse(M.m,1);

% solve for fixed point (i.e., 64ms burn in)
%--------------------------------------------------------------------------
S    = M;
S.g  = {};
V.u  = sparse(8,M.m);
V.dt = 8/1000;
x    = spm_int_L(P,S,V);
x    = spm_unvec(x(end,:),S.x);
M.x  = x;

% get delay operator
%--------------------------------------------------------------------------
try
    [fx dfdx D] = feval(M.f,M.x,M.u,P,M);
catch
    D = 1;
end

% spectrum of innovations (Gu)
%--------------------------------------------------------------------------
Gu   = exp(P.a)*f.^(-1)*2;                % spectral density of (AR) input
Gu   = Gu + exp(P.b);                     % spectral density of IID input

% channel noise (specific and non-specific)
%--------------------------------------------------------------------------
Gs   = (exp(P.c(1))*f.^(-1) + exp(P.d(1)))/8;
Gn   = (exp(P.c(2))*f.^(-1) + exp(P.d(2)))/16;


% trial-specific effects
%==========================================================================
try, X = U.X; catch, X = sparse(1,0); end


% cycle over trials
%--------------------------------------------------------------------------
for  c = 1:size(X,1)

    % basline parameters
    %----------------------------------------------------------------------
    Q  = P;

    % trial-specific effective connectivity
    %----------------------------------------------------------------------
    for i = 1:size(X,2)
        Q.A{1} = Q.A{1} + X(c,i)*P.B{i};         % forward   connections
        Q.A{2} = Q.A{2} + X(c,i)*P.B{i};         % backward  connections
        Q.A{3} = Q.A{3} + X(c,i)*P.B{i};         % lateral   connections
       try
            Q.H = Q.H + X(c,i)*diag(P.B{i});      % intrinsic connections
        catch
            Q.G = Q.G + X(c,i)*diag(P.B{i});
        end
    end

    % augment and bi-linearise (with delays)
    %----------------------------------------------------------------------
    [M0,M1,L] = spm_bireduce(M,Q,D);

    % project onto spatial modes
    %----------------------------------------------------------------------
    try
        L = M.U'*L;
    end

    % compute modulation transfer function using FFT of the kernels
    %----------------------------------------------------------------------
    [K0,K1]   = spm_kernels(M0,M1,L,N,dt);


    % [cross]-spectral density
    %----------------------------------------------------------------------
    [N,nc,nu] = size(K1);
    G     = zeros(N/2,nc,nc);
    for i = 1:nc
        for j = i:nc

            % cross-spectral density from neuronal interactions
            %--------------------------------------------------------------
            for k = 1:nu
                Si       = fft(K1(:,i,k));
                Sj       = fft(K1(:,j,k));
                Gij      = Si.*conj(Sj);
                Gij      = abs(Gij([1:N/2] + 1)).*Gu;
                G(:,i,j) = G(:,i,j) + Gij;
            end

            % cross-spectral density from channel noise
            %--------------------------------------------------------------
            G(:,i,j) = G(:,i,j) + Gn;            % common noise
            if i == j
                G(:,i,i) = G(:,i,i) + Gs;        % and channel specific

            else
                % fill in lower half of CSD matrix
                %------------------------------------------------------
                G(:,j,i) = G(:,i,j);
            end
        end
    end

    % save frequencies of interest
    %----------------------------------------------------------------------
    y{c} = G(If,:,:);

end

