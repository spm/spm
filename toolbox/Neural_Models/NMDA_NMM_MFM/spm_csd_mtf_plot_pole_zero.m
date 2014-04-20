function [b,a] = spm_csd_mtf_plot_pole_zero(P,M,U,region_stab)
% Spectral response of a NMM (transfer function x noise spectrum)
% FORMAT [b,a] = spm_csd_mtf_plot_pole_zero(P,M,U,region_stab)
%
% P - parameters
% M - neural mass model structure
% U - trial-specific effects
% regions stab: which region in the DCM (per source list) to examine
% stability
%
% Returns poles and zeros and plots them 
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Rosalyn Moran
% $Id: spm_csd_mtf_plot_pole_zero.m 5964 2014-04-20 09:48:58Z karl $


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
        for j = 1:length(Q.A)
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
    
    % augment and bi-linearise
    %----------------------------------------------------------------------
    M.D       = 1;
    [M0,M1,L] = spm_bireduce(M,Q);
    
    % project onto spatial modes
    %----------------------------------------------------------------------
    try, L = M.U'*L; end
    
    % kernels
    %----------------------------------------------------------------------
    [~,K1] = spm_kernels(M0,M1,L,N,dt);
    
    % State Space Representation
    %----------------------------------------------------------------------
    A = full(M0(2:end,2:end));
    for i = 1:nc
        B(:,i) = full(M1{i}(2:end,1));
    end
    
    C     = full(L(:,2:end));
    D     = zeros(nc,nc);
    [b,a] = spm_ss2tf(A,B,C(region_stab,:),D(region_stab,:),1);
    zplane(b,a)    

end


