function [y,w] = spm_nfm_mtf(P,M,U)
% Spectral response of a NFM (transfer function x noise spectrum)
% FORMAT [y,w] = spm_nfm_mtf(P,M,U)
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
 
% Dimitris Pinotsis, Karl Friston
% $Id: spm_nfm_mtf.m 4305 2011-04-12 18:15:32Z karl $
 
 
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
 
 
% spectrum of innovations (Gu) and noise (Gs and Gn)
%--------------------------------------------------------------------------
[Gu,Gs,Gn,w] = spm_csd_mtf_gu(P,M);
Gs           = Gs/8;
Gn           = Gn/8;
 
 
% trial-specific effects
%==========================================================================
try, X = U.X; catch, X = sparse(1,0); end
 
% cycle over trials
%--------------------------------------------------------------------------
for  t = 1:size(X,1)
    
    % baseline parameters
    %----------------------------------------------------------------------
    Q  = P;
    
    % trial-specific effective connectivity
    %----------------------------------------------------------------------
    for i = 1:size(X,2)
        
        
        % extrinsic connections
        %------------------------------------------------------------------
        for j = 1:length(A)
            Q.A{j} = Q.A{j} + X(t,i)*P.B{i};
        end
        
        % intrinsic connections
        %----------------------------------------------------------------------
        Q.H(:,1) = Q.H(:,1) + X(t,i)*diag(P.B{i});
        
    end
    
    % extrinsic coupling
    %----------------------------------------------------------------------
    E    = [32 16 4];
    A{1} = exp(Q.A{1})*E(1);
    A{2} = exp(Q.A{2})*E(2);
    A{3} = exp(Q.A{3})*E(3);
    C    = exp(Q.C);
    
    % intrinsic coupling
    %----------------------------------------------------------------------
    I    = [2 -8 2 1]*1000;
    A31  = exp(Q.G(:,1))*I(1);
    A32  = exp(Q.G(:,2))*I(2);
    A13  = exp(Q.G(:,3))*I(3);
    A23  = exp(Q.G(:,4))*I(4);
        
    % synaptic coupling
    %----------------------------------------------------------------------
    H    = [8 32];                  % receptor densities (excitatory, inhibitory)
    T    = [4 28]/1000;             % synaptic constants (excitatory, inhibitory)
    
    ke   = exp(Q.T(:,1))./T(1);     % excitatory rate constants
    me   = exp(Q.H(:,1)).*H(1);     % excitatory receptor density
    ki   = 1/T(2);                  % inhibitory rate constants
    mi   = H(2);                    % inhibitory receptor density
    g    = exp(P.R(1))/8;           % postsynaptic gain
    
    % spatial coupling
    %----------------------------------------------------------------------
    vel  = 16*exp(Q.vel)/1000;      % inverse velocity (transit time)
    c    = [1 2 1 1]*exp(Q.ext)/16; % extent (units of patch size)
    C23  = c(1);
    C13  = c(2);
    C32  = c(3);
    C31  = c(4);
    
  
    % analytic form of the spectrum
    %----------------------------------------------------------------------
    phi1  = 100;                    % gain of (Gaussian) lead field
    phi2  = 100;                    % precision of (Gaussian) lead field
    for m = 1:32
        
        k      = pi.*m;
        
        D23    = A23*(C23+1i*vel*w)./(C23^2-vel^2*w.^2+2*1i*vel*w+k^2);
        D13    = A13*(C13+1i*vel*w)./(C13^2-vel^2*w.^2+2*1i*vel*w+k^2);
        D32    = A32*(C32+1i*vel*w)./(C32^2-vel^2*w.^2+2*1i*vel*w+k^2);
        D31    = A31*(C31+1i*vel*w)./(C31^2-vel^2*w.^2+2*1i*vel*w+k^2);
        
        E3D    = D31*g*(ke^2)*(me^2).*((ki + 1i*w).^2);
        E3N    = (-D23.*D32.*(g^2)*ke*ki*me*mi).*((ke+1i*w).^2)+((ki+1i*w).^2).*(ke^4+4*ke^3*(1i*w)-4*1i*ke*(w).^3 + (w).^4 +...
                  -ke^2*(D13.*D31.*g^2*me^2+6*(w.^2)));
        
        E3     = pi*phi1*exp(-pi^2*k^2/phi2)*(E3D)./(E3N);
        G(:,m) = E3.*conj(E3);
        
    end
 
    % save trial-specific frequencies of interest
    %----------------------------------------------------------------------
    G    = M.U*sum(G,2);
    y{t} = full(G(If).*Gu(If));
    
end
 
% and add channel noise
%--------------------------------------------------------------------------
nc    = 1;
for t = 1:length(y)
    
    G     = y{t};
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
    y{t} = G(If,:);
    
end
 


