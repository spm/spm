function [PF] = spm_dcm_csd_source_optimise
% Stochastic optimisation of single source neural mass model
% FORMAT [PF] = spm_dcm_csd_source_optimise
%
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_csd_source_optimise.m 4348 2011-06-10 20:50:23Z karl $
 
 
% Initaislie fixed paramters
%==========================================================================
% E  = [1 1/2 1 1/2]*200;             % extrinsic (forward and backward)  
% G  = [4 4 4 4 4 4 4 4 2 1]*200;     % intrinsic connections
% D  = [1 16];                        % delays (intrinsic, extrinsic)
% T  = [2 2 16 28];                   % synaptic time constants
% R  = 1;                           % slope of sigmoid activation function

PF.G  = [4 4 4 4 4 4 4 4 2 1]*200;
PF.D  = [1 16];
PF.T  = [2 2 16 28];
HZ    = [60 16]; model = 'CMC'; s = [3 7]; Hz = [1:128]';


% PF.G  = [1 1 1 1/2]*64;         % intrinsic rates (g1 g2 g3 g4)
% PF.H  = [4 64];                 % receptor densities (excitatory, inhibitory)
% PF.T  = [4 8];                  % synaptic constants (excitatory, inhibitory)
% HZ    = [20 16]; model = 'SEP'; s = [7 9]; Hz = [1:64]';

N     = 512;
n     = length(spm_vec(PF));
sC    = speye(n,n)/32;
for k = 1:4
    for i = 1:N
        
        % parameters
        %----------------------------------------------------------------------
        pF     = spm_vec(PF);
        pF     = pF.*exp(sC*randn(n,1)/k);
        P(:,i) = pF;
        pF     = spm_unvec(pF,PF);
                
        % CSD
        %------------------------------------------------------------------
        G      = spm_dcm_csd_source_plot(model,s,pF,2*Hz(end));
        R(:,i) = spm_vec(G);
        
        % score - mean and variance over Hz
        %------------------------------------------------------------------
        for j = 1:length(s)
            p     = G(:,j,j);
            S(j)  = norm(p,1);
            p     = p/sum(p);
            E1(j) = sum(Hz.*p);
            E2(j) = sum((Hz - E1(j)).^2.*p);
            E3(j) = sum((Hz - E1(j)).^3.*p)/(E2(j)^(3/2));
            E4(j) = sum((Hz - E1(j)).^4.*p)/(E2(j)^2) - 3;
        end
        
        % entropy over sources
        %------------------------------------------------------------------
        S      = S/sum(S);
        S      = S*log(S');
        
        % cost function
        %------------------------------------------------------------------
        c(1,i) = sum((E1 - HZ).^2);
        c(2,i) = sum((sqrt(E2) - sqrt([64 64])).^2);
        c(3,i) = sum((E3 - [0 0]).^2);
        c(4,i) = sum((E4 - [0 0]).^2);
        c(5,i) = S*32;
        
    end
    
    % total cost and minimise
    %----------------------------------------------------------------------
    F     = sum(c)';
    i     = isfinite(F);
    F     = F(i); P = P(:,i); R = R(:,i);
    i     = find(F < (min(F) + 8*std(F)));
    F     = F(i); P = P(:,i); R = R(:,i);

    
    % Laplace approximation: p(P)
    %======================================================================

    % temperature
    %----------------------------------------------------------------------
    T     = std(F)/128;

    % mean
    %----------------------------------------------------------------------
    q     = exp(-(F - mean(F))/T);
    q     = q/sum(q);
    Lq    = P*q;

    % dispersion
    %----------------------------------------------------------------------
    for i = 1:n
        for j = 1:n
            C(i,j) = ((P(i,:) - Lq(i)).*(P(j,:) - Lq(j)))*q;
        end
    end   
    PF    = spm_unvec(Lq,PF);
    
    spm_dcm_csd_source_plot(model,s,PF,2*Hz(end));
    disp(PF)
    
end

return

% SEP models
%--------------------------------------------------------------------------
% E  = [32 16 4];              % extrinsic rates (forward, backward, lateral)
% G  = [1 1 1/4 1/4]*128;      % intrinsic rates (g1 g2 g3 g4)
% D  = [2 32];                 % propagation delays (intrinsic, extrinsic)
% H  = [4 32];                 % receptor densities (excitatory, inhibitory)
% T  = [4 8];                  % synaptic constants (excitatory, inhibitory)
% R  = [1 0];                  % parameters of static nonlinearity


% SVD of mapping from parameters to CSD
%--------------------------------------------------------------------------
X       = log(abs(R'));
X       = X(:,sum(X) > -8);
[U S V] = spm_svd(X'*P',0);

