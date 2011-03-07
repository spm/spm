function [PE] = spm_dcm_csd_source_optimise
% Stochastic optimisation of single source neural mass model
% FORMAT [PE] = spm_dcm_csd_source_optimise
%
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_csd_source_optimise.m 4232 2011-03-07 21:01:16Z karl $
 
 
% Initaislie fixed paramters
%==========================================================================
% PF.G  = [1 1 -1 -1 1 1 1 1 -1 1]*512;
PF.T  = [2 2 10 5];
PF.D  = [2 16];
 
N     = 512;
s     = [3 7];
Hz    = [1:128]';
n     = length(spm_vec(PF));
sC    = speye(n,n)/2;
for k = 1:8
    for i = 1:N
        
        % parameters
        %----------------------------------------------------------------------
        pF     = spm_vec(PF);
        pF     = pF.*exp(sC*randn(n,1)/k);
        P(:,i) = pF;
        pF     = spm_unvec(pF,PF);
                
        % CSD
        %------------------------------------------------------------------
        G      = spm_dcm_csd_source_plot('CMC',s,pF);
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
        c(1,i) = sum((E1 - [60 16]).^2);
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
    T     = std(F)/64;

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
    
    spm_dcm_csd_source_plot('CMC',s,PF);
    disp(PF)
    
end

return

% SVD of mapping from parameters to CSD
%--------------------------------------------------------------------------
X       = log(abs(R'));
X       = X(:,sum(X) > -8);
[U S V] = spm_svd(X'*P',0);

