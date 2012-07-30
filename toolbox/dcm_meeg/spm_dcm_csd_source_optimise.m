function [PE] = spm_dcm_csd_source_optimise
% Stochastic optimisation of single source neural mass model
% FORMAT [PE] = spm_dcm_csd_source_optimise
%
% Edit the set up variable in the main body of this routine to specify 
% desired frequency responses (in selected populations)
%
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_csd_source_optimise.m 4814 2012-07-30 19:56:05Z karl $
 
 
% setup
%==========================================================================
HZ    = [48 16];                  % target frequencies
model = 'CMM';                    % model
s     = [2 4];                    % indices of hidden states
Hz    = (1:128)';                 % frequency range
param = {'G'};                    % parameters to vary
 
% priors
%==========================================================================
[PE PC] = spm_dcm_neural_priors({0 0 0},{},1,model);
 
i     = spm_fieldindices(PE,param{:});
n     = length(spm_vec(PE));
sC    = sparse(i,i,1/2,n,n);
 
 
% desired CSD for (2) sources (s)
%--------------------------------------------------------------------------
D     = [(exp(-(Hz - HZ(1)).^2/(2*64))*2); ...
         (exp(-(Hz - HZ(2)).^2/(2*32)))];
D     = D/max(D);
D     = D + 1/32;
D     = D/sum(D);
 
% optimise using KL divergence between desired and parameterised CSD
%==========================================================================
CF    = [];
N     = 256;
for k = linspace(1,8,8)
    for i = 1:N
        try
            
            % parameters
            %--------------------------------------------------------------
            pE     = spm_vec(PE);
            pE     = pE + sC*randn(n,1)/k;
            
            % CSD
            %--------------------------------------------------------------
            G  = spm_dcm_csd_source_plot(model,s,spm_unvec(pE,PE),2*Hz(end));
 
            % score - mean and variance over Hz
            %--------------------------------------------------------------
            G  = [G(:,1,1); G(:,2,2)];
            G  = G/max(G);
            G  = G + 1/32;
            G  = G/sum(D);
            
        end
        P(:,i) = pE;
        F(i,1) = sum(D.*log(D./G));
        
    end
    
    % (KL) cost and minimise
    %----------------------------------------------------------------------
    i     = isfinite(F);
    F     = F(i); P = P(:,i);
    i     = find(F < mean(exp(log(F))));
    F     = F(i); P = P(:,i);
    
    
    % Laplace approximation: p(P)
    %======================================================================
    
    % temperature
    %----------------------------------------------------------------------
    T     = std(F)/8;
    
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
    PE    = spm_unvec(Lq,PE);
    
    
    % plot update
    %----------------------------------------------------------------------
    spm_dcm_csd_source_plot(model,s,PE,2*Hz(end));
    disp(PE)
    
    % minimal cost
    %----------------------------------------------------------------------
    CF(end + 1) = min(F)
    
    subplot(2,2,3)
    plot(CF)
    xlabel('iteration')
    title('minimal cost function','FontSize',16)
    axis square
    
end
