function spm_induced_optimise
% Demo routine that computes transfer functions for free parameters
%==========================================================================
%
% This an exploratory routine that computes the modulation transfer function
% for a range of parameters and states to enable the spectral responses to 
% be optimised with respect to the model parameters of neural mass models 
% under different hidden states.
%
% By editing the script, one can change the neuronal model or the hidden
% neuronal states that are characterised in terms of induced responses
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_induced_optimise.m 4928 2012-09-14 21:40:18Z karl $
 
 
% Model specification
%==========================================================================
 
 
% number of regions in coupled map lattice
%--------------------------------------------------------------------------
Nc    = 1;
Ns    = 1;
options.spatial  = 'LFP';
options.model    = 'CMM';
M.dipfit.model = options.model;
M.dipfit.type  = options.spatial;
M.dipfit.Nc    = Nc;
M.dipfit.Ns    = Ns;
 
 
% get priors
%--------------------------------------------------------------------------
[pE pC] = spm_dcm_neural_priors({0 0 0},{},1,options.model);
P       = fieldnames(pE);
[pE pC] = spm_L_priors(M.dipfit,pE,pC);
[pE pC] = spm_ssr_priors(pE,pC);
[x,f]   = spm_dcm_x_neural(pE,options.model);
 
% hidden neuronal states of interest
%--------------------------------------------------------------------------
pE.J(1:4) = [0 1 0 0];
 
 
% orders and model
%==========================================================================
nx      = length(spm_vec(x ));
nu      = size(pE.C,2);
u       = sparse(1,nu);
 
% create LFP model
%--------------------------------------------------------------------------
M.f     = f;
M.g     = 'spm_gx_erp';
M.x     = x;
M.n     = nx;
M.pE    = pE;
M.pC    = pC;
M.m     = nu;
M.l     = Nc;
 
% solve for steady state
%--------------------------------------------------------------------------
M.x     = spm_dcm_neural_x(pE,M);
 
 
% Dependency on parameters in terms of Modulation transfer functions
%==========================================================================
M.u     = u;
M.Hz    = 4:64;
 
% compute transfer functions for different parameters
%--------------------------------------------------------------------------
iplot = 1;
ifig  = 1;
D     = 2;
for k = 1:length(P)
    
    % check parameter exists
    %----------------------------------------------------------------------
    spm_figure('GetWin',sprintf('Dependency on parameters %i',ifig));
    
    Q = getfield(pE,P{k});
    
    if isnumeric(Q)
        for i = 1:size(Q,1)
            for j = 1:size(Q,2);
                
                % line search
                %----------------------------------------------------------
                dQ    = linspace(Q(i,j) - D,Q(i,j) + D,32);
                for q = 1:length(dQ)
                    qE      = pE;
                    qE      = setfield(qE,P{k},{i,j},dQ(q));
                    [G w]   = spm_csd_mtf(qE,M);
                    GW(:,q) = G{1};
                end
                
                % plot
                %----------------------------------------------------------
                subplot(4,2,2*iplot - 1)
                plot(w,GW)
                xlabel('frequency {Hz}')
                title(sprintf('Param: %s(%i,%i)',P{k},i,j),'FontSize',16)
                
                
                subplot(4,2,2*iplot - 0)
                imagesc(dQ,w,log(abs(GW)))
                title('Transfer functions','FontSize',16)
                ylabel('Frequency')
                xlabel('(log) parameter scaling')
                axis xy; drawnow
                
                % update graphics
                %------------------------------------------------------
                iplot     = iplot + 1;
                if iplot > 4
                    iplot = 1;
                    ifig  = ifig + 1;
                    spm_figure('GetWin',sprintf('Dependency on parameters %i',ifig));
                end
                
            end
        end
    end
end
 
% Dependency on hidden states in terms of Modulation transfer functions
%==========================================================================

% new figure
%----------------------------------------------------------------------
spm_figure('GetWin',sprintf('Dependency on states %i',1));
iplot = 1;
ifig  = 1;
D     = 16;
for i = 1:size(x,1)
    for j = 1:size(x,2);
        for k = 1:size(x,3);
            
            % line search
            %--------------------------------------------------------------
            dQ    = linspace(-D,D,32);
            for q = 1:length(dQ)
                
                
                % MTF
                %----------------------------------------------------------
                qx        = x;
                qx(i,j,k) = qx(i,j,k) + dQ(q);
                QX(q)     = qx(i,j,k);
                M.x       = qx;
                [G w]     = spm_csd_mtf(pE,M);
                GW(:,q)   = G{1};
                
                % spectral decompostion
                %----------------------------------------------------------
                S        = spm_ssm2s(pE,M,[]);
                W(:,q)   = abs(imag(S)/(2*pi));
                A(:,q)   = min(4, (exp(real(S)) - 1)./(real(S)) );
                
                
            end
            
            % plot
            %----------------------------------------------------------
            subplot(4,3,3*iplot - 2)
            plot(w,GW)
            xlabel('frequency {Hz}')
            title(sprintf('Hidden State: (%i,%i,%i)',i,j,k),'FontSize',16)
            
            
            subplot(4,3,3*iplot - 1)
            imagesc(x(i,j,k) + dQ,w,log(abs(GW)))
            title('Transfer functions','FontSize',16)
            ylabel('Frequency')
            xlabel('(log) deviation')
            axis xy; drawnow
            
            subplot(4,3,3*iplot - 0)
            plot(W',log(A'),'o',W',log(A'),'k')
            title('Eigenmodes','FontSize',16)
            xlabel('Frequency')
            ylabel('log-power')
            axis xy; drawnow
            
            
            % update graphics
            %------------------------------------------------------
            iplot     = iplot + 1;
            if iplot > 4
                iplot = 1;
                ifig  = ifig + 1;
                spm_figure('GetWin',sprintf('Dependency on states %i',ifig));
            end
            
        end
    end
end

return


% Demonstration of how to optise a single source mean filed model in terms
% of its spectrial respsones
%==========================================================================

% Target spectrun (Y)
%--------------------------------------------------------------------------

H      = 2*pi*[48 24 12]';
Y      = -H/4 + 1i*H;
Y(4:6) = -4;

spm_figure('GetWin','Oprimzation');
subplot(2,1,1)
[g,w]  = spm_s2csd(Y);
plot(w,g)
title('Target spectral density','FontSize',16)
xlabel('Frequency')
ylabel('Power')

% create generative mdoel of the eigenspectrum (M) and invert
%--------------------------------------------------------------------------
M.IS  = 'spm_ssm2s'
for i = 1:32
    Ep   = spm_nlsi_GN(M,[],Y);
    M.pE = Ep;
end

% Show results with iotmised parameters
%--------------------------------------------------------------------------
spm_figure('GetWin','Oprimzation');
subplot(2,1,2)
[g,w]  = spm_s2csd(spm_ssm2s(Ep,M,[]));
plot(w,g)
title('Target spectral density','FontSize',16)
xlabel('Frequency')
ylabel('Power')







