function spm_induced_optimise(Ep,model)
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
% $Id: spm_induced_optimise.m 6937 2016-11-20 12:30:40Z karl $
 
 
% Model specification
%==========================================================================
if nargin < 2
    model = 'CMC';
end

% number of regions in coupled map lattice
%--------------------------------------------------------------------------
Nc    = 1;
Ns    = 1;
options.spatial  = 'LFP';
options.model    = upper(model);
M.dipfit.model = options.model;
M.dipfit.type  = options.spatial;
M.dipfit.Nc    = Nc;
M.dipfit.Ns    = Ns;

% within-trial effects: adjust onset relative to pst
%----------------------------------------------------------------------
M.ns   = 64;
M.ons  = 64;
M.dur  = 16;
U.dt   = 1/256;

 
 
% get priors
%--------------------------------------------------------------------------
pE  = spm_dcm_neural_priors({0 0 0},{},1,options.model);
pE  = spm_L_priors(M.dipfit,pE);
pE  = spm_ssr_priors(pE);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% specify hidden state
%--------------------------------------------------------------------------
% pE.J    = pE.J*0;
% pE.J(3) = 1;

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
P   = fieldnames(pE);
P   = {'G','T','S'};
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
% use input argument if specified
%--------------------------------------------------------------------------
if nargin, pE = Ep; end
 
 
% hidden neuronal states of interest
%--------------------------------------------------------------------------
[x,f]   = spm_dcm_x_neural(pE,options.model);
 
 
% orders and model
%==========================================================================
nx      = length(spm_vec(x));
nu      = size(pE.C,2);
u       = sparse(1,nu);
 
% create LFP model
%--------------------------------------------------------------------------
M.f     = f;
M.g     = 'spm_gx_erp';
M.x     = x;
M.n     = nx;
M.pE    = pE;
M.m     = nu;
M.l     = Nc;
 
% solve for steady state
%--------------------------------------------------------------------------
M.x     = spm_dcm_neural_x(pE,M);
 
 
% Dependency on parameters in terms of Modulation transfer functions
%==========================================================================
M.u     = u;
M.Hz    = 4:96;
 
% compute transfer functions for different parameters
%--------------------------------------------------------------------------
iplot = 1;
ifig  = 1;
D     = 2;
for k = 1:length(P)
    
    % check parameter exists
    %----------------------------------------------------------------------
    sfig = sprintf('%s: Parameter dependency - %i',model,ifig);
    spm_figure('GetWin',sfig);
    
    Q = pE.(P{k});
    
    if isnumeric(Q)
        for i = 1:size(Q,1)
            for j = 1:size(Q,2);
                
                % line search (with solution for steady state)
                %----------------------------------------------------------
                dQ    = linspace(Q(i,j) - D,Q(i,j) + D,32);
                for q = 1:length(dQ)
                    qE      = pE;
                    qE      = setfield(qE,P{k},{i,j},dQ(q));
                    [G,w]   = spm_csd_mtf(qE,M,[]);
                    GW(:,q) = G{1};
                end
                
                % plot
                %----------------------------------------------------------
                subplot(4,2,2*iplot - 1)
                plot(w,abs(GW))
                xlabel('frequency {Hz}')
                title(sprintf('Param: %s(%i,%i)',P{k},i,j),'FontSize',16)
                set(gca,'XLim',[0 w(end)]);
                
                
                subplot(4,2,2*iplot - 0)
                imagesc(dQ,w,log(abs(GW)))
                title('Transfer functions','FontSize',16)
                ylabel('Frequency')
                xlabel('(log) parameter scaling')
                axis xy; drawnow
                
                % update graphics
                %----------------------------------------------------------
                iplot     = iplot + 1;
                if iplot > 4
                    iplot = 1;
                    ifig  = ifig + 1;
                    sfig = sprintf('%s: Parameter dependency - %i',model,ifig);
                    spm_figure('GetWin',sfig);
                end
                
            end
        end
    end
end
 
return

% Dependency on hidden states in terms of Modulation transfer functions
%==========================================================================
 
% new figure
%--------------------------------------------------------------------------
iplot = 1;
ifig  = 1;
D     = 4;
M.Nm  = 3;
sfig  = sprintf('%s: State dependency - %i',model,ifig);
spm_figure('GetWin',sfig);

% Steady state solution and number of eigenmodes
%--------------------------------------------------------------------------
M.Nm  = 3;
x     = full(M.x);
 
 
% MTF, expanding around perturbed states
%==========================================================================

% evoked flucutations in hidden states
%--------------------------------------------------------------------------
G    = M;
G.g  = @(x,u,P,M) x;
erp  = spm_gen_erp(pE,G,U);
xmax = spm_unvec(max(erp{1}),x);
xmin = spm_unvec(min(erp{1}),x);
    
for i = 1:size(x,1)
    for j = 1:size(x,2);
        for k = 1:size(x,3);
            
            % line search
            %--------------------------------------------------------------
            dQ    = linspace(xmin(i,j,k),xmax(i,j,k),32);
            for q = 1:length(dQ)
                
                
                % MTF
                %----------------------------------------------------------
                qx        = x;
                qx(i,j,k) = qx(i,j,k) + dQ(q);
                M.x       = qx;
                [G w]     = spm_csd_mtf(pE,M);
                GW(:,q)   = G{1};
                
                % spectral decomposition
                %----------------------------------------------------------
                S         = spm_ssm2s(pE,M);
                S         = S(1:M.Nm);
                W(:,q)    = abs(imag(S)/(2*pi));
                A(:,q)    = min(4, (exp(real(S)) - 1)./(real(S)) );
                  
            end
            
            % plot
            %----------------------------------------------------------
            subplot(4,3,3*iplot - 2)
            plot(w,abs(GW))
            xlabel('frequency {Hz}')
            title(sprintf('Hidden State: (%i,%i,%i)',i,j,k),'FontSize',16)
            set(gca,'XLim',[0 w(end)]);            
            
            subplot(4,3,3*iplot - 1)
            imagesc(x(i,j,k) + dQ,w,log(abs(GW)))
            title('Transfer functions','FontSize',16)
            ylabel('Frequency')
            xlabel('Deviation')
            axis xy; drawnow
            
            subplot(4,3,3*iplot - 0)
            plot(W',log(A'),'.',W',log(A'),':','MarkerSize',16)
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
                sfig = sprintf('%s: State dependency - %i',model,ifig);
                spm_figure('GetWin',sfig);
            end
            
        end
    end
end
