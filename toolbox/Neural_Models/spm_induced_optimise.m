function spm_induced_optimise
% Demo routine that computes transfer functions for free parameters
%==========================================================================
%
% This an exploratory routine that computes the modulation transfer function
% for a range of parameters to enable the spectral responses to be optimised
% with respect to the model parameters of neural mass models.
%
% By editing the script, one can change the neuronal model or the hidden
% neuronal states that are characterised in terms of induced responses
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_induced_optimise.m 4866 2012-08-28 12:47:34Z karl $


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
pE      = spm_dcm_neural_priors({0 0 0},{},1,options.model);
P       = fieldnames(pE);
pE      = spm_L_priors(M.dipfit,pE);
pE      = spm_ssr_priors(pE);
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
M.m     = nu;
M.l     = Nc;

% solve for steady state
%--------------------------------------------------------------------------
M.x     = spm_dcm_neural_x(pE,M);


% Modulation transfer functions
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
    spm_figure('GetWin',sprintf('Panel %i',ifig));
    
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
                xlabel('(log) parameter scaling','FontSize',16)
                axis xy; drawnow
                
                % update graphics
                %------------------------------------------------------
                iplot     = iplot + 1;
                if iplot > 4
                    iplot = 1;
                    ifig  = ifig + 1;
                    spm_figure('GetWin',sprintf('Panel %i',ifig));
                end
                
            end
        end
    end
end
