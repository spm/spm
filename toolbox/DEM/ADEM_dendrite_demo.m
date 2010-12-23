% This demo illustrates the use of Lotka-Volterra form SHCs (Stable
% heteroclinic channel) to prescribe active sampling (inference). In this
% example each (unstable) fixed point in the SHC attracts the agent to
% points on the circumference of a circle.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: ADEM_dendrite_demo.m 4146 2010-12-23 21:01:39Z karl $

% preliminaries
%==========================================================================
DEMO = 0;     % switch for demo
clear

% generative model - Stable heterclinic channel of pre-synaptic neurons
%==========================================================================
fx    = inline('spm_lotka_volterra(x,v,P)','x','v','P');
gx    = inline('x(P.w)','x','v','P');

% cpnnection weights
%--------------------------------------------------------------------------
np    = 5;                            % number of pre-synaptic neurons
ns    = 4;                            % number of denitric segments
ny    = ns*np;                        % number of synapses

w     = kron([1:np],ones(1,ns));
g     = rem(randperm(ns*np),np) + 1;  % g(1:ns) = 1; ###

P.w   = w';                           % synaptic weights
Q.w   = g';                           % synaptic weights


% level 1
%--------------------------------------------------------------------------
M(1).x  = sparse(1,1,np,np,1) - np;
M(1).f  = fx;
M(1).g  = gx;
M(1).pE = P;
M(1).Q  = spm_Ce(ones(1,ny));          % error components
M(1).hE = -2;                          % log-precision prior mean
M(1).hC = 2;                           % log-precision prior covaraince
M(1).W  = exp(16);                     % error precision

% level 2
%--------------------------------------------------------------------------
M(2).v  = 0;
M(2).V  = exp(16);


% generte presynaptic inputs
%--------------------------------------------------------------------------
M(1).E.nE = 4;
M(1).E.s  = 1;

N       = 128;
U       = sparse(1,N) + 1/2;
DEM     = spm_DEM_generate(M,U,Q,{2 16},{16});
DEM.U   = U;

spm_DEM_qU(DEM.pU)


% Synaptci pruning
%==========================================================================
if DEMO
    load DEM_Dendrites
else
    
    % Synaptci pruning
    %======================================================================
    for i = 1:64
        
        % integrate
        %------------------------------------------------------------------
        DEM   = spm_DEM(DEM);
        
        % Get priors and posteriors
        % -----------------------------------------------------------------
        qE    = DEM.qH.h{1};
        qC    = DEM.qH.C;
        pE    = DEM.M(1).hE;
        pC    = DEM.M(1).hC;
        
        
        % record
        %------------------------------------------------------------------
        G(:,i) = g';
        F(1,i) = DEM.F(end);
        H(:,i) = qE;
        
        % model search over new prior without the i-th parameter
        % -----------------------------------------------------------------
        for k = 1:ny
            R     = speye(ny,ny) - sparse(k,k,1,ny,ny);
            rE    = pE;
            rE(k) = 2;
            rC    = R*pC*R;
            Z(k)  = spm_log_evidence(qE,qC,pE,pC,rE,pC);
        end
        
        % change synatic locations based on optmised precision
        %------------------------------------------------------------------
        p     = (1 + tanh(Z))/2;
        j     = find(rand(size(p)) > p);
        j     = Z < 0;
        
        % and generate new presynatic inputs
        %------------------------------------------------------------------
        r     = rem(randperm(ns*np),np) + 1;
        g(j)  = r(j);
        Q.w   = g';
        GEN   = spm_DEM_generate(DEM.M,U,Q,{2 16},{8});
        DEM.Y = GEN.Y;
        
        
        % report
        %==================================================================
        spm_figure('Getwin','Figure 1');
        
        D     = sparse(1:ny,G(:,i),exp(H(:,i)),ny,np);
        subplot(2,1,1)
        imagesc(-D)
        axis image, axis off
        title('synaptic connections','FontSize',16)
        MM(i) = getframe(gca);
        
        subplot(2,2,3)
        plot(F - F(1))
        xlabel('iterations')
        axis square, box off
        title('Free-energy','FontSize',16)
        
        subplot(2,2,4)
        imagesc(-exp(H))
        xlabel('iterations')
        ylabel('synaptic strength')
        axis square
        title('Precision','FontSize',16)
        
    end
    
    % save 
    %---------------------------------------------------------------------
    save DEM_Dendrites

end

return

% Illustare slection of synatic connections
%==========================================================================
spm_figure('Getwin','Figure 1');

D     = sparse(1:ny,G(:,i),exp(H(:,i)),ny,np);
subplot(2,1,1)
imagesc(-D)
axis image, axis off
title('synaptic connections','FontSize',16)
MM(i) = getframe(gca);

subplot(2,2,1)
imagesc(-exp(H))
xlabel('iterations')
ylabel('synaptic strength')
axis square
title('Precision','FontSize',16)

subplot(2,2,2)
plot(F - F(1))
xlabel('iterations')
axis square, box off
title('Free-energy','FontSize',16)

T     = [2 4 32 64]:
for i = 1:length(T)
    
    subplot(2,4,3 + i)
    
end
    








% Illustare seuence specificity
%==========================================================================
SIM   = DEM;
U     = sparse(1,N) + 1;
SIM.M(1).E.nE = 1;
SIM.M(1).W    = 8;
SIM.M(1).hE   = 2;
SIM.M(2).v    = 1;
SIM.M(2).V    = 8;
n     = 4;
t     = 1:N;
j     = 1:np;
Q     = P;
for i = 1:n
    
    
    % Intgerate
    %----------------------------------------------------------------------
    spm_figure('Getwin','DEM');
    
    SIM   = spm_DEM_generate(SIM.M,U,Q,{2 16},{8});
    SIM   = spm_DEM(SIM);
    
    % plot
    %----------------------------------------------------------------------
    spm_figure('Getwin','Figure 2');
    
    subplot(n,2,(i - 1)*2 + 1)
    imagesc(SIM.Y)
    xlabel('time','FontSize',12)
    ylabel('activity','FontSize',12)
    title('pre-synaptic input','FontSize',16)
    
    subplot(n,2,(i - 1)*2 + 2)
    spm_plot_ci(t,SIM.qU.v{2},SIM.qU.C),hold on
    plot(t*0,'LineWidth',4), hold off
    xlabel('time','FontSize',12)
    ylabel('activity','FontSize',12)
    title('post-synaptic response','FontSize',16)
    axis([1 N -2 2])
    
    % change firing order
    %----------------------------------------------------------------------
    j     = randperm(np);
    Q.w   = kron(j,ones(1,ns));
    
    
end

