function DEM_demo_dendrite
% Free-energy and the single neuron:
%__________________________________________________________________________
% This demo illustrates the use of Lotka-Volterra form SHCs (Stable
% heteroclinic channels) to prescribe active sampling (inference). In this
% example, we assume that neurons self-organise to minimise a free energy
% bound on the informational surprise in the pre-synaptic inputs that are 
% sampled.  We model this as a selective pruning of post-synaptic spines 
% that are expressed on the dendritic tree.  This pruning occurs when the 
% (optimised) post-synaptic gain falls to small values.  Crucially, post-
% synaptic gain (encoding the precision of the neuron’s prediction errors 
% about its pre-synaptic inputs) is itself optimised with respect to free-
% energy.  Furthermore, the pruning itself suppresses free-energy as the 
% neuron selects post-synaptic specialisations that conform to its 
% expectations.  This provide a principled account of how neurons organise 
% and selectively sample the myriad of potential pre-synaptic inputs they 
% are exposed to, but it also connects elemental neuronal (dendritic) 
% processing to generic schemes in statistics and machine learning:
% such as Bayesian model selection and automatic relevance determination. 
% The demonstration of this scheme simulates direction selectivity in post
% synaptic transients and (see notes after 'return') spike-timing dependent
% plasticity.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_dendrite.m 6290 2014-12-20 22:11:50Z karl $
 
% preliminaries
%==========================================================================
rng('default')
DemoMode = 0;
 
% generative model - Stable heteroclinic channel of pre-synaptic neurons
%==========================================================================
fx  = inline('spm_lotka_volterra(x,exp(v))','x','v','P');
gx  = inline('x(P.w)','x','v','P');
 
% connection weights: NB synapses are changed by switching Y (not P.w)
%--------------------------------------------------------------------------
np  = 5;                              % number of pre-synaptic neurons
ns  = 4;                              % number of dendritic segments
ny  = ns*np;                          % number of synapses

w   = kron((1:np),ones(1,ns));        % prior contacts
g   = rem(randperm(ns*np),np) + 1;    % initial contacts
 
P.w = w';                             % synaptic weights 
Q.w = g';                             % synaptic weights
 
 
% level 1
%--------------------------------------------------------------------------
M(1).x  = sparse(1,1,np,np,1) - np;
M(1).f  = fx;
M(1).g  = gx;
M(1).pE = P;
M(1).Q  = spm_Ce(ones(1,ny));          % error components
M(1).hE = -2;                          % log-precision prior mean
M(1).hC = 2;                           % log-precision prior covariance
M(1).W  = exp(16);                     % error precision
 
% level 2
%--------------------------------------------------------------------------
M(2).v  = 0;
M(2).V  = exp(16);
 
 
% generate presynaptic inputs
%--------------------------------------------------------------------------
M(1).E.nE = 4;
M(1).E.s  = 1;
 
N       = 128;
U       = log(sparse(1,N) + 1/2);
DEM     = spm_DEM_generate(M,U,Q,{2 16},{8});
DEM.U   = U; 
 
% Synaptic pruning
%==========================================================================
if DemoMode
    
    load DEM_Dendrites
    
else
    
    % Synaptic pruning
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
        G(:,i) = Q.w;
        F(1,i) = DEM.F(end);
        H(:,i) = qE;
        
        % model search over new prior without the i-th parameter
        % -----------------------------------------------------------------
        for k = 1:ny
            rE    = pE;
            rE(k) = 4;
            Z(k)  = spm_log_evidence(qE,qC,pE,pC,rE,pC);
        end
        
        % change synaptic locations based on optimized precision
        %------------------------------------------------------------------
        j     = Z < 0;
        
        % and generate new presynaptic inputs
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
        drawnow
        
    end
    
    % save 
    %----------------------------------------------------------------------
    save DEM_Dendrites
 
end
 
 
% Illustrate selection of synaptic connections
%==========================================================================
spm_figure('Getwin','Figure 1');
 
subplot(2,2,1)
imagesc(-exp(H))
xlabel('cycles')
ylabel('synaptic strength')
axis square
title('Precision','FontSize',16)
 
subplot(2,2,2)
plot(F - F(1))
xlabel('cycles')
axis square, box off
title('Free-energy','FontSize',16)
spm_axis tight
 
% show connectivity at selected time points
%--------------------------------------------------------------------------
T     = [1 2 32 48];
for i = 1:length(T)
    subplot(2,4,4 + i)
    imagesc(-sparse(1:ny,G(:,T(i)),exp(H(:,T(i))),ny,np))
    axis image
    title(sprintf('Cycle %i',T(i)),'FontSize',16)
end
    
 
% Illustrate synaptic selection function (cf post hoc model selection)
%==========================================================================
spm_figure('Getwin','Figure 2');
 
% Get priors and posteriors
% -------------------------------------------------------------------------
qE    = DEM.qH.h{1};
qC    = DEM.qH.C;
pE    = DEM.M(1).hE;
pC    = DEM.M(1).hC;
 
% model search over new prior without the i-th parameter
% -------------------------------------------------------------------------
QE    = linspace(-2,4,32);
for k = 1:length(QE)
    rE    = pE;
    rE(1) = 4;
    qE(1) = QE(k);
    L     = spm_log_evidence(qE,qC,pE,pC,rE,pC);
    p(k)  = exp(L)/(exp(L) + 1);
end
 
spm_figure('Getwin','Figure 2');
subplot(2,1,1)
plot(QE,p,QE,QE*0 + 1/2,':');
ylabel('conditional probability')
xlabel('synaptic strength (log-precision)')
title('Synaptic selection function','FontSize',16)
axis square
 
 
% make move of synaptic re-organization
%--------------------------------------------------------------------------
subplot(2,1,2),cla
for i = 1:size(G,2)
    imagesc(-sparse(1:ny,G(:,i),exp(H(:,i)),ny,np))
    axis image
    MM(i) = getframe(gca);
end
 
title('synaptic connections','FontSize',16)
xlabel('(left) click image for movie','FontSize',16)
 
% set ButtonDownFcn
%--------------------------------------------------------------------------
h = findobj(gca,'type','image');
set(h(1),'Userdata',{MM,16})
set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')
 
return


% Illustrate sequence specificity
%==========================================================================
SIM   = DEM;
U     = -((1:N) - N/2).^2/(2*(N/8)^2);
SIM.M(1).E.nE =  1;
SIM.M(1).W    =  8;
SIM.M(1).hE   =  2;
SIM.M(2).v    = -1;
SIM.M(2).V    =  8;
n     = 3;
t     = 1:N;
Q     = P;
for i = 1:n
    
    
    % Integrate
    %----------------------------------------------------------------------  
    SIM   = spm_DEM_generate(SIM.M,U,Q,{2 16},{8});
    SIM   = spm_DEM(SIM);
    
    % plot
    %----------------------------------------------------------------------
    spm_figure('Getwin','Figure 3');
    
    subplot(n,2,(i - 1)*2 + 1)
    imagesc(SIM.Y)
    xlabel('time','FontSize',12)
    ylabel('synapse','FontSize',12)
    title('presynaptic input','FontSize',16)
    
    subplot(n,2,(i - 1)*2 + 2)
   
    m    = SIM.qU.v{2};
    ci   = 1.64* sqrt(spm_vec(SIM.qU.C)');
    cu   = exp(m + ci);
    cl   = exp(m - ci);
    m    = exp(m);
    fill([t fliplr(t)],[full(cu) fliplr(full(cl))],...
         [1 1 1]*.8,'EdgeColor',[1 1 1]*.8),hold on
    plot(t,m), hold off
    xlabel('time','FontSize',12)
    ylabel('activity','FontSize',12)
    title('postsynaptic response','FontSize',16)
    axis([1 N 0 2])
    
    % change presynaptic firing order
    %----------------------------------------------------------------------
    j     = randperm(np);
    Q.w   = kron(j,ones(1,ns));
    
end
 
% show response
%==========================================================================
spm_figure('Getwin','DEM');
spm_DEM_qU(SIM.qU,SIM.pU)
 

% Changing the speed (u) using IN and OUT sequences
%==========================================================================
SIM           = DEM;
 
SIM.M(1).E.nE =  1;
SIM.M(1).W    =  8;
SIM.M(1).hE   =  2;
SIM.M(2).v    = -1;
SIM.M(2).V    =  8;
 
IN    = 1:np;
OUT   = fliplr(IN);
u     = [1 2 4 6];
for i = 1:length(u)
    
    % Speed
    %----------------------------------------------------------------------
    U      = log(zeros(1,round(64/u(i))) + u(i));
        
    % Integrate IN
    %---------------------------------------------------------------------- 
    Q.w    = kron(IN,ones(1,ns));
    SIM    = spm_DEM_generate(SIM.M,U,Q,{2 16},{8});
    SIM    = spm_DEM(SIM);
    V(1,i) = mean(exp(SIM.qU.v{2}));
    
    % Integrate OUT
    %---------------------------------------------------------------------- 
    Q.w    = kron(OUT,ones(1,ns));
    SIM    = spm_DEM_generate(SIM.M,U,Q,{2 16},{8});
    SIM    = spm_DEM(SIM);
    V(2,i) = mean(exp(SIM.qU.v{2}));
    
end
 
% plot direction-selectivity
%--------------------------------------------------------------------------
spm_figure('Getwin','Figure 4');
 
subplot(2,1,1)
plot(u,V(1,:),'r',u,V(1,:),'or',u,V(2,:),'b',u,V(2,:),'ob')
xlabel('input velocity','FontSize',12)
ylabel('mean response','FontSize',12)
title('Velocity-dependent responses','FontSize',16)
axis square, box off

 
return


% Spike-timing dependent plasticity
%==========================================================================
U             = -((1:N) - N/2).^2/(2*(N/8)^2);
SIM           = spm_DEM_generate(DEM.M,U,P,{2 16},{8});
SIM.M(1).E.nE =  4;
SIM.M(1).W    =  8;
SIM.M(1).hE   =  0;
SIM.M(2).v    = -1;
SIM.M(2).V    =  8;
 
% time delay and synapse
%--------------------------------------------------------------------------
delay = -32:8:32;                   
j     = 10;
Y     = SIM.Y;
for i = 1:length(delay)
    
    % delay j-th synaptic input
    %----------------------------------------------------------------------
    y          = Y(j,:);
    y          = [(zeros(1,128) + y(1)) y (zeros(1,128) + y(end))];
    y          = y((1:N) + 128 - delay(i));
    SIM.Y(j,:) = y;
    
    % optimise precision
    %----------------------------------------------------------------------
    SIM        = spm_DEM(SIM);
    
    % record plasticity
    %----------------------------------------------------------------------
    qH(i)      = SIM.qH.h{1}(j);
        
end
 
% plot Spike-timing dependent plasticity
%--------------------------------------------------------------------------
spm_figure('Getwin','Figure 4');
 
subplot(2,1,2)
plot(delay,qH)
xlabel('delay','FontSize',12)
ylabel('change in synaptic precicion','FontSize',12)
title('Spike-timing dependent plasticity','FontSize',16)
axis square, box off
