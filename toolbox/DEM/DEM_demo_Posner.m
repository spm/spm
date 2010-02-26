% This demonstration routine simulates the Posner paradigm to show that 
% some of the characteristic speed-accuracy trade-offs associated with 
% valid and invalid cueing can be explained easily in terms of optimizing 
% precisions during hierarchical inference. This demonstration uses 
% generalised filtering and state-space model that includes state-dependent
% noise. Here, this dependency is used to set the attentional gain or bias 
% using a cue, which modulates the prediction errors induced by subsequent 
% targets. The phenomena that emerge from this scheme include a competition
% for attentional resources; given that only one context can exist at any 
% time and this probabilistic context is encoded by state-dependent 
% precisions on the causes of sensory input. Attended stimuli have greater 
% precision and greater penetration of their prediction errors in the 
% hierarchy. We will also see characteristic differences between perceptual 
% inference, under valid and invalid cues. This is illustrated using 
% simulated psychophysical and electrophysiological responses. Biased 
% competition is simulated by presenting both valid and invalid targets 
% simultaneously.

%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_Posner.m 3740 2010-02-26 13:13:14Z karl $
clear
 

% Create a generative model: To keep thing simple we will model just 2 
% target locations; on the left (invalid) and the right (valid).
%==========================================================================

% set dimensions for generalised coordinates
%--------------------------------------------------------------------------
G(1).E.d        = 2;                   % approximation order
G(1).E.n        = 4;                   % embedding order
G(1).E.s        = 1/2;                 % temporal smoothness
G(1).E.method.x = 1;                   % state-dependent noise
G(1).E.sP       = 1e6;                 % smoothness precision
N               = 2;                   % number of stimuli (locations)
 
                                       
% level 1; with state-dependent variance
%--------------------------------------------------------------------------
G(1).m  = N + 1;                       % causes (inputs)
G(1).n  = 2;                           % hidden states (context)
G(1).l  = N + 1;                       % output channels (stimuli)
G(1).xP = 1/32;                        % prior precision on hidden state
 
 
G(1).f  = inline('[1 -1 -2;-1 1 2]*v/4 - x/32','x','v','P');
G(1).g  = inline('v','x','v','P');
G(1).V  = exp(16);                     % error variances (noise)
G(1).W  = exp(8);                      % error variances (states)
 
 
% level 2; causes
%--------------------------------------------------------------------------
G(2).l  = N + 1;                       % output channels (locations)
G(2).V  = exp(16);
 
 
% state-dependent precision (attentional bias) in generative model (M):
%--------------------------------------------------------------------------
M       = G;
M(1).ph = inline('[2 + [1 -1;-1 1]*x; 4]','x','v','h','M');
M(1).W  = exp(4);                      % error variances (states)
M(1).V  = [];
M(2).V  = 0;
 
 
% Data (stimuli] are created by integrating the model for some input. The 
% input here comprises cue and target stimuli modeled with bump (Gaussian)
% functions. 
%==========================================================================
 
% create inputs and cue 
%--------------------------------------------------------------------------
T   = 64;                                   % length of data sequence
dt  = 640/T;                                % ms time bins; 512 ms trials
pst = [1:T]*dt;                             % peristimulus time (ms bins)
TS  = spm_Npdf([1:T + 1],2*T/3,T/8);        % this is the Gaussian target
CS  = spm_Npdf([1:T],T/4,T/16);             % this is the Gaussian cue
TS  = diff(TS);
TS  = TS/max(TS);
CS  = CS/max(CS);

% valid target (right)
%--------------------------------------------------------------------------
V(N,    :) = TS;
V(N + 1,:) = CS;
 
% invalid target (left)
%--------------------------------------------------------------------------
U(1,    :) = TS;
U(N + 1,:) = CS;
 
% both targets
%--------------------------------------------------------------------------
B(1,    :) = TS;
B(N,    :) = TS;
B(N + 1,:) = CS;
 

% integrate G to obtain stimuli (DEM.Y)
%--------------------------------------------------------------------------
LAPV     = spm_DEM_generate(G,V);
LAPU     = spm_DEM_generate(G,U);
LAPB     = spm_DEM_generate(G,B);
 

% Filtering: Here. we expose the model M to the data and record the
% responses.  The scheme is essentially a form of Variational Learning
% that provides an upper bound on perceptual inference and learning.
% We use this bound to simulate neuronal responses, under the assumption
% they are near-optimal.
%==========================================================================
 
% place the generative model in LAP
%--------------------------------------------------------------------------
LAPV.M = M;
LAPU.M = M;
LAPB.M = M;
 
% simulate a valid trial
%==========================================================================
spm_figure('GetWin','DEM');

LAPV   = spm_LAP(LAPV);

spm_figure('GetWin','Figure 1');
spm_DEM_qU(LAPV.qU,LAPV.pU)
 
 
% show stimulus and percept
%--------------------------------------------------------------------------
subplot(2,2,4)
load Posner

qU    = LAPV.qU.v{2};
for i = 1:length(qU)
    image(64*spm_unvec(S.V*qU(:,i),S.F))
    axis image off
    mov(i) = getframe;
end
 
% set ButtonDownFcn for movie of perceived (inferred) stimuli
%--------------------------------------------------------------------------
h = get(gca,'Children');
set(h(1),'Userdata',{mov,T})
set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title({'conditional percept','click on image'},'FontSize',16)
 
 
 
% simulate an invalid trial
%==========================================================================
spm_figure('GetWin','DEM');

LAPU   = spm_LAP(LAPU);

spm_figure('GetWin','Figure 2');
spm_DEM_qU(LAPU.qU,LAPU.pU)
 
 
% compare targets with and without valid cues
%--------------------------------------------------------------------------
subplot(2,2,4)
 
spm_plot_ci(pst,LAPV.qU.v{2},LAPV.qU.C,2,'g'), hold on
spm_plot_ci(pst,LAPU.qU.v{2},LAPU.qU.C,1,'b'), hold off
 
set(gca,'XLim',[pst(1) pst(end)])
xlabel('time (ms)','FontSize',12)
title({'target responses with','valid and invalid cues'},'FontSize',16)
axis square
box off
drawnow
 
 
% Psychophysics - speed accuracy trade-off
%==========================================================================
spm_figure('GetWin','Figure 3');

% valid trial - evaluate p(v > q,t|Y,M) using conditional density
%--------------------------------------------------------------------------
q   = 1/4;                                % threshold
c   = [0 1 0];                            % contrast (valid target)
E   = LAPV.qU.v{2};                       % conditional expectations
C   = LAPV.qU.C;                          % conditional covariance
 
% Conditional moments of time-averaged parameters
%--------------------------------------------------------------------------
for i = 1:T
    Ec    = c*E(:,i);                     % Expectation of contrast
    Cc    = c*C{i}*c';                    % Covariance  of contra
    PV(i) = 1 - spm_Ncdf(q,Ec,Cc);        % Posterior p(Ec > q)
end
 
% Repeat for invalid trial
%--------------------------------------------------------------------------
c   = [1 0 0];                            % contrast (invalid target)
E   = LAPU.qU.v{2};                       % conditional expectations
C   = LAPU.qU.C;                          % conditional covariance
for i = 1:T
    Ec    = c*E(:,i);                     % Expectation of contrast
    Cc    = c*C{i}*c';                    % Covariance  of contra
    PU(i) = 1 - spm_Ncdf(q,Ec,Cc);        % Posterior p(Ec > q)
end
 
% Graphics
%--------------------------------------------------------------------------
subplot(2,1,1)
plot(pst,100*PV,pst,100*PU,'-.',pst,TS*100,':')

[m i] = max(PU);
set(gca,'XLim',[pst(T/2) pst(i)])
xlabel('peristimulus time (ms)','FontSize',12)
ylabel('posterior confidence that target is present','FontSize',12)
title({'speed-accuracy trade-off'},'FontSize',16)
axis square
box off
legend('valid','invalid','true intensity')
drawnow
 
% Biased competition
%==========================================================================
spm_figure('GetWin','DEM');

% simulate dual presentation
%--------------------------------------------------------------------------
LAPB   = spm_LAP(LAPB);

spm_figure('GetWin','Figure 4');
spm_DEM_qU(LAPB.qU,LAPB.pU)
 
 
% compare non-attended stimulus with and without attended stimulus
%--------------------------------------------------------------------------
subplot(2,2,4)
 
spm_plot_ci(pst,LAPB.qU.v{2},LAPB.qU.C,1,'b-.'), hold on
spm_plot_ci(pst,LAPU.qU.v{2},LAPU.qU.C,1,'b'), hold off
 
set(gca,'XLim',[pst(1) pst(end)])
xlabel('time (ms)','FontSize',12)
title({'unattended stimulus','with and without attended'},'FontSize',16)
axis square
box off
drawnow
 
% Electrophysiology
%==========================================================================
 
% plot the simulated ERPs (see spm_DEM_ERP)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 5'); clf
color = {'r','b'};
qU    = {LAPU.qU,LAPV.qU};
PST   = pst - T*dt/3;
for i = 1:length(qU)
    
    
    % loop over levels
    %----------------------------------------------------------------------
    for  j = 1:2
        
        % PST (assuming 32 ms times bins)
        %------------------------------------------------------------------
        try
            EEG = qU{i}.Z{j}(1:2,:);
        catch
            EEG = qU{i}.z{j}(1:2,:);
        end

        % ERPs
        %------------------------------------------------------------------
        subplot(2,2,1)
        plot(PST,EEG,'Color',color{i}),hold on
    end
    title('LFPs (Causal) Content','FontSize',16)
    xlabel('pst (ms)')
    axis square
    set(gca,'XLim',[PST(1) PST(end)])
    
    % Hidden states
    %------------------------------------------------------------------
    j  = 1;
    try
        EEG = qU{i}.W{j}(2,:);
    catch
        EEG = qU{i}.w{j}(2,:);
    end
    
    % ERPs
    %------------------------------------------------------------------
    subplot(2,2,2)
    plot(PST,EEG,'Color',color{i}),hold on
    title('LFPs (Hidden) Context','FontSize',16)
    xlabel('pst (ms)')
    axis square 
    set(gca,'XLim',[PST(1) PST(end)])
    
end

legend('invalid','valid')
drawnow
