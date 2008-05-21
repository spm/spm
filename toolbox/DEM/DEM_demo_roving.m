% This Demo uses the linear convolution model of previous examples to
% simulate chirps (c.f., the bird song demos). By presenting a train of
% chirps and changing the stimulus after a couple of presentations, we can
% simulate a roving oddball paradigm used in ERP research. Critically, we
% hope to see a more exuberant response to the first presentation of a
% novel chirp (oddball) relative to the same stimulus after learning
% (standard).  The simulation shows that although veridical percepts obtain
% from variational de-convolution, the prediction error continues to fall
% with repetition (as the parameters are optimised). This repetition
% suppression subtends a mismatch response that has many of the
% characteristics of the well-known mismatch negativity (MMN).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_demo_roving.m 1703 2008-05-21 13:59:23Z karl $
 
% figure
%--------------------------------------------------------------------------
spm_figure('GetWin','Graphics'); 
colormap('pink')
 
% level 1
%--------------------------------------------------------------------------
M(1).m   = 1;                               % 1 input or cause
M(1).n   = 2;                               % 2 hidden states
M(1).l   = 2;                               % 3 outputs (coefficients for face)
 
% first stimulus parameters
%--------------------------------------------------------------------------
P.f     = [-1  4 ;
           -2 -1]/16;                        % The Jacobian
P.g     = [ 0  1 ;
            4  1]*2;                        % The mixing parameters
 
% second stimulus parameters
%--------------------------------------------------------------------------
B.f     = [-1  2 ;
           -2 -1]/16;                        % The Jacobian
B.g     = [ 1  0 ;
            3 -1]*2;                        % The mixing parameters
 
M(1).f  = inline('P.f*x + [v; 0]','x','v','P');
M(1).g  = inline('P.g*x + [0; 16]','x','v','P');
M(1).pE = P;                                 % The prior expectation
M(1).pC = speye(length(spm_vec(P)))*exp(8);  % The prior covariance
M(1).V  = exp(2);                            % error precision (data)
M(1).W  = exp(8);                            % error precision (states)                                          % with a low level of noise
 
% level 2
%--------------------------------------------------------------------------
M(2).l  = 1;                                 % 1 output 
M(2).V  = exp(2);                            % error precision (cause)
 
% The input here is simply a bump [Gaussian] function.
%==========================================================================
N       = 64;                                % length of data sequence
T       = 1/2;                               % duration of chirp (sec)
dt      = T/N;                               % time bin (sec)
c       = exp(-([1:N] - 16).^2/(2.^2));      % this is the Gaussian cause


 
% DEM estimation:  Here we expose the model M to the data and record the
% responses.  The DEM scheme is essential a form of Variational Learning
% that provides an upper bound on perceptual inference and learning.
% We use this bound to simulate neuronal responses, under the assumption
% they are near-optimal.
%==========================================================================
M(1).E.s  = 1;                    % temporal smoothness
 
M(1).E.nD = 1;                    % D-steps per time bin (1 for dynamic systems)
M(1).E.nE = 1;                    % E-steps per iteration
M(1).E.nM = 1;                    % M-steps per iteration
M(1).E.nN = 1;                    % number of iterations
 
 
% Because we want to record the response of the model over time, to each
% stimulus we will proceed one iteration at a time and replace the starting
% values of the parameters (initialised with the prior expectation) with
% the conditional estimates of the previous trial.  We will consider 8
% presentations
%--------------------------------------------------------------------------
Na    = 2;                        % number of chirp A
Nb    = 4;                        % number of chirp B
Nc    = Na + Nb;                  % number of chirps
qR    = {};
qE    = [];
for i = 1:Nc
    
    % Change stimulus
    %----------------------------------------------------------------------
    if i <= Na
        DEM = spm_DEM_generate(M,c,P,{4 16},{16 []});
    else
        DEM = spm_DEM_generate(M,c,B,{4 16},{16 []});
    end
 
    
    % Invert
    %---------------------------------------------------------------------- 
    DEM     = spm_DEM(DEM);           % compute conditional densities
    M(1).pE = DEM.qP.P{1};            % update parameter estimates
    
    % plot
    %----------------------------------------------------------------------
    spm_figure('GetWin','Graphics'); 
    subplot(Nc,3,(i - 1)*3 + 1)
    spm_plot_ci([1:N]*dt,DEM.qU.x{1},DEM.qU.S)
    hold on
    plot([1:N]*dt,DEM.pU.x{1},':')
    axis tight
    set(gca,'YLim',[-4 4])
    
    subplot(Nc,3,(i - 1)*3 + 2)
    spm_DEM_play_song(DEM.qU,T);
    axis square
    
    subplot(Nc,3,(i - 1)*3 + 3)
    qR{i} = spm_DEM_EEG(DEM,dt);
    qE(i) = sum(spm_vec(qR{i}).^2);   % record the prediction error
    axis tight
    set(gca,'YLim',[-8 12])
    drawnow
 
end
 
subplot(Nc,3,1),title('Hidden states','FontSize',16)
subplot(Nc,3,2),title('percept','FontSize',16)
subplot(Nc,3,3),title('prediction error','FontSize',16)

% Show song in DEM window
%--------------------------------------------------------------------------
spm_figure('GetWin','DEM');
colormap('pink')
subplot(2,2,4)
spm_DEM_play_song(DEM.qU,T);
axis square
 
 
% ERPs to [first] oddball) and standard [last Oddball]
%--------------------------------------------------------------------------
spm_figure('GetWin','MFM'); 
t   = [1:N]*dt*1000;
 
subplot(3,1,1)
bar(qE)
title('SSQ prediction error','FontSize',16)
xlabel('repetition')
axis square
 
subplot(3,2,4)
plot(t,qR{Na + 1}{1},':r',t,qR{Na + 1}{2},'r')
title('LFP: Oddball','FontSize',16)
xlabel('peristimulus time (ms)')
axis square
set(gca,'YLim',[-8 8])
 
subplot(3,2,3)
plot(t,qR{Na + Nb}{1},':r',t,qR{Na + Nb}{2},'r')
title('LFP: Standard (P1)','FontSize',16)
xlabel('peristimulus time (ms)')
axis square
set(gca,'YLim',[-8 8])
 
subplot(3,1,3)
plot(t,qR{Na + 1}{1} - qR{Na + Nb}{1},':r',t,qR{Na + 1}{2} - qR{Na + Nb}{2},'r')
title('Difference waveform (MMN)','FontSize',16)
xlabel('peristimulus time (ms)')
axis square
set(gca,'YLim',[-8 8])

legend({'primary area - amplitude','primary area - frequency','secondary area - chirp'})
