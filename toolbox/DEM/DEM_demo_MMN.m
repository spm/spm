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
% $Id: DEM_demo_MMN.m 4146 2010-12-23 21:01:39Z karl $
 

% level 1
%--------------------------------------------------------------------------
M(1).m  = 1;                                 % 1 input or cause
M(1).n  = 2;                                 % 2 hidden states
M(1).l  = 2;                                 % 3 outputs (coefficients for face)
 
% first stimulus parameters
%--------------------------------------------------------------------------
A.g     = [ 0  1 ;                           % amplitude
            4  0];                           % frequency
A.f     = [-1  4;
           -2 -1]/16;                        % The Jacobian
ip      = [2];
pC      = sparse(ip,ip,1,8,8);               % covariance
 
% second stimulus parameters
%--------------------------------------------------------------------------
B.g     = [ 0  1;
            1  0];                           % The mixing parameters
B.f     = A.f;
 
        
M(1).f  = inline('P.f*x + [v; 0]','x','v','P');
M(1).g  = inline('P.g*x + [0; 16]','x','v','P');
M(1).pE = A;                                 % The prior expectation
M(1).pC = pC;                                % The prior covariance
M(1).Q  = {eye(2)};                          % error precision (data)
M(1).hE = 4;                                 % error log-precision prior
M(1).hC = 1/32;                              % error log-precision prior
M(1).W  = exp(16);                           % error log-precision prior
 
% level 2
%--------------------------------------------------------------------------
M(2).l  = 1;                                 % 1 output 
M(2).V  = exp(16);                           % error precision (cause)
 
% The input here is simply a bump [Gaussian] function.
%==========================================================================
n       = 5;                                 % number of trials
N       = 64;                                % length of data sequence
dt      = 0.005;                             % time bin (sec)
T       = N*dt;                              % duration of chirp (sec)
t       = [1:N]*dt;                          % time bin (sec)
c       = exp(-([1:N] - 20).^2/(4.^2));      % this is the Gaussian cause
 
 
% DEM estimation:  Here we expose the model M to the data and record the
% responses.  The DEM scheme is essential a form of Variational Learning
% that provides an upper bound on perceptual inference and learning.
% We use this bound to simulate neuronal responses, under the assumption
% they are near-optimal.
%==========================================================================
M(1).E.s  = 1;                    % temporal smoothness
 
M(1).E.nD = 1;                    % D-steps per time bin (1 for dynamic systems)
M(1).E.nE = 1;                    % E-steps per iteration
M(1).E.nM = 8;                    % M-steps per iteration
 
 
% Because we want to record the response of the model over time, to each
% stimulus we will proceed one iteration at a time and replace the starting
% values of the parameters (initialised with the prior expectation) with
% the conditional estimates of the previous trial.  We will consider 8
% presentations
%--------------------------------------------------------------------------
DEM   = {};
for i = 1:(n + 2);
    
    % Change stimulus
    %----------------------------------------------------------------------
    if i <= 2
        DEM{i} = spm_DEM_generate(M,c,A,{16 16},{16 []});
    else
        DEM{i} = spm_DEM_generate(M,c,B,{16 16},{16 []});
    end
    DEM{i}.U = c;
    
    % Invert
    %---------------------------------------------------------------------- 
    DEM{i}  = spm_DEM(DEM{i});           % compute conditional densities
    M(1).pE = DEM{i}.qP.P{1};            % update parameter estimates
    M(1).hE = DEM{i}.qH.h{1};            % update hyperparameter estimates
end
DEM   = DEM(2:end);                      % discard burn-in trial
 
 
for i = 1:(n + 1)
    
    % gather trial specific parameters and precisions changes
    %----------------------------------------------------------------------
    spm_figure('GetWin','Figure 1');
    subplot(n + 1,3,(i - 1)*3 + 3)
    
    dP{i} = spm_vec(DEM{i}.M(1).pE) - spm_vec(DEM{end}.M(1).pE);
    qR{i} = spm_DEM_MEG(DEM{i},dt,1,1);       % prediction error (LFP)
    qR{i} = spm_DEM_EEG(DEM{i},dt,[1 2],1);   % prediction error (LFP)
    qH(i) = DEM{i}.M(1).hE;                   % and precision
    drawnow
    if i == 2, a = axis; end
 
    % plot recognition density and prediction
    %----------------------------------------------------------------------
    subplot(n + 1,3,(i - 1)*3 + 1)
    spm_plot_ci([1:N]*dt,DEM{i}.qU.x{1},DEM{i}.qU.S)
    hold on
    plot([1:N]*dt,DEM{i}.pU.x{1},':')
    hold off
    
    subplot(n + 1,3,(i - 1)*3 + 2)
    spm_DEM_play_song(DEM{i}.qU,T);
    axis square
    
end

% scale axes and title
%--------------------------------------------------------------------------
for i = 1:(n + 1)
    subplot(n + 1,3,(i - 1)*3 + 3)
    axis(a);
    axis square
end
subplot(n + 1,3,1),title('Hidden states','FontSize',16)
subplot(n + 1,3,2),title('percept','FontSize',16)
subplot(n + 1,3,3),title('prediction error','FontSize',16)
    
% Show song in DEM window
%--------------------------------------------------------------------------
spm_figure('GetWin','DEM');
colormap('pink')
subplot(2,2,4)
spm_DEM_play_song(DEM{n + 1}.qU,T);
axis square
 
% Repetition effects
%==========================================================================
spm_figure('GetWin','Figure 2'); 
t   = [1:N]*dt*1000;
 
% Changes in Parameters
%--------------------------------------------------------------------------
subplot(3,2,1)
qP  = spm_cat(dP)';
bar(abs(qP(2:end,ip)))
title( 'extrinsic connectivity','FontSize',16)
xlabel('presentation','FontSize',12)
ylabel('changes in parameters','FontSize',12)
axis square
 
% Changes in log-precision
%--------------------------------------------------------------------------
subplot(3,2,2)
bar(qH(2:end))
title( 'intrinsic connectivity','FontSize',16)
xlabel('presentation','FontSize',12)
ylabel('hyperparameters','FontSize',12)
axis square
 
% ERPs (precision weighted) prediction error
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2');
a     = [t(1) t(end) -30 30];
b     = [t(1) t(end) -1 1];
for i = 1:n
    
    subplot(6,n,i + 2*n)
    plot(t,qR{i + 1}{1},'r')
    axis square tight
    axis(a)
    box off
        
    subplot(6,n,i + 3*n)
    plot(t,qR{i + 1}{2},'r')
    axis square tight
    axis(b)
    box off
    
end
 
% Differences
%--------------------------------------------------------------------------
subplot(3,2,5)
plot(t,qR{2}{1} - qR{n + 1}{1},'r')
title('primary level (N1/P1)','FontSize',16)
xlabel('peristimulus time (ms)')
ylabel('Difference waveform','FontSize',12)
axis square
 
subplot(3,2,6)
plot(t,qR{2}{2} - qR{n + 1}{2},'r')
title('secondary level (MMN)','FontSize',16)
xlabel('peristimulus time (ms)','FontSize',12)
ylabel('Difference waveform','FontSize',12)
axis square


drawnow, disp(' '),disp('Click sonograms to play songs'),disp(' ')
