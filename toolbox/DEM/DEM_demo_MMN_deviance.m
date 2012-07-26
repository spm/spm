function DEM_demo_MMN_deviance
% This Demo is a more refined illustration of the mismatch negativity that
% uses a plausible (nonlinear) mapping from hidden states to sensory
% samples; namely, a discrete spectral density (with a fixed frequency or
% pitch tuning). This is modelled using a radial basis function over
% frequencies, whose location parameter is a hidden state encoding pitch and
% whose amplitude is modulated dynamically by hidden states encoding
% amplitude modulation (the sum of two squared hidden states showing damped
% linear oscillation when perturbed by a hidden cause). The recognition
% dynamics illustrate a dissociation between the effects of changing the
% (i) the level of deviancy of a stimulus (modelled here as a deviation of 
% the hidden state (cause) encoding pitch from zero) and (ii) the
% probability of encountering a pitch deviant. This is modelled by changing
% the precision on the hidden (pitch) state. Crucially, the nonlinearities
% in this more plausible generative model of pure tone stimuli induce a
% latency difference in the mismatch response and increase the amplitude
% of the mismatch (i.e., MMN). Conversely, changing the precision only
% affects the amplitude. In these simulations the MMN is modelled simply as
% the difference between  prediction errors evoked by an expected stimulus
% (the standard) and a deviant stimulus. Prior expectations are encoded
% in terms of hidden causes, where the onset of the stimulus is known.
% This means the agent has correct prior expectations about amplitude
% modulation (i.e., knows when to expect a stimulus) but can have
% incorrect expectations about its pitch (of varying confidence or
% precision).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_MMN_deviance.m 4804 2012-07-26 13:14:18Z karl $
 
 
% Create the generative model
%==========================================================================
 
% level 1
%--------------------------------------------------------------------------
M.E.s   = 1;                                 % temporal smoothness
P       = [-1  4;                            % parameters of motion or
           -2 -1]/16;                        % Jacobian
 
M(1).f  = inline('[P*x(1:2); (v(2) - x(3))/32] + [v(1); 0; 0]','x','v','P');
M(1).g  = inline('sum(x(1:2).^2)*exp(-([-2;-1;0;1;2] - x(3)*2).^2/2)/32','x','v','P');
M(1).pE = P;                                 % The prior expectation
M(1).x  = zeros(3,1);                        % hidden states
M(1).V  = exp(4);                            % precision (data)
M(1).W  = exp(4);                            % precision (motion)
M(1).xP = exp(-4);                           % precision on states
 
% level 2
%--------------------------------------------------------------------------
M(2).v  = zeros(2,1);                        % hidden causes
M(2).V  = diag(exp([2,0]));                  % precision
 
 
% Create hidden causes that determine amplitude and pitch
%--------------------------------------------------------------------------
N       = 64;                                % length of data sequence
dt      = 0.005;                             % time bin (sec)
T       = N*dt;                              % duration of chirp (sec)
pst     = [1:N]*dt*1000;                     % time bin (ms)
C       = exp(-([1:N] - 20).^2/(8.^2));      % amplitude
U       = [C; (zeros(1,N))];                 % pitch
 

% Difference waveforms under different levels of deviance (from zero)
%==========================================================================
DEM   = {};
D     = [0 0.75 1 1.25 1.5];
n     = length(D);
for i = 1:n
    
    % generate stimulus
    %----------------------------------------------------------------------
    C(1,:)   = U(1,:);
    C(2,:)   = U(2,:) + D(i);
    DEM{i}   = spm_DEM_generate(M,C,P,{16 16},{16});
    
    % solve (under prior expectations about a standard stimulus)
    %----------------------------------------------------------------------
    DEM{i}.U = U;
    DEM{i}   = spm_DEM(DEM{i});
    
end
 
% graphics (ERPs at different levels)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');

for i = 1:n
    
    subplot(n,2,(i - 1)*2 + 1)
    qx{i} = spm_DEM_MEG(DEM{i},dt,1,1);
    str = sprintf('(N1/P1) deviance: %.2f',D(i));
    title(str,'FontSize',12)
    
    subplot(n,2,(i - 1)*2 + 2)
    qv{i} = spm_DEM_MEG(DEM{i},dt,2,1);
    title('secondary level (MMN)','FontSize',12)
    drawnow
end
 

% graphics (difference waveforms)
%--------------------------------------------------------------------------
for i = 1:n
    MMN(i,:) = sum(qx{i}{1}) - sum(qx{1}{1});
    leg{i}   = sprintf('deviance: %.2f',D(i));
end

spm_figure('GetWin','Figure 2');

subplot(2,1,1), plot(pst,MMN), legend(leg), drawnow
title('Difference waveforms under different levels of deviance','FontSize',16)
xlabel('prestimulus time (ms)','FontSize',12)
 
 
% Difference waveforms under different levels of confidence (D)
%==========================================================================
DEM   = {};
D     = [0:4] - 4;
n     = length(D);
for i = 1:n
    
    % change probabilistic context - log-precision of hidden pitch
    %----------------------------------------------------------------------
    M(1).xP   = diag(exp([-8 -8 D(i)]));
    
    % generate standard
    %----------------------------------------------------------------------
    C(1,:)   = U(1,:);
    C(2,:)   = U(2,:);
    DEM{i}   = spm_DEM_generate(M,C,P,{16 16},{16});
    
    % solve and record sensory prediction errors
    %----------------------------------------------------------------------
    DEM{i}.U = U;
    DEM{i}   = spm_DEM(DEM{i});
    qs{i}    = spm_DEM_MEG(DEM{i},dt,1);
    
    % generate deviant
    %----------------------------------------------------------------------
    C(1,:)   = U(1,:);
    C(2,:)   = U(2,:) + 2;
    DEM{i}   = spm_DEM_generate(M,C,P,{16 16},{16});
    
    
    % solve and record sensory prediction errors
    %----------------------------------------------------------------------
    DEM{i}.U = U;
    DEM{i}   = spm_DEM(DEM{i});
    qd{i}    = spm_DEM_MEG(DEM{i},dt,1);
    
end
 
 
% graphics (differences)
%--------------------------------------------------------------------------
for i = 1:n
    MMN(i,:) = sum(qd{i}{1}) - sum(qs{1}{1});
    leg{i}   = sprintf('log-precision: %.2f',D(i));
end

spm_figure('GetWin','Figure 2');

subplot(2,1,2), plot(pst,MMN), legend(leg)
title('Difference waveforms under different levels of confidence','FontSize',16)
xlabel('prestimulus time (ms)','FontSize',12)
