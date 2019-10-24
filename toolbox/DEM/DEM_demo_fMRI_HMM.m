function DEM_demo_fMRI_HMM
% Demonstration of Hidden Markov models for fMRI
%__________________________________________________________________________
%  This demonstration routine illustrates the modelling of state
%  transitions generating resting state fMRI timeseries. The hidden states
%  are modelled as a hidden Markov model, where each state corresponds to a
%  particular  point in the parameter space of effective connectivity. This
%  effective connectivity then generates complex cross spectral data
%  features  of the observed timeseries. Model specification requires prior
%  constraints on the probability transition matrix among hidden states,
%  which implicitly specifies the number of hidden states. The user also
%  has to specify the number of windows for epochs to apply to the
%  timeseries, where each epoch  places a lower bound on the duration of
%  each (discrete) state.
%     We first generate synthetic data using regular transitions among
%  three hidden states  (C.F., a discrete version of a heteroclinic
%  cycle  for orbit). The data are then converted by a routine that
%  combines a parametric empirical Bayesian model and a hidden Markov model
%  (as implemented as a special case of a Markov decision process). This
%  inversion is repeated for each model specified in terms of the
%  transition matrices (as prior Dirichlet concentration parameters).
%  Setting a prior transition parameter to 0 precludes that transition. In
%  this way, several different models of transitions and number of hidden
%  states can be  scored in terms of  the variational free energy.
%     Following inversion, the results are plotted in terms of expected
%  state transitions, fluctuations in connections that are allowed to
%  change (specified in the usual way by DCM.b), the deviations in
%  connectivity associated with each hidden state and the expected
%  probability transition matrix.
%     Finally, we consider Bayesian model comparison in terms of group
%  differences (here, simply the difference between the first and second
%  simulated subject).  Bayesian model comparison is simple to do in this
%  context  by comparing the free energy of a hidden Markov model in which
%  both groups share the same state dependent connections and transition
%  probabilities, with two independent models. These can be evaluated
%  efficiently using Bayesian model reduction implicit in PEB. in this
%  example, we did not introduce any differences between the two groups
%  (i.e., subjects) and therefore expected to infer no group effect.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_demo_fMRI_HMM.m 7679 2019-10-24 15:54:07Z spm $



% (I) Simulate fMRI timeseries
%==========================================================================
rng('default')

%  Assume we have P sessions with N epochs of T scans with a TR of 2 secs:
% -------------------------------------------------------------------------
%  These epochs could be from a single subject or result from the
%  concatenation of multiple sessions (under the assumption that they share
%  the same modes of connectivity).
% -------------------------------------------------------------------------
P  = 2;                               % number of sessions (e.g., subjects)
S  = 3;                               % number of latent (hidden) states
N  = 9;                               % number of epochs or windows
T  = 128;                             % number of observations (per epoch)
TR = 2;                               % repetition time or timing
t  = (1:(T*N))*TR;                    % observation times (seconds)
n  = 3;                               % number of regions or nodes

% setup model for generating timeseries
% -------------------------------------------------------------------------
options.nonlinear  = 0;
options.two_state  = 0;
options.stochastic = 0;
options.induced    = 1;

% get priors to generate simulated data
% -------------------------------------------------------------------------
a   = ones(n,n);
b   = zeros(n,n,0);
c   = zeros(n,n);
d   = zeros(n,n,0);
pP  = spm_dcm_fmri_priors(a,b,c,d,options);

% average parameters - a simple hierarchy of three nodes
% -------------------------------------------------------------------------
pP.A = [  0  -.3    0;
         .4    0  -.1;
          0   .3    0];
pP.C = eye(n,n);
pP.transit = randn(n,1)/64;

% generate spectral density of neuronal fluctuations and observation noise
% -------------------------------------------------------------------------
[Gu,Gn,Hz,dt] = spm_csd_fmri_gu(pP,TR);
Gu    = Gu(:,1,1)*ones(1,n);
Gn    = Gn(:,1,1)*ones(1,n);


% specify and generate Markovian succession of hidden states
%==========================================================================

% connections associated with hidden states: here, intrinsic connectivity
% -------------------------------------------------------------------------
B        = zeros(n,n,S);
B(1,1,1) = 1/2;
B(2,2,2) = 1/2;
B(3,3,3) = 1/2;
pP.B     = any(B,3);                     % state-dependent connections

% generate sequence of hidden states: here, a simple orbit
% -------------------------------------------------------------------------
bb    = spm_speye(S,S,-1); bb(1,S) = 1;
o     = kron(ones(1,N),1:S);
o     = o(1:N);

% State-dependent deviations in connectivity from average
% -------------------------------------------------------------------------
for i = 1:N
    tB(:,:,i) = B(:,:,o(i));
end
for i = 1:S
    for j = 1:S
        tB(i,j,:) = tB(i,j,:) - mean(tB(i,j,:));
    end
end

%% simulate epoch-specific responses to endogenous fluctuations
%==========================================================================
M.x   = sparse(n,5);
M.f   = 'spm_fx_fmri';
X     = cell(N,1);
Y     = cell(N,1);
E     = cell(N,1);
u     = cell(N,1);
for p = 1:P
    
    % parameters for this epoch, plus a small random effect
    % ---------------------------------------------------------------------
    gu    = spm_rand_power_law(Gu,Hz,dt,N*T);
    ge    = spm_rand_power_law(Gn,Hz,dt,N*T);
    for s = 1:N
        
        % parameters for this epoch, plus a small random effect
        % -----------------------------------------------------------------
        tP   = pP;
        tP.A = tP.A + tB(:,:,s);
        tP.A = tP.A.*(1 + randn(n,n)/64);
        tP.C = eye(n,n);
        
        % integrate states with endogenous fluctuations (gu)
        % -----------------------------------------------------------------
        j    = (1:T) + (s - 1)*T;
        M.f  = 'spm_fx_fmri';
        U.u  = gu(j,:);
        U.dt = TR;
        x    = spm_int_J(tP,M,U);
        M.x  = spm_unvec(x(end,:),M.x);
        
        % haemodynamic observer function to produce BOLD signal
        % -----------------------------------------------------------------
        for i = 1:T
            y(i,:) = spm_gx_fmri(spm_unvec(x(i,:),M.x),[],tP)';
        end
        
        % response with observation noise (ge)
        % -----------------------------------------------------------------
        e       = ge(j,:);
        X{s}    = x;
        Y{s}    = y + e;
        E{s}    = e;
        u{s}    = U.u;
        TP(s,p) = tP;
        
    end
    
    
    % concatenate epochs into a single timeseries
    %----------------------------------------------------------------------
    xY.dt = TR;
    xY.y  = spm_cat(Y);
    xY.u  = spm_cat(u);
    xY.X  = spm_cat(X);
    xY.E  = spm_cat(E);
    
    % and create DCM cell array
    %----------------------------------------------------------------------
    DCM{1,p}.options = options;
    DCM{1,p}.a = logical(pP.A);
    DCM{1,p}.b = logical(pP.B);
    DCM{1,p}.c = zeros(n,0);
    DCM{1,p}.d = zeros(n,n,0);
    DCM{1,p}.Y = xY;
    
end

%% show simulated responses and windows
%--------------------------------------------------------------------------
spm_figure('Getwin','Figure 1'); clf

subplot(3,2,1), plot(t,xY.u)
title('Endogenous fluctuations','FontSize',16)
xlabel('Time (seconds)'), ylabel('Amplitude'), axis square, spm_axis tight

subplot(3,2,2), hold off
plot(t,xY.X(:,(n + 1):end),'c'), hold on
plot(t,xY.X(:,1:n)),             hold off
title('Hidden states','FontSize',16)
xlabel('Time (seconds)'), ylabel('Amplitude'), axis square, spm_axis tight

subplot(3,2,3)
plot(t,xY.y,t,xY.E,':')
title('Hemodynamic response and noise','FontSize',16)
xlabel('Time (seconds)'), ylabel('Amplitude'), axis square, spm_axis tight


%  This completes the simulation of the data. We now turn to inverting the
%  data to see if one can recover the number of hidden states, the form of 
%  the state transitions and the connectivity modes associated with each
%  state:
%--------------------------------------------------------------------------


% (II) Inversion under a hidden Markov model
%==========================================================================
%  Specify model space as a cell array of probability transition matrices:
%  here, the model space at the level of the HMM  will allow all
%  transitions among one to 4 hidden states. These models are specified in
%  terms of Dirichlet priors; starting with a small value of allowable
%  transitions (1/16)
%--------------------------------------------------------------------------
for i = 1:4
    b{i} = ones(i,i)/16;
end

% invert hidden Markov model: this is the routine demonstrated
%--------------------------------------------------------------------------
[HMM,CSD] = spm_dcm_HMM(DCM,N,b);

% This completes the inversion. We now just need to look at the results:
%--------------------------------------------------------------------------



%% (III) report analysis
%==========================================================================
spm_figure('Getwin','Figure 1');

%  plot windows
% -------------------------------------------------------------------------
subplot(3,2,3), hold on
for i = 1:N, plot(t,CSD{i,end}.W - 1), end, hold off

% show estimates for a single session
% -------------------------------------------------------------------------
subplot(3,2,4)
spm_plot_ci(CSD{end}.Ep,CSD{end}.Cp), hold on
bar(TP(end).A(:),1/4), hold off, axis square
title('True and MAP connections (Deterministic)','FontSize',16)

% show state-dependent changes in connectivity over sessions
% -------------------------------------------------------------------------
for i = 1:numel(CSD)
    tp(i,:) = spm_vec(TP(i).A);
    qp(i,:) = spm_vec(CSD{i}.Ep.A);
    pp(i,:) = spm_vec(HMM(S).Ep{i}.A);
end
subplot(3,3,7); imagesc(tp)
title('True connections','FontSize',16), axis square
subplot(3,3,8); imagesc(qp)
title('MAP estimates',   'FontSize',16), axis square
subplot(3,3,9); imagesc(pp)
title('PEB estimates',   'FontSize',16), axis square

% report hidden Markov model
%==========================================================================
spm_dcm_HMM_plot(HMM,S)

% And overlay true values, as cyan dots
%==========================================================================

% true state transitions
%--------------------------------------------------------------------------
x     = sparse(o,1:N,1,S,N);
x     = kron(ones(1,P),x);
N     = size(x,2);

% associate true and discovered states - and reorder
%--------------------------------------------------------------------------
r     = x*HMM(S).X';
j     = zeros(S,1);
for i = 1:S
    [d,m]  = max(r(:,i));
    j(i)   = m;
    r(m,:) = 0;
end
[o,i] = find(x(j,:));
B     = B(:,:,j);
bb    = bb(j,j);

% superimpose true values
%--------------------------------------------------------------------------
spm_figure('Getwin','HMM')

%  hidden states
%--------------------------------------------------------------------------
subplot(4,1,1), hold on
for i = 1:N, plot(i,o(i),'.c','MarkerSize',32), end, hold off

% state-dependent parameters - fluctuations
%--------------------------------------------------------------------------
subplot(4,1,2), hold on
for i = 1:N
    [j,k] = max(spm_vec(TP(i).A - pP.A));
    plot(i,k,'.c','MarkerSize',32)
end, hold off

subplot(4,1,3), hold on
for i = 1:N, pA(:,i) = spm_vec(TP(i).A); end
plot(1:N,pA(HMM(S).iP,:),'-.'), hold off, spm_axis tight

% state-dependent parameters - expectations
%--------------------------------------------------------------------------
subplot(4,2,7), hold on
for i = 1:S
    c     = spm_vec(B(:,:,i));
    [j,k] = max(c(HMM(S).iP));
    plot(i,k,'.c','MarkerSize',32)
end, hold off

% expected transition probabilities
%--------------------------------------------------------------------------
subplot(4,2,8), hold on
for i = 1:S, [j,k] = max(bb(:,i)); plot(i,k,'.c','MarkerSize',32), end
hold off

%% Bayesian model comparison in terms of group (i.e., subject) differences
%==========================================================================

% model as a single group or two separate groups
%--------------------------------------------------------------------------
hmm0 = HMM(S);

hmm1 = spm_dcm_HMM(CSD(:,1),b(S));
hmm2 = spm_dcm_HMM(CSD(:,2),b(S));

% compare the free energy of the combined groups with the combined
% free energy:
%--------------------------------------------------------------------------
F    = [hmm0.F; hmm1.F + hmm2.F];
F    = F - min(F);

% report model comparison in terms of free energy (i.e., log evidence)
%--------------------------------------------------------------------------
spm_figure('Getwin','HMM-F')
subplot(2,2,3)
bar(F,'c'),  title('Group difference','FontSize',16)
xlabel('Effect'), ylabel('Log evidence'), axis square
set(gca,'XTickLabel',{'None','Effect'})

% and show the independent maximum a posteriori estimates of state
% dependent connectivity
%--------------------------------------------------------------------------
subplot(4,2,6)
bar(hmm1.qP),  title('Group 1','FontSize',16)
xlabel('Parameter'), ylabel('Connectivity (log)'), axis square
subplot(4,2,8)
bar(hmm2.qP),  title('Group 2','FontSize',16)
xlabel('Parameter'), ylabel('Connectivity (log)'), axis square

return





