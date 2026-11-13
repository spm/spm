function DEMO_tvec_erp_sim
% DEMO_tvec_erp_sim - Demo script for modelling time-varying effective 
% connectivity in a DCM for ERP.
%
% This script demonstrates key aspects of the modelling approach
% described in the paper. Specifically, it simulates and recovers
% a dynamic causal model (DCM) with time-varying connectivity,
% showcasing the use of temporal basis functions to model slow
% fluctuations in synaptic efficacy.
%
% Overview of this script:
% 1. **Model Setup:**
%    - Defines a simple two-region DCM with forward and backward connections.
%    - Uses a cosine basis set to represent time-varying connectivity modulations.
% 
% 2. **Data Simulation:**
%    - Generates synthetic neural responses based on predefined modulations.
%    - Adds noise to simulate observed data.
% 
% 3. **Parameter Recovery:**
%    - Implements Bayesian model inversion to estimate connectivity changes.
%
% 4. **Visualization:**
%    - Plots the true and recovered connectivity modulations.
%    - Compares simulated, observed, and recovered neuronal responses.
%
%
% Outputs:
% - Visualization of true vs. recovered connectivity modulations.
% - Signal-to-noise ratio (SNR) of synthetic data.
% - Simulated, observed, and recovered neuronal responses.
%
% For further details, refer to the paper:
%   Medrano, J., Friston, K. J., & Zeidman, P. (2024). 
%   Dynamic Causal Models of Time-Varying Connectivity.
%
%__________________________________________________________________________
 
% Johan Medrano
% Copyright (C) 2024-2024 Wellcome Centre for Human Neuroimaging


%% Initialization
%==========================================================================

% Define colormap for plots
%--------------------------------------------------------------------------
c = {'4BE1EC', 'CB5EEE'};
colorgrad = @(N) interp1([1; N], [html2rgb(c{1}); html2rgb(c{2})], 1:N);


%% Region and connectivity configuration
%==========================================================================

% Setup regions
%--------------------------------------------------------------------------
Sname = {'Region 1', 'Region 2'}; % Region names
nr = numel(Sname);               % Number of regions
nsim = 256;                      % Number of simulations

% Base connectivity
%--------------------------------------------------------------------------
A = [];
A{1} = zeros(nr, nr);    % Intrinsic connectivity
A{1}(2, 1) = 1;          % Forward connection (Region 1 -> Region 2)
A{2} = 0.5 * A{1}';      % Backward connection (Region 2 -> Region 1)
A{3} = 0 .* A{1};        % No lateral connections

% Input matrix (driving input into Region 1)
%--------------------------------------------------------------------------
C = zeros(2, 1);
C(1, 1) = 1;

%% Basis set for condition-dependent modulations
%==========================================================================

% Get basis set
%--------------------------------------------------------------------------
nconditions = 16;        % Number of conditions
bf = [];
bf.order = 6;           % Number of basis functions
bf.name = "Cosine set";  % Basis set type
bf.dt = 0.5;             % Time step (seconds)
bf.length = nconditions * bf.dt - bf.dt / 2; % Duration
bf = spm_get_bf(bf);     % Generate basis set
bf.t = (0:bf.dt:bf.length)'; % Time vector

% Process basis set (remove intercept, mean-center)
%--------------------------------------------------------------------------
X = bf.bf(:, 2:end); 
X = X - mean(X, 1);
X = X(:, any(X, 1)); % Remove empty columns

%% Modulation of connectivity parameters
%==========================================================================

% Define modulation functions
%--------------------------------------------------------------------------
fii = @(t) -4 * t .* exp(-0.7 * t) / gamma(2); % Self-inhibition
fij = @(t) exp(-(t - 3).^2 / 4) / gamma(2) / 2; % Forward connectivity

% Generate modulations
%--------------------------------------------------------------------------
uii = fii(bf.t); % Self-inhibition modulation
uii = uii - flip(uii);
uii = (uii - mean(uii)) * 1;

uij = fij(bf.t); % Forward connectivity modulation
uij = (uij - mean(uij)) * 4;

% Project modulations onto basis set
%--------------------------------------------------------------------------
Bii = X \ uii; % Self-inhibition coefficients
Bij = X \ uij; % Forward connectivity coefficients

%% Synthetic data generation
%==========================================================================

% Generate condition-specific modulations
%--------------------------------------------------------------------------
Btrue = cell(1, nconditions);
[Btrue{:}] = deal(zeros(nr, nr));
for i = 1:nconditions
    Btrue{i}(1, 1) = uii(i); % Self-inhibition
    Btrue{i}(2, 1) = uij(i); % Forward connectivity
end

% Add noise to modulations
%--------------------------------------------------------------------------
Bvar = exp(-1 + cos(2 * 0.1 * pi * bf.t)) * 2; % Variance
Bnois = Btrue;
for i = 1:nconditions
    Bnois{i}(1, 1) = uii(i) + normrnd(0, Bvar(i)) / sqrt(nsim);
    Bnois{i}(2, 1) = uij(i) + normrnd(0, Bvar(1)) / sqrt(nsim);
end

% Fill in generating B matrix
% -------------------------------------------------------------------------
Bval = cell(1, size(X, 2)); 
[Bval{:}] = deal(zeros(nr, nr));
for i = 1:size(X, 2)
    Bval{i}(1, 1) = Bii(i);
    Bval{i}(2, 1) = Bij(i); 
end

% Modulated parameters 
%--------------------------------------------------------------------------
Btemplate = [1 0; 1 0]; 

%% Model specification - a simple two source model with two electrodes
%==========================================================================

rng('default')

% specify model
%==========================================================================
 
% Setup model options
%--------------------------------------------------------------------------
Nc    = 2;  % number of channels
Ns    = 2;  % number of sources
options.spatial  = 'LFP'; % model source-level data
options.model    = 'CMC'; % use canonical microcircuit
options.analysis = 'ERP'; % ERP analysis

% Setup source model
%--------------------------------------------------------------------------
M = [];
M.dipfit.model   = options.model;
M.dipfit.type    = options.spatial;
M.dipfit.Nc      = Nc;
M.dipfit.Ns      = Ns;

% Repeat B matrix for the number of ERPs
%--------------------------------------------------------------------------
B = repelem({Btemplate}, nconditions);

% Get priors 
%--------------------------------------------------------------------------
[pE,pC] = spm_dcm_neural_priors(A,B,C,options.model);   % neuronal priors
[gE,gC] = spm_L_priors(M.dipfit);                       % spatial  priors
[x,f,h] = spm_dcm_x_neural(pE,options.model);

% hyperpriors (assuming a high signal to noise)
%--------------------------------------------------------------------------
hE      = 6;
hC      = 1/64;

% create model
%--------------------------------------------------------------------------
M.IS   = 'spm_gen_erp';
M.G    = 'spm_lx_erp';
M.FS   = 'spm_fy_erp'; 
M.f    = f;
M.x    = x;
M.h    = h; 
M.pE   = pE;
M.pC   = pC;
M.gE   = gE;
M.gC   = gC;
M.hE   = hE;
M.hC   = hC;
M.m    = size(C, 2);
M.n    = length(spm_vec(M.x));
M.l    = Nc;
M.ns   = 128;

% create input structure
%--------------------------------------------------------------------------
dt     = 2/1000;
pst    = (1:M.ns)*dt;
M.ons  = 64;
M.dur  = 16;
U.dt   = dt;
U.X    = eye(nconditions);

% specified true connectivity (P) and spatial parameters (G) - with
% condition specific effects on the intrinsic connectivity of the first
% source and its forward extrinsic connection
%--------------------------------------------------------------------------
P      = pE;
G      = gE;
P.B    = Btrue; 

%% Simulate ERPs
%==========================================================================

% generate neuronal response and data
%--------------------------------------------------------------------------
x     = spm_gen_erp(P,M,U);                 % neuronal response
L     = spm_lx_erp(G,M.dipfit);             % lead field

% Add noise
%--------------------------------------------------------------------------
s = {}; y = {}; 
V     = spm_sqrtm(spm_inv(spm_Q(1/2,M.ns))); % square root of noise covariance
snr = [];
for i = 1:length(x)
    n    = exp(-hE/2)*V*randn(M.ns,Nc);     % noise
    s{i} = spm_fy_erp(x{i}*L');             % signal
    y{i} = s{i} + n;                        % data (signal plus noise)
    snr = [snr; mean(s{i}.^2, 'all') / mean(n.^2, 'all')];
end

% Report SNR
%--------------------------------------------------------------------------
fprintf('SNR: %.2fdB\n', 10*log10(mean(snr)));

% Plot generated (noise-less) ERPs
%--------------------------------------------------------------------------
spm_figure('GetWin', 'Simulated and recovered ERPs');
clf; 
subplot(3,1,1)
plot_erps(pst, s, [-0.2, 0.4]); 
title('Generated')

% Plot generated (noisy) ERPs
%--------------------------------------------------------------------------
subplot(3,1,2)
plot_erps(pst, y, [-0.3, 0.5]);
title('Observed')

%% Recover parameter estimates
%==========================================================================

% Reset priors
%--------------------------------------------------------------------------
B = repelem({Btemplate}, size(X, 2));
[M.pE,M.pC] = spm_dcm_neural_priors(A,B,C,options.model);
[M.gE,M.gC] = spm_L_priors(M.dipfit);
[M.x,M.f,M.h]   = spm_dcm_x_neural(M.pE,options.model);

% Setup model input
%--------------------------------------------------------------------------
xU = U; 
xU.X = X;

% Setup data
%--------------------------------------------------------------------------
xY = []; 
xY.y = y; 
xY.nt = ones(nconditions, 1) * nsim; 
xY.Q = {spm_Q(1/2,M.ns,1)};

% Increase prior covariance
%--------------------------------------------------------------------------
for i = 1:size(X, 2)
    M.pC.B{i} = M.pC.B{i} * 8; 
end

% % Setup nice initial state to speed up inversion
% %------------------------------------------------------------------------
P   = M.pE;
P.B = spm_unvec(normrnd(0,1,size(spm_vec(M.pE.B))), M.pE.B); 
% G   = M.gE;
M.P = P; 
% M.Q = G;

% EM: inversion
%==========================================================================
[Qp,Qg,Cp,Cg,Ce,F,LE] = spm_nlsi_N(M,xU,xY);


%% Get results in the time domain
%==========================================================================

% Get an upsampled basis set
%--------------------------------------------------------------------------
bfup    = bf; 
bfup.dt = bf.dt / 10; 
bfup    = spm_get_bf(bfup); 

% Setup an upsampled time-vector (for plotting)
%--------------------------------------------------------------------------
bfup.t = (0:bfup.dt:bfup.length)'; 
bfup.t = bfup.t(1:end-1); 

% Get the upsampled design matrix
%--------------------------------------------------------------------------
Xup = bfup.bf(:, 2:end); 

% Mean- and scale-correct the design matrix to match the previous one
%--------------------------------------------------------------------------
Xup = Xup - mean(X, 1); 
Xup = Xup ./ (max(Xup, [], 1) - min(Xup, [], 1)); 
Xup = Xup .* (max(X, [], 1) - min(X, [], 1));

% Project the estimated posteriors onto the basis set
%--------------------------------------------------------------------------
[Bproj, Cbproj] = spm_dcm_project_B(M, Qp.B, Cp, Xup); 

% Upsample the modulation of self-inhibition
%--------------------------------------------------------------------------
Uii = fii(bfup.t);
Uii = Uii - flip(Uii); 
Uii = (Uii - mean(Uii)) * 1; 

% Upsample the modulation of forward connectivity
%--------------------------------------------------------------------------
Uij = fij(bfup.t); 
Uij = (Uij - mean(Uij))  * 4; 

%% Plot the results 
%==========================================================================

% Rescale time vector
%--------------------------------------------------------------------------
t = bfup.t / bfup.t(end) * ((nconditions-1)*M.ns*dt + pst(end));

% Setup plot
%--------------------------------------------------------------------------
color = [ html2rgb('cb5eee'); html2rgb('4be1ec')]; 

spm_figure('GetWin', 'Original and recovered time-varying connectivity');
clf;
subplot(2,1,1)
% axes('NextPlot','replacechildren', 'ColorOrder',color);
hold on; 

% Compute std. dev. band (self-inhibition)
%--------------------------------------------------------------------------
low  = Bproj{1} - 1.*sqrt(diag(Cbproj{1,1}));
high = Bproj{1} + 1.*sqrt(diag(Cbproj{1,1}));

% Plot std. dev. band (self-inhibition)
%--------------------------------------------------------------------------
fill([t; flip(t)],[low; flipud(high)], ...
    color(1, :), 'FaceAlpha', 0.15, 'EdgeAlpha', 0);

% Plot recovered and original TVEC (self-inhibition)
%--------------------------------------------------------------------------
plot(t, Bproj{1}, 'Color', color(1, :), 'LineWidth', 2)
plot(t, Uii , '--', 'Color',  html2rgb('9947b4'), 'LineWidth', 2) 

% Configure plot
%--------------------------------------------------------------------------
ax = gca;
ax.XAxisLocation = "origin";
ax.YAxisLocation = "origin";
ylim([-3 4]);
xlim([0, t(end)])
xlabel('Time (s)')
axis tight; 
ylabel('Log scaling factor')
legend('', 'True', 'Recovered', ...
    'Location', 'southoutside', 'Orientation', 'horizontal')
legend('boxoff')
title('Modulation of self-inhibition')


spm_figure('GetWin', 'Original and recovered time-varying connectivity');
subplot(2,1,2)
hold on; 

% Compute std. dev. band (forward connection)
%--------------------------------------------------------------------------
low  = Bproj{2} - 1.*sqrt(diag(Cbproj{2,2})); 
high = Bproj{2} + 1.*sqrt(diag(Cbproj{2,2}));

% Plot std. dev. band (forward connection)
%--------------------------------------------------------------------------
fill([t; flip(t)], [low; flipud(high)], color(2, :), ...
            'FaceAlpha', 0.15, 'EdgeAlpha', 0);

% Plot recovered and original TVEC (forward connection)
%--------------------------------------------------------------------------
plot(t, Bproj{2}, 'Color', color(2, :), 'LineWidth', 2)
plot(t, Uij , '--', 'Color',html2rgb('3fbdc6'), 'LineWidth', 2) 

% Configure plot
%--------------------------------------------------------------------------
ax = gca;
ax.XAxisLocation = "origin";
ax.YAxisLocation = "origin";
ylim([-3 4]);
xlim([0, t(end)])
xlabel('Time (s)')
axis tight; 
ylabel('Log scaling factor')
legend('', 'True', 'Recovered', ...
    'Location', 'southoutside', 'Orientation', 'horizontal')
legend('boxoff')
title('Modulation of forward connection')

%% Plot the ERPs
%==========================================================================

% generate neuronal response and data
%--------------------------------------------------------------------------
x     = spm_gen_erp(Qp,M,xU);                % neuronal response
L     = spm_lx_erp(Qg,M.dipfit);             % lead field

r = {}; 
for i = 1:length(x)
    r{i} = spm_fy_erp(x{i}*L');             % signal
end

% generate neuronal response and data
%--------------------------------------------------------------------------
spm_figure('GetWin', 'Simulated and recovered ERPs');
subplot(3,1,3)
plot_erps(pst, r, [-0.2, 0.4], s)
title('Recovered')

end

function plot_erps(pst, x, yrange, g)
    nconditions = numel(x); 
    ns = length(pst); 
    dt = mean(diff(pst)); 
    hasgroundtruth = nargin > 3; 

    % Setup colors
    %----------------------------------------------------------------------
    colors = [ html2rgb("eddd8c"); html2rgb("8ADD94")]; 

    for i = 1:nconditions
        % Plot ERP stick
        %------------------------------------------------------------------
        scatter(64/dt/1e6 + (i-1)*ns*dt, ...
            yrange(1) + 0.1, 15, 'k', '^', 'filled');
        if i == 1
            hold on; 
        end
        plot([0 0] + 64/dt/1e6 + (i-1)*ns*dt, ...
            [yrange(1) yrange(1) + 0.1], '-k', 'LineWidth', 1)
        
        % Plot ERP
        %------------------------------------------------------------------
        e = x{i}; 
        plot(pst + (i-1)*ns*dt, e(:, 1), '-', 'Color',colors(1,:), 'LineWidth', 2)
        plot(pst + (i-1)*ns*dt, e(:, 2), '-', 'Color',colors(2,:), 'LineWidth', 2)
    
        if hasgroundtruth
            % Plot ground truth
            %--------------------------------------------------------------
            plot(pst + (i-1)*ns*dt, g{i}, [':k'], 'LineWidth', 1.5)
        end
    end

    % Setup legend
    %----------------------------------------------------------------------
    nplots = 4;
    if hasgroundtruth
        nplots = nplots + 1;
    end
    lgd = cell(nconditions*nplots,1);
    [lgd{:}] = deal(''); 
    lgd{3} = 'LFP1'; 
    lgd{4} = 'LFP2'; 
    if hasgroundtruth
        lgd{5} = 'True'; 
    end
    legend(lgd{:}, 'Location', 'northeast', 'Orientation', 'horizontal')
    legend('boxoff')

    % Setup axes 
    %----------------------------------------------------------------------
    xlabel('Time (s)')
    ylabel('Amplitude (a.u.)')
    ylim(yrange)
    xlim([pst(1) pst(end) + (i-1)*ns*dt])
end


function r = html2rgb(c)
    r = sscanf(c, '%2x%2x%2x', [1 3]) / 255;
end