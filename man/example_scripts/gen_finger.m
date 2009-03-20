function [yy,P,M,U] = gen_finger (sim)
% Generate finger movement data
% FORMAT [yy,P,M,U] = gen_finger (sim)
%
% sim       Simulation data structure:
%
%           .Nt         number of trials
%           .m          first or 2nd order PIF
%           .init       'partial': initial phase diff
%                       restricted to small range
%                       'full': initial phase diff
%                       uniform in 0 to 2 pi
%           .noise_dev  STD of additive noise
%           .do_plot    plot data (1)
%
% yy        yy{n} for nth trial data
% P         model parameters
% M         model structure
% U         input structure

Nr=2; % Number of regions
Nu=0; % Number of between-trial effects

Nt=sim.Nt; % Number of trials

f=6;
fb=2;
fs=100;
dt=1/fs;
secs=1;
N=secs*fs;
U.tims=[1:N]'*dt;

a=0.5;
switch sim.m
    case 1,
        b=0;
    case 2,
        b=0.75*a;
end
P.As(:,:,1)=[0 0; a 0];
P.As(:,:,2)=[0 0; b 0];
        
P.B=[];
P.df=zeros(Nr,1);
P.L=[1 1]';

% Loop over trials
for n=1:Nt,
    % Initial phases
    switch sim.init
        case 'full',
            phi_0=2*pi*rand(Nr,1);
        case 'partial'
            phi_0(1)=2*pi*rand(1,1);
            delta_phi=4*rand(1,1)-2;
            phi_0(2)=phi_0(1)+delta_phi;
    end
    
    M.freq=f;
    M.fb=fb;

    % Initial phase
    trial{n}.x=phi_0;
    x(:,1)=phi_0;

    U.u=zeros(N,1);
    U.dt=dt;
    U.X=sparse(1,0);
    M.f='spm_fx_phase';
    M.G='spm_lx_phase';
    M.x=x(:,1);

    M.ns=N; % number of time samples
    
    % Get state variables
    states=spm_gen_phase(P,M,U);

    % Get observables
    G=spm_lx_phase(P,M);
    y{n}=(states{1}*G');
end

% Noisy data
for n=1:Nt,
    yy{n}=y{n}+sim.noise_dev*randn(N,Nr);
end

M.trial=trial;
M.x=zeros(Nr,1);

if sim.do_plot
    figure;
    for n=1:Nt,
        subplot(Nt,1,n);
        plot(sin(yy{n}));
        title(sprintf('Generated data, trial %d',n));
    end
end