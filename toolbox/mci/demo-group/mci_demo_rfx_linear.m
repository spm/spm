
clear all
close all

disp('Linear RFX example');

% Number of subjects
N=8;

% Design matrix 
des='dct';
%des='quadratic';
%des='linear-offset';

% Time points per subject
Nobs=20;

% Observation noise precision
sd=0.2;
lambda=1/(sd^2);

% Population-level parameters
[M{1},U{1}] = mci_linear_struct (Nobs,lambda,des);
w_pop=spm_normrnd(M{1}.pE,M{1}.pC,1);

% Subject parameters and data
for n=1:N,
    [M{n},U{n},Xfull] = mci_linear_struct (Nobs,lambda,des);
    w_true(:,n) = spm_normrnd(w_pop,M{n}.pC,1);
    Y{n}.y = U{n}.X*w_true(:,n)+sqrt(M{n}.Ce)*randn(Nobs,1);
end

% PEB analysis
do_PEB=1;
if do_PEB
    disp('Parametric Empirical Bayes ...');
    P{1}.X=[];y=[];
    for n=1:N,
        P{1}.X=blkdiag(P{1}.X,U{n}.X);
        y=[y;Y{n}.y];
    end
    Np=length(M{1}.pE);
    Nx=size(P{1}.X,1);
    P{1}.C{1}=eye(Nx);
    P{2}.X=kron(ones(N,1),eye(Np));
    % Prior variance components for each parameter
    for j=1:Np,
        z=zeros(Np,Np);
        z(j,j)=1;
        P{2}.C{j}=kron(eye(N),z);
    end
    tic;
    [C,P,F] = spm_PEB(y,P);
    toc
    m_PEB=C{3}.E; % Population - posterior mean
    C_PEB=C{3}.C; % Population - posterior covariance
    h_PEB=C{2}.h; % Estimate prior vars
    w_PEB=reshape(C{2}.E,Np,N); % Subject effects
    names={'True','MCI','PEB'};
else
    names={'True','MCI'};
end

% MCI
disp('MCI ...');
MCI.verbose=0;
MCI.M=M; MCI.U=U; MCI.Y=Y;
MCI.update_obs_noise=1;

% Assume all variables are random effects
Nrand=length(M{1}.pE);
S.prior.P=Nrand;
S.prior.a=Nrand/2;
S.prior.B=eye(Nrand);
S.prior.beta=1;
S.prior.m=zeros(Nrand,1);
S.N=length(M);
MCI.S=S;
    
tic;
MCI = spm_mci_mfx (MCI);
toc

% Population-level estimates
dist{1}.P=MCI.sm;
dist{1}.Ep=MCI.sm_mean;
dist{1}.ind=MCI.post_ind;
dist{1}.type='sample';
dist{1}.color='k';
dist{1}.names=U{1}.names;
mci_plot_dist_multi (dist,'MCI: Pop-Level',w_pop);

% Subject-level estimates
h=figure;
set(h,'Name',sprintf('Subject Level'));
for j=1:Np,
    subplot(Np,1,j);
    plot(squeeze(MCI.sw(j,:,:))');
    grid on
    xlabel('RFX iteration');
    ylabel(sprintf('w(%d)',j));
end

if do_PEB
    dist{1}.Ep=m_PEB;
    dist{1}.Cp=C_PEB;
    dist{1}.type='Gaussian';
    dist{1}.color='b';
    dist{1}.names=U{1}.names;
    mci_plot_dist_multi (dist,'PEB: Pop-Level',w_pop);
end

mci_plot_noiseSD (MCI.Ce,MCI.post_ind);

plot_model_fits=0;
if plot_model_fits
    % Plot model fits
    Nplot=min(N,8);
    for n=1:Nplot,
        yt{1}=Xfull*w_true(:,n);
        yt{2}=Xfull*MCI.sw_mean(:,n);
        if do_PEB
            yt{3}=Xfull*w_PEB(:,n);
        end
        h=figure;
        set(h,'Name',sprintf('Subject %d',n));
        mci_linear_plot_fit (M{n},Y{n}.y,yt,names);
    end
end
