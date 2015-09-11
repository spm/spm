
clear all
close all

disp('Logistic RFX example');

% Number of subjects
N=16;

[M1,U1] = mci_logistic_struct ('dct');
w_pop=[0.2,8,6]';
Np=length(w_pop);

% [g,y]=mci_logistic_gen(P1,M1,U1);
% figure;plot(g);
% figure;imagesc(U1.X);
% return

w_true=spm_normrnd(w_pop,0.001*eye(3),N);

% Subject parameters and data
for n=1:N,
    M{n}=M1;
    U{n}=U1;
    [tmp,Y{n}.y] = mci_logistic_gen(w_true(:,n),M1,U1);
end

% MCI
disp('MCI ...');
MCI.verbose=0;
MCI.M=M; MCI.U=U; MCI.Y=Y;

MCI.total_its=1024;
MCI.rinit=0.25;

MCI.inference='lgv';
%MCI.inference='amc';

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
grid on

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

SE=(MCI.sw_mean-w_true).^2;
RMSE=sqrt(mean(SE,2));
aRMSE=mean(RMSE);
disp(sprintf('Average RMSE over parameters=%1.2f',aRMSE));
