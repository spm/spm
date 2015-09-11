
clear all
close all

disp('Random Effects demo with');
disp('Phase Coupling Models');
disp(' ');
disp('Initial conditions are known');
disp('Flow parameters are treated as random effects');

disp('Synthetic data from d-region weakly coupled oscillator');
disp('model with reduced connectivity');

% Dimension of state space
d=3;

conn(1).regions=[1,2];
conn(1).name='vis2mtl';
conn(2).regions=[3,2];
conn(2).name='vis2ifg';

[M1,U1] = mci_rphase_struct(d,conn);
M1.freq=6;
M1.Npflow=M1.Np;
M1.Npout=0;
M1.x0=[-pi/2,0,pi/2]';

% Observation noise SD
sd=0.01;

% Number of subjects
N=12;
m=[-0.6,-0.6,0,0];
Np=length(m);
C=0.1^2*eye(Np);
w_true = spm_normrnd(m,C,N);
for n=1:N,
    M{n}=M1;U{n}=U1;
    wn=w_true(:,n);
    pars.aconn=[wn(1), wn(2)];
    pars.bconn=[wn(3), wn(4)]';
    %pars.bconn=[0,0]';
    P=spm_vec(pars);
    Y{n}.y = mci_rphase_gen(P,M{n},U{n});
end

% Assign init, flow and output parameters
assign.init_par='known';
assign.flow_par='random';
assign.out_par='known';
MCI.pout0=[]; % there are no 'output' parameters
MCI.assign=assign;

MCI.M=M; MCI.U=U; MCI.Y=Y;
MCI.verbose=1;
MCI.total_its=16;
MCI.rinit=0.25;

S.prior.P=Np;
S.prior.a=Np/2;
S.prior.B=M{1}.pC;
S.prior.beta=1;
S.prior.m=zeros(Np,1);
S.N=N;
MCI.S=S;
    
tic;
MCI = spm_mci_mfx_dynamic (MCI);
toc

disp('Mean true params:');
mwt=mean(w_true,2)

disp('MCI second level posterior:');
MCI.sm_mean

Np=M{1}.Npflow;
for j=1:Np,
    h=figure;
    set(h,'Name',sprintf('w(%d)',j));
    
    min_w=min(w_true(j,:));
    max_w=max(w_true(j,:));
    
    plot(w_true(j,:),MCI.sw_mean(j,:),'kx','MarkerSize',10);
    hold on
    plot([min_w max_w],[min_w max_w],'k-','LineWidth',2);
    set(gca,'FontSize',18);
    xlabel('True');
    ylabel('Estimated');
    grid on
end

e=MCI.sw_mean-w_true;
sse2=trace(e'*e);
disp(sprintf('Final first level SSE=%1.2f',sse2));

figure;
plot(MCI.sm');
grid on
xlabel('RFX iteration');
ylabel('Population Level Parameters');

h=figure;
set(h,'Name',sprintf('Subject Level'));
for j=1:Np,
    subplot(Np,1,j);
    plot(squeeze(MCI.sw(j,:,:))');
    grid on
    xlabel('RFX iteration');
    ylabel(sprintf('w(%d)',j));
end


