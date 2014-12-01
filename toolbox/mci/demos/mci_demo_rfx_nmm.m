
clear all
close all

disp('Random Effects demo with');
disp('2 region Neural Mass Models');
disp(' ');
disp('Initial conditions are known');
disp('Flow parameters are treated as random effects');

% Backward connection ?
back=1;

% Observation noise SD
sd=0.01;

% Number of subjects
N=3;
m=[1,1]';
C=0.1^2*eye(2);
w_true = spm_normrnd(m,C,N);
for n=1:N,
    [M{n},U{n}] = mci_nmm_struct(back,sd,2);
    P=w_true(:,n);
    Y{n}.y = mci_nmm_gen(M{n},U{n},P);
    %[M{n},U{n},Y{n}] = nmm_r2_init (back,sd,w_true(:,n));
end

% Assign init conditions as known
% Assign flow as random effects
assign.init_par='known';
assign.flow_par='random';
MCI.assign=assign;

MCI.verbose=1;
MCI.M=M; MCI.U=U; MCI.Y=Y;
MCI.update_rfx=1;
MCI.update_ffx=0;

tic;
MCI = spm_mci_mfx (MCI);
toc

% MCI
MCI.M=M;MCI.U=U;MCI.Y=Y;
MCI.verbose=1;

disp('Mean true params:');
mwt=mean(w_true,2)

disp('MCI second level posterior:');
MCI.sm_mean

for j=1:d,
    h=figure;
    set(h,'Name',sprintf('w(%d)',j));
    
    min_w=min(w_true(j,:));
    max_w=max(w_true(j,:));
    
    plot(w_true(j,:),MCI.w_init(j,:),'r.');
    hold on
    plot(w_true(j,:),MCI.sw_mean(j,:),'k.');
    plot([min_w max_w],[min_w max_w],'k-','LineWidth',2);
    set(gca,'FontSize',18);
    xlabel('True');
    ylabel('Estimated');
    grid on
    
end

e=MCI.w_init-w_true;
sse1=trace(e'*e);
e=MCI.sw_mean-w_true;
sse2=trace(e'*e);

disp(sprintf('Initial first level SSE=%1.2f',sse1));
disp(sprintf('Final first level SSE=%1.2f',sse2));

figure;
plot(MCI.sm');
grid on
title(sprintf('%s Inference: Population Level',MCI.sample));
xlabel('MFX iteration');
ylabel('Parameters');

h=figure;
set(h,'Name',sprintf('%s Inference: Subject Level',MCI.sample));
for j=1:d,
    subplot(d,1,j);
    plot(squeeze(MCI.sw(j,:,:))');
    grid on
    xlabel('RFX iteration');
    ylabel(sprintf('w(%d)',j));
end


