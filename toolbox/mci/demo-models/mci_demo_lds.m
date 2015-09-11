
clear all
close all

disp('Synthetic data from d-region LDS');
disp('Model has constrained connectivity');
disp(' ');

% Number of regions
M.d=4;

% Between region connectivity [to, from]
M.name='forward';
%M.name='backward';
%M.name='bidirectional';

% Observation noise SD
M.sd=0.1;

% Initial states
M.R=linspace(3,1.5,M.d)';
M.drop=0.5;
%M.t=[1:400]'/4;  
M.t=[1:25]'/0.25;  
    
[M,U,names] = mci_lds_struct (M);

Pt = [-0.04,-0.01,-0.005,-0.01,0.01,0.01,0.01]';
P = mci_lds_par2lat (Pt,M);

%P = lds_params (M,U);
turn_off_connectivity=0;
if turn_off_connectivity
    P(M.d+1:end)=0;
end
Y = mci_lds_gen (M,U,P);

figure;
plot(M.t,Y);

%mcmc.inference='amc';
%mcmc.inference='vl';
mcmc.inference='langevin';

mcmc.maxits=1024;
mcmc.verbose=0;

post = spm_mci_post (mcmc,M,U,Y);

disp('True dynamics:');
[f,A] = mci_lds_fx (M.x0,U,P,M);
disp(A);

disp('Estimated dynamics:');
[f,A] = mci_lds_fx (M.x0,U,post.Ep,M);
disp(A);

plot_prior=0;
if strcmp(mcmc.inference,'langevin')
    scale=100;
    dist{1}=post;
    dist{1}.color='k';
    dist{1}.names=names;
    S=size(dist{1}.P,2);
    for s=1:S,
        Pt_post(:,s)=mci_lds_lat2par (dist{1}.P(:,s),M);
    end
    dist{1}.P=Pt_post*scale;
    
    if plot_prior
        Ns=1000;
        P=spm_normrnd(M.pE,M.pC,Ns);
        dist{2}.ind=[1:Ns];
        dist{2}.color='b';
        dist{2}.names=dist{1}.names;
        dist{2}.type='sample';
        for s=1:Ns,
            Pt_prior(:,s)=mci_lds_lat2par (P(:,s),M);
        end
        dist{2}.P=Pt_prior*scale;
    end
    mci_plot_dist_multi (dist,'Flow',Pt*scale);
end

switch mcmc.inference
    case 'langevin',
        ind=[1:size(post.Ce,3)];
        mci_plot_noiseSD (post.Ce,ind);
    case 'vl',
        disp('Error SD:');
        disp(sqrt(diag(post.Ce)));
end

stats = spm_mci_mvnpost (post,'ESS')
stats = spm_mci_mvnpost (post,'thinning')
