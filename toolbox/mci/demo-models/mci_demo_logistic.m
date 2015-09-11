
clear all
close all

disp('Logistic regression model');

%[M,U,Y] = mci_logistic_struct ('ripley');
%[M,U,Y] = mci_logistic_struct ('pima');

[M,U] = mci_logistic_struct ('dct');
w_true=[0.2,8,6]';
[tmp,Y.y] = mci_logistic_gen(w_true,M,U);

%mcmc.inference='amc';
%mcmc.inference='vl';
mcmc.inference='langevin';

mcmc.verbose=0;
mcmc.maxits=1024;

post = spm_mci_post (mcmc,M,U,Y);

% q=post.quantiles;
% Np=size(q,1);
% figure
% lw=2;
% errorbar([1:Np],q(:,2),q(:,1),q(:,3),'k','LineWidth',lw);
% set(gca,'FontSize',16);
% grid on
% xlabel('Parameter');
% title('Posterior median and 95% intervals');

if ~strcmp(mcmc.inference,'vl')
    diag.traceplot=1;
    diag.eplot=1;
    if strcmp(mcmc.inference,'langevin')
        diag.bplot=1;
    else
        diag.bplot=0;
    end
    spm_mci_diag(post,diag);
    
    stats = spm_mci_mvnpost (post,'ESS')
    stats = spm_mci_mvnpost (post,'thinning')
    
    spm_mci_quantiles (post,2);
end
