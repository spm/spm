function [log_prob,S,E] = mci_plot_surface (P,M,U,Y,S,dist)
% Plot log probability surface
% FORMAT [log_prob,S,E] = mci_plot_surface (P,M,U,Y,S,dist)
%
% P         Parameters
% M         Model structure
% U         Inputs
% Y         Data
% S         Surface data structure
% dist      'prior', 'like' or 'post'
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: mci_plot_surface.m 6548 2015-09-11 12:39:47Z will $

try 
    S.surf;
catch
    S.surf=1;
end

try 
    S.image;
catch
    S.image=1;
end

% Number of bins defining surface in each dimension
Nbins=S.Nbins;
pxy=S.pxy;

[bup,bdown] = meshgrid(pxy(1,:),pxy(2,:));
        
% For computing log prior term
M = spm_mci_priors (M);
M = spm_mci_minit (M);

% Which params to vary
P1str=[S.param{1},'=P1;'];
P2str=[S.param{2},'=P2;'];
        
[I,J]=size(bup);
mm=1;
for i=1:I,
    for j=1:J,
        %disp(sprintf('Model %d out of %d',mm,I*J));
        P1=bup(i,j);
        P2=bdown(i,j);
        eval(P1str);
        eval(P2str);
        
        % Get parameters in eigenspace of prior 
        Pv = spm_vec(P);
        M.vpE=spm_vec(M.pE);
        p = M.V'*(Pv-M.vpE);
        switch lower(dist),
            case 'post',
                log_prob(i,j)= spm_mci_joint (p,M,U,Y);
                
            case 'prior',
                log_prob(i,j)=-p'*M.ipC*p/2 + M.log_prior_t2;
                
            case 'like',
                log_prob(i,j)= feval(M.L,Pv,M,U,Y);
                
            otherwise
                disp('Unknown distribution type');
                return
        end
        mm=mm+1;
    end
end

if S.surf
    figure
    surf(bup,bdown,log_prob);
    set(gca,'FontSize',18);
    xlabel(S.name{1});
    ylabel(S.name{2});
    title(['log ',dist]);
end

if S.image
    figure;
    imagesc(pxy(1,:),pxy(2,:),log_prob);
    set(gca,'FontSize',18);
    axis xy
    hold on
    xlabel(S.name{1});
    ylabel(S.name{2});
    title(['log ',dist]);
end

S.x=bup;
S.y=bdown;
S.pxy=pxy;

