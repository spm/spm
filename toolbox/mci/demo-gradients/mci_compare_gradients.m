function [els,names] = mci_compare_gradients (model,cost,methods)
% Compare methods for gradient computation
% FORMAT [els,names] = mci_compare_gradients (model,cost,methods)
%
% model     'phase', 'nmm-r2p2'
% cost      'loglike', 'spm_mci_joint' (default)
% methods   vector of integers indicating which methods to
%           compare eg. [1,2,3,4,5] (default) for 1. SensMat, 
%           2. SensSun, 3. AdjMat, 4. AdjSun, 5. FD
%
% els       Computation times
% names     Names of compared methods
%
% Note: 4. AdjSun may not work for nmm2-r2p2.
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_compare_gradients.m 6548 2015-09-11 12:39:47Z will $

if nargin < 2 | isempty(cost)
    cost='spm_mci_joint';
end

if nargin < 3 | isempty(methods)
    methods=[1:5];
end

for i=1:5,
    if length(intersect(methods,i))>0
        m(i)=1;
    else
        m(i)=0;
    end
end

% Set dp for Finite Differences
sdp=1e-1;

[P,M,U,Y,ind] = mci_compare_setup (model);

mci_plot_outputs(M,Y);

M.reltol=1e-2;
M.abstol=1e-4;

M = spm_mci_minit (M);
if isstruct(M.pC)
    pC=full(diag(spm_vec(M.pC)));
else
    pC = M.pC;
end

% Evalute gradients at new point
P=spm_normrnd(M.vpE,pC,1);

% Parameters in reduced space
Pr = M.V'*(P-M.vpE);

names={'SensMat','SensSun','AdjMat','AdjSun','FD'};
names=names(methods);
cm=1;

if strcmp(cost,'loglike')   
    
    if m(1)
        M.int='ode15';
        tic;
        [dLdp,tmp,st] = spm_mci_glike_deriv (P,M,U,Y);
        els(cm)=toc;
        cm=cm+1;
    end
    
    if m(2)
        M.int='sundials';
        tic;
        [dLdp_sun,tmp,st] = spm_mci_glike_deriv (P,M,U,Y);
        els(cm)=toc;
        cm=cm+1;
    end
    
    if m(3)
        M.adjlike=1;
        tic;
        [dLdp_adj] = spm_mci_adjoint (Pr,M,U,Y);
        els(cm)=toc;
        cm=cm+1;
    end
    
    if m(4)
        %M.backint=1;
        tic;
        [dLdp_adj_sun] = spm_mci_adjoint_sun (Pr,M,U,Y);
        els(cm)=toc;
        cm=cm+1;
    end
    
    if m(5)
        tic;
        P0=P;
        L0=spm_mci_glike(P0,M,U,Y);
        
        for i=1:M.Np,
            Pd=P0;
            dp=sdp*Pd(i);
            Pd(i)=Pd(i)+dp;
            Ld = spm_mci_glike (Pd,M,U,Y);
            dLdp_finite(i)=(Ld-L0)/dp;
        end
        els(cm)=toc;
    end
    
else
    if m(1)
        tic;
        M.int='ode15';
        [dLdp,tmp,st] = spm_mci_joint_grad (Pr,M,U,Y);
        els(cm)=toc;
        cm=cm+1;
    end
    
    if m(2)
        tic;
        M.int='sundials';
        [dLdp_sun,tmp,st] = spm_mci_joint_grad (Pr,M,U,Y);
        els(cm)=toc;
        cm=cm+1;
    end
    
    if m(3)
        tic;
        [dLdp_adj] = spm_mci_adjoint (Pr,M,U,Y);
        els(cm)=toc;
        cm=cm+1;
    end
    
    if m(4)
        tic;
        [dLdp_adj_sun] = spm_mci_adjoint_sun (Pr,M,U,Y);
        els(cm)=toc;
        cm=cm+1;
    end
    
    if m(5)
        tic;
        P0=Pr;
        L0=spm_mci_joint(P0,M,U,Y);
        for i=1:M.Np,
            Pd=P0;
            dp=sdp*Pd(i);
            Pd(i)=Pd(i)+dp;
            Ld = spm_mci_joint (Pd,M,U,Y);
            dLdp_finite(i)=(Ld-L0)/dp;
        end
        els(cm)=toc;
    end
end

if st==-1
    disp('Problem with integration');
    return
end
    
figure;
if m(1) plot(dLdp(ind)); end
hold on
if m(2) plot(dLdp_sun(ind),'c'); end
if m(3) plot(dLdp_adj(ind),'r'); end
if m(4) plot(dLdp_adj_sun(ind),'m'); end
if m(5) plot(dLdp_finite(ind),'k'); end
xlabel('p');
ylabel('dLdp');
grid on
legend(names);

figure
bar(els,'k');
grid on
ylabel('Seconds');
set(gca,'XTickLabel',names);

disp(' ');
if isfield(M,'dfdx')
    disp('State Jacobian M.dfdx is specified for this model');
else
    disp('State Jacobian M.dfdx not specified for this model');
end
if isfield(M,'dfdp')
    disp('Parameter Jacobian M.dfdp is specified for this model');
else
    disp('Parameter Jacobian M.dfdp not specified for this model');
end

disp(' ');
disp('SensMat, AdjMat and AdjSun will use');
disp('M.dfdx and M.dfdp if these fields are specified.');
