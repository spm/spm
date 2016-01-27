function [M,stats] = spm_mci_lgv_vl (mcmc,M,U,Y,vl,beta)
% Sampling using Langevin Monte Carlo on path from VL solution
% FORMAT [M,stats] = spm_mci_lgv_vl (mcmc,M,U,Y,vl,beta)
%
% mcmc  Sampling parameters
%       .verbose            display progress
%       .maxits             maximum number of total samples 
%       .init               initial sample values (start of chain)
%       .h                  step size
%
% M     Model Structure
%       .dL                 Gradients and curvatures are computed using 
%                           this user-specified function. If this is absent
%                           they will be computed using (i) the forward
%                           sensitivity method for dynamical models 
%                           (ie. if M.f exists) or (ii) finite differences
%                           otherwise
%                           
% U     Inputs
% Y     Data
% vl    Variational Laplace solution
%       .Ep                 Posterior Mean
%       .Cp                 Posterior Covariance
% beta  Inverse Temperature (0 at VL solution, 1 at posterior)
%
% M     Updated model structure
% stats Structure with fields:
%
% .P     Samples, [maxits x M.Np] 
% .E     Negative log joint prob, [maxits x 1]
%
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: spm_mci_lgv_vl.m 6697 2016-01-27 14:57:28Z spm $

% Defaults
try, verbose=mcmc.verbose; catch, verbose=0; end
try, maxits=mcmc.maxits; catch, maxits=64; end
try, h=mcmc.h; catch, h=0.5; end 
try, init=mcmc.init; catch, init=spm_vec(M.pE); end

% Compute eigen-parameterisation
M = spm_mci_minit (M);
V  = M.V;

% Initial param in eigenspace
xinit = M.V'*(init-M.vpE);
Nr = length(xinit);

% Sample matrix
x = zeros(maxits,M.Np);
x(1,:) = xinit';         

if verbose figure; end
    
acc=zeros(maxits,1);

% Main Loop
i=1; Ce=[];
while (i <= maxits),
    
    if verbose
        if mod(i,plot_int) == 0 & i > 2
            spm_mci_progress (x,E,i);
        end
    end
    
    % Proposal (first proposal always accepted)
    if i==1
        pos=x(1,:)';
        curr=[];
    else
        pos=spm_normrnd(curr.mu,curr.Cp,1);
    end
        
    % Quantities re proposal
    prop.pos=pos;
    [j,iCpY,st,prop.L,prop.L2] = spm_mci_joint_grad (pos,M,U,Y);
    if st==-1
        disp('Integration problem in spm_mci_lgv.m');
        keyboard
    end
    
    % Compute gradient for this temperature
    pos_orgspace=V*pos+M.vpE;
    j2 = (1-beta)*vl.Lambdap*(vl.Ep-pos_orgspace);
    %j2 = (1-beta)*vl.Lambdap*(pos_orgspace-vl.Ep);
    j = beta*j + j2';
    
    % Covariance for this temperature
    %prop.Cp = (h/beta)*eye(Nr);
    %prop.Cp = (h/beta)*inv(iCpY+M.ipC);
    prop.Cp = h*inv((1-beta)*vl.Lambdap+beta*iCpY+beta*M.ipC);
    
    % Posterior mean under local linear approximation
    prop.mu = pos+0.5*h*prop.Cp*j(:);
    
    prop.iCp = inv(prop.Cp);
    prop.logdetCp = spm_logdet(prop.Cp);
    
    % Accept proposal ?
    [curr,accepted,bayes_fb(i),dL(i)] = spm_mci_mh_update(curr,prop,verbose);
    acc(i)=accepted;
    E(i) = -curr.L;
    L2(i) = curr.L2;
    if i > 1
        dEdit(i-1)=100*(E(i)-E(i-1))/E(i-1);
    end
    x(i,:)=curr.pos;
    
    i=i+1;
end

if verbose
    disp(sprintf('Total accepted samples = %d', sum(acc)));
end

% Project parameters back from eigenspace into original space
x=x(1:i-1,:);
nj=size(x,1);
stats.P=M.vpE*ones(1,nj)+V*x';

stats.E=E;
stats.dEdit=dEdit;
stats.acc=acc;
stats.Ce=Ce;
stats.bayes_fb=bayes_fb;
stats.dL=dL;
stats.L2=L2;
