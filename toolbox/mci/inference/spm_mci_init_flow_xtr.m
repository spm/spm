function [MCI] = spm_mci_init_flow_xtr (MCI)
% Extract init and flow parameters for dynamical models
% FORMAT [MCI] = spm_mci_init_flow_xtr (MCI)
%
% This routine adds the following fields:
%
% Initial states
% .pinit.P          [Nsamples x Ninit] samples
% .pinit.ind        indices for posterior (ie. excluding burn-in)
% .pinit.Ep         [Ninit x 1] posterior mean
% .pinitK = .pinit.Ep ??
%
% Flow parameters
% .pflow.P          [Nsamples x Nflow] samples
% .pflow.ind        indices for posterior (ie. excluding burn-in)
% .pflow.Ep         [Nflow x 1] posterior mean
%
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_init_flow_xtr.m 6275 2014-12-01 08:41:18Z will $

assign=MCI.assign;
sw=MCI.sw;sm=MCI.sm;v=MCI.sv';
post_ind=MCI.post_ind;
S=MCI.S;

if strcmp(assign.init_par,'random')
    MCI.pinitK=MCI.sw_mean(assign.w_init,:);
    for n=1:S.N,
        MCI.pinitK_cov(:,:,n)=cov(squeeze(sw(assign.w_init,n,post_ind))');
    end
    
    % Actual Samples (group level)
    MCI.pinit.P=sm(assign.w_init,:);
    MCI.pinit.ind=post_ind;
    MCI.pinit.Ep=mean(sm(assign.w_init,post_ind),2);
else
    MCI.pinitK=MCI.sv_mean(assign.v_init);
    MCI.pinitK_cov=cov(v(post_ind,assign.v_init));
    
    % Actual Samples:
    MCI.pinit.P=v(:,assign.v_init)';
    MCI.pinit.ind=post_ind;
    MCI.pinit.Ep=MCI.pinitK;
end
if strcmp(assign.flow_par,'random')
    MCI.pflowK=MCI.sw_mean(assign.w_flow,:);
    for n=1:S.N,
        MCI.pflowK_cov(:,:,n)=cov(squeeze(sw(assign.w_flow,n,post_ind))');
    end
    
    % Actual Samples (group level)
    MCI.pflow.P=sm(assign.w_flow,:);
    MCI.pflow.ind=post_ind;
    MCI.pflow.Ep=mean(sm(assign.w_flow,post_ind),2);
else
    MCI.pflowK=MCI.sv_mean(assign.v_flow);
    MCI.pflowK_cov=cov(v(post_ind,assign.v_flow));
    
    % Actual Samples:
    MCI.pflow.P=v(:,assign.v_flow)';
    MCI.pflow.ind=post_ind;
    MCI.pflow.Ep=MCI.pflowK;
end

MCI.pflow.type='sample';
MCI.pinit.type='sample';