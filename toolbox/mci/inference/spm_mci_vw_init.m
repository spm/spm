function [w_init,v_init,assign] = spm_mci_vw_init (MCI)
% Initialise fixed and random effects
% FORMAT [w_init,v_init,assign] = spm_mci_vw_init (MCI)
%
% MCI       MCI data structure
%
% w_init    initial rfx values
% v_init    initial ffx values
% assign    data structure describing how rfx/ffx are assigned 
%           to initial conditions and flow params
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_vw_init.m 6275 2014-12-01 08:41:18Z will $

assign=MCI.assign;
fixed=MCI.fixed;
S=MCI.S;

Np=length(spm_vec(MCI.M{1}.pE));
if isfield(MCI.M{1},'x0')
    % For dynamical systems
    d=size(MCI.M{1}.x0,1);
end

% First part of w/v vector is for initial condition parameters
% Second part for flow parameters
w_init=[];v_init=[];
switch assign.init_par,
    case 'random',
        assign.w_init=[1:d];
        assign.v_init=[];
        if isfield(MCI,'pinit0');
            w_init=MCI.pinit0;
        end
    case 'fixed',
        assign.v_init=[1:d];
        assign.w_init=[];
        if isfield(MCI,'pinit0');
            v_init=MCI.pinit0;
        end
    otherwise
        % Assume 'known'
        assign.v_init=[];
        assign.w_init=[];
end

if strcmp(assign.flow_par,'random')
    assign.w_flow=length(assign.w_init)+[1:Np];
    assign.v_flow=0;
    if isfield(MCI,'pflow0');
        w_init=[w_init;MCI.pflow0];
    end
else
    assign.v_flow=length(assign.v_init)+[1:Np];
    if isfield(MCI,'pflow0');
        v_init=[v_init;MCI.pflow0];
    end
end

% Initialise fixed effects
if isempty(v_init)
    v_init=spm_normrnd(fixed.vpE,fixed.pC,1);
end

% Initialise random effects
if isempty(w_init)
    m=S.prior.m;
    C=S.prior.B/S.prior.a;
    w_init=spm_normrnd(m,C,S.N);
end