function [w_init,v_init,assign,update_ffx,update_rfx] = spm_mci_vw_init (MCI)
% Initialise fixed and random effects
% FORMAT [w_init,v_init,assign,update_ffx,update_rfx] = spm_mci_vw_init (MCI)
%
% MCI       MCI data structure
%
% w_init        initial rfx values
% v_init        initial ffx values
% assign        data structure describing how rfx/ffx are assigned 
%               to initial conditions, flow and output params
% update_ffx    (1/0)
% update_rfx    (1/0)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_vw_init.m 6697 2016-01-27 14:57:28Z spm $

assign=MCI.assign;
try, fixed=MCI.fixed; catch, fixed=[]; end
S=MCI.S;

Np=length(spm_vec(MCI.M{1}.pE));
Nflow=MCI.M{1}.Npflow;
Nout=MCI.M{1}.Npout;

if isfield(MCI.M{1},'x0')
    % For dynamical systems
    d=size(MCI.M{1}.x0,1);
end

% First part of w/v vector is for initial condition parameters
% Second part for flow parameters
% Third part for output parameters
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

switch assign.flow_par,
    case 'random',
        assign.w_flow=length(assign.w_init)+[1:Nflow];
        assign.v_flow=[];
        if isfield(MCI,'pflow0');
            w_init=[w_init;MCI.pflow0];
        end
    case 'fixed',
        assign.v_flow=length(assign.v_init)+[1:Nflow];
        if isfield(MCI,'pflow0');
            v_init=[v_init;MCI.pflow0];
        end
        assign.w_flow=[];
    otherwise
        % Assume 'known'
        assign.v_flow=[];
        assign.w_flow=[];
end

switch assign.out_par,
    case 'random',
        assign.w_out=length(assign.w_init)+length(assign.w_flow)+[1:Nout];
        assign.v_out=[];
        if isfield(MCI,'pout0');
            w_init=[w_init;MCI.pout0];
        end
    case 'fixed',
        assign.v_out=length(assign.v_init)+length(assign.v_flow)+[1:Nout];
        if isfield(MCI,'pout0');
            v_init=[v_init;MCI.pout0];
        end
        assign.w_out=[];
    otherwise,
        % Assume 'known'
        assign.v_out=[];
        assign.w_out=[];
end

if ~isempty(fixed)
    % Initialise fixed effects
    if isempty(v_init)
        v_init=spm_normrnd(fixed.vpE,fixed.pC,1);
    end
end

% Initialise random effects
if isempty(w_init)
    m=S.prior.m;
    C=S.prior.B/S.prior.a;
    w_init=spm_normrnd(m,C,S.N);
end


% Update RFX/FFX
update_ffx=0;
update_rfx=0;
if strcmp(assign.init_par,'fixed')
    update_ffx=1;
end
if strcmp(assign.init_par,'random')
    update_rfx=1;
end
if strcmp(assign.flow_par,'fixed')
    update_ffx=1;
end
if strcmp(assign.flow_par,'random')
    update_rfx=1;
end
if strcmp(assign.out_par,'fixed')
    update_ffx=1;
end
if strcmp(assign.out_par,'random')
    update_rfx=1;
end