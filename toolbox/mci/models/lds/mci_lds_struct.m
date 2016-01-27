function [M,U,names] = mci_lds_struct (M)
% LDS constrained: Initialise model structure
% FORMAT [M,U,names] = mci_lds_struct (M)
%
% M.d       Number of regions
% M.sd      Observation noise SD
% M.name    'uncoupled','forward','backward','bidirectional'
% M.R       Initial state
% M.t       Vector of Times
% M.drop    final value as proportion of initial value
%           eg. 0.5 indicates typical state at M.t(end) is
%           half of M.t(1). Used to set M.a_typical, typical
%           self connection values
%
% M         Model structure
% U         Inputs
% names     Names of variables
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_lds_struct.m 6697 2016-01-27 14:57:28Z spm $

d=M.d;
try, R=M.R; catch, R=ones(d,1); end

M.N=length(M.t);
M.T=M.t(end);
M.dt=M.t(2)-M.t(1);

switch M.name
    case 'uncoupled',
        disp('No coupling among regions');
        Aconn=[];
    case 'forward',
        for i=1:d-1,
            Aconn(i,:)=[i+1,i];
        end
    case 'backward',
        for i=1:d-1,
            Aconn(i,:)=[i,i+1];
        end
    case 'bidirectional',
        % Forward
        for i=1:d-1,
            Aconn(i,:)=[i+1,i];
        end
        % Backward
        for i=1:d-1,
            Aconn(d+i-1,:)=[i,i+1];
        end
    otherwise
        disp('Unknown Hypothesis');
end

% Number of connections between regions
M.Nb=size(Aconn,1);
M.Aconn=Aconn;

M.f='mci_lds_fx';  
M.dfdx='mci_lds_dfdx';
M.g='mci_lds_gx';

% Initial states
M.x0=R;
M.x=zeros(d,1);
M.m=0;
M.n=d;
M.l=d;

% Typical diagonal entry in A matrix
if ~isfield(M,'drop')
    M.drop=0.5;
end
M.a_typical=log(M.drop)/M.T;
    
% Self connections
M.pEs.self=zeros(d,1);

% Strength of priors
%M.sd_self=1/4;
%M.sd_between=exp(-5);
M.sd_self=1;
M.sd_between=0.05;

% Coupling among regions
if M.Nb>0
    M.pEs.between=zeros(M.Nb,1);
    M.pC=diag([M.sd_self^2*ones(d,1);M.sd_between^2*ones(M.Nb,1)]);
else
    M.pC=M.sd_self^2*eye(d);
end

M.pE=spm_vec(M.pEs);
M.Npflow=length(M.pE);
M.Npout=0;
M.ipC=inv(M.pC);

%M.int='euler';
%M.int='ode15';
M.int='sundials';

% Inputs
U=zeros(1,M.N);
M.m=1;

% Likelihood function
M.L='spm_mci_glike';
M.Ce=M.sd^2*eye(d);

if strcmp(M.name,'forward')
    
    for j=1:d,
        jn=int2str(j);
        names{j}=['a_{',jn,jn,'}'];
    end
    for i=1:d-1,
        jn=int2str(i);
        j1n=int2str(i+1);
        names{i+j}=['a_{',j1n,jn,'}'];
    end
    
end

