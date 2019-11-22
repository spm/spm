function  [DCM]= spm_dcm_mm_specify(SPM,xY_fMRI, MEEG, Model,N_exclude,Sess_exclude,options)
% Specify un-estimated structure for (multimodal) DCM for fMRI and M/EEG
% FORMAT DCM = spm_dcm_specify(SPM,xY_fMRI, MEEG, Model, N_exclude, Sess_exclude)
%
%Input
%--------------------------------------------------------------------------
% SPM          -  Address of SPM or SPM.mat
% xY_fMRI      -  Address of VOIs (in the same order as neuronal sources in DCM for M/EEG)
% MEEG         -  Address of DCM or DCM_???.mat
% Model        -  Model space definition
% N_exclude    -  Excluding any neuronal drives to be used for DCM for fMRI (optional)
% Sess_exclude -  Excluding any sessions from SPM.mat (optional)
%Output
%--------------------------------------------------------------------------

% DCM          -  unestiamted DCM
%__________________________________________________________________________
% Jafarian, A., Litvak, V., Cagnan, H., Friston, K.J. and Zeidman, P., 2019.
% Neurovascular coupling: insights from multi-modal dynamic causal modelling
% of fMRI and MEG. arXiv preprint arXiv:1903.07478.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Amirhossein Jafarian
% $Id $

%-Get SPM file and DCM for MEG
%--------------------------------------------------------------------------
if nargin <=4
    N_exclude     = ones (1,4)   ;
    Sess_exclude  = 'not defined';
end
if nargin <=5
    Sess_exclude  = 'not defined';
end

if ischar(SPM)
    swd = spm_file(SPM);
    try
        load(fullfile(swd,'SPM.mat'))
    catch
        error('Cannot read %s.',fullfile(swd,'SPM.mat'));
    end
    SPM.swd = swd;
else
    SPM = SPM.SPM;
    SPM.swd = pwd;
end

if ischar(MEEG)
    try
        EEG_DCM =load(MEEG); 
    catch
        error('Cannot read DCM for M/EEG');
    end
else
    EEG_DCM    = MEEG;
end
%------------------------------------------------------------------------
P     = xY_fMRI ;
m     = numel(P);
xY    = [];
if m  == 0
    error('Cannot read VOIs');
end
for i = 1:m
    p  = load(P{i},'xY');
    xY = spm_cat_struct(xY,p.xY);
end
% Inputs
%==========================================================================
% Experimental fMRI inputs U
%--------------------------------------------------------------------------
Sess   = SPM.Sess(xY(1).Sess);
U.dt   = Sess.U(1).dt;
U.name = {};
U.u    = [];
U.idx  = [];

if     ~ ischar(Sess_exclude)
    u      = find (Sess_exclude == 1);
else
    u      = 1: length(Sess.U);
end

for i = u
    for j = 1:length(Sess.U(i).name)
        U.u             = [U.u Sess.U(i).u(33:end,j)];
        U.name{end + 1} = Sess.U(i).name{j};
        U.idx           = [U.idx; i j];
    end
end
% Timings
%==========================================================================
%-VOI timings
%--------------------------------------------------------------------------
RT     = SPM.xY.RT;
t0     = spm_get_defaults('stats.fmri.t0');
t      = spm_get_defaults('stats.fmri.t');
T0     = RT * t0 / t;
DCM.delays = repmat(SPM.xY.RT/2,m,1);

%-Echo time (TE) of data acquisition
%--------------------------------------------------------------------------
TE    = 0.04;
%==========================================================================
% Response
%==========================================================================
n     = length(xY);                      % number of regions
v     = length(xY(1).u);                 % number of time points
Y.dt  = SPM.xY.RT;
Y.X0  = xY(1).X0;
for i = 1:n
    Y.y(:,i)  = xY(i).u;
    Y.name{i} = xY(i).name;
end

%-Error precision components (one for each region) - i.i.d. (because of W)
%--------------------------------------------------------------------------
Y.Q        = spm_Ce(ones(1,n)*v);
%==========================================================================
% DCM structure
%==========================================================================
%-Store all variables in DCM structure
%--------------------------------------------------------------------------
DCM.U                   =  U;
DCM.Y                   =  Y;
DCM.xY                  =  xY;
DCM.v                   =  v;
DCM.n                   =  n;
DCM.TE                  =  TE;
DCM.options.nmm         = 'TFM';
DCM.options.centre      =  options.centre;   
DCM.options.hE          =  options.hE;        
DCM.options.hC          =  options.hC;      
DCM.options.maxit       =  options.maxit;      
DCM.model               =  Model;
DCM.N                   =  N_exclude;
DCM.MEEG                =  EEG_DCM.DCM ;
end


