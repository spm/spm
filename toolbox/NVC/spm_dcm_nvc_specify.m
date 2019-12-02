function DCM = spm_dcm_nvc_specify(SPM,xY_fMRI,MEEG,model,n_exclude,fmri_cond,options)
% Specify unestimated structure for (multimodal) DCM for fMRI and M/EEG
% FORMAT DCM = spm_dcm_nvc_specify(SPM,xY_fMRI, MEEG, Model,N_exclude,fmri_cond,options)
%
% See spm_dcm_nvc.m for detailed descriptions of the parameters
%
% Inputs:
% -------------------------------------------------------------------------
% SPM          -  SPM structure or location of SPM.mat
% xY_fMRI      -  Cell array of VOI filenames (the same order as sources in EEG DCM)
% MEEG         -  Location of DCM for M/EEG .mat file or DCM structure
% model        -  Model space definition (see spm_dcm_nvc.m)
% n_exclude    -  Which neuronal populations should drive haemodynamics (optional)
% fmri_cond    -  Which fMRI conditions to include (optional)
% options      -  DCM options
%
% Evaluates:
% -------------------------------------------------------------------------
%
% DCM          -  unestimated DCM
%__________________________________________________________________________
% Jafarian, A., Litvak, V., Cagnan, H., Friston, K.J. and Zeidman, P., 2019.
% Neurovascular coupling: insights from multi-modal dynamic causal modelling
% of fMRI and MEG. arXiv preprint arXiv:1903.07478.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Amirhossein Jafarian
% $Id: spm_dcm_nvc_specify.m 7735 2019-12-02 09:15:27Z peter $

%-Get SPM file and DCM for MEG
%--------------------------------------------------------------------------
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

% Prepare fMRI VOIs
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

if isvector(fmri_cond)
    u = find(fmri_cond == 1);
else
    u = 1:length(Sess.U);
end

for i = u
    for j = 1:length(Sess.U(i).name)
        U.u             = [U.u Sess.U(i).u(33:end,j)];
        U.name{end + 1} = Sess.U(i).name{j};
        U.idx           = [U.idx; i j];
    end
end

%-VOI timings
%--------------------------------------------------------------------------
RT     = SPM.xY.RT;
t0     = spm_get_defaults('stats.fmri.t0');
t      = spm_get_defaults('stats.fmri.t');
T0     = RT * t0 / t;
DCM.delays = repmat(T0,1,m);

%-Echo time (TE) of data acquisition
%--------------------------------------------------------------------------
TE     = options.TE;

% Prepare fMRI data
%--------------------------------------------------------------------------
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

%-Store all variables in DCM structure
%--------------------------------------------------------------------------
DCM.U                   =  U;
DCM.Y                   =  Y;
DCM.xY                  =  xY;
DCM.v                   =  v;
DCM.n                   =  n;
DCM.TE                  =  TE;
DCM.options             =  options;
DCM.options.nmm         =  'TFM';
DCM.model               =  model;
DCM.N                   =  n_exclude;
DCM.MEEG                =  EEG_DCM.DCM ;
