function meeg = spm_cfg_dcm_meeg
% SPM Configuration file for DCM for M/EEG
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_cfg_dcm_meeg.m 6000 2014-05-20 17:16:38Z guillaume $

% -------------------------------------------------------------------------
% dcmmat Select DCM_*.mat
% -------------------------------------------------------------------------
dcmmat         = cfg_files;
dcmmat.tag     = 'dcmmat';
dcmmat.name    = 'Select DCM_*.mat';
dcmmat.help    = {'Select DCM_*.mat files.'};
dcmmat.filter  = 'mat';
dcmmat.ufilter = '^DCM_.*\.mat$';
dcmmat.num     = [1 Inf];

% -------------------------------------------------------------------------
% estimate Estimate
% -------------------------------------------------------------------------
estimate      = cfg_exbranch;
estimate.tag  = 'estimate';
estimate.name = 'Estimate';
estimate.val  = { dcmmat };
estimate.help = {'Estimate parameters of a DCM.'};
estimate.prog = @spm_run_dcm_meeg_est;
estimate.vout = @vout_dcm_meeg;

% ---------------------------------------------------------------------
% meeg Dynamic Causal Model for M/EEG
% ---------------------------------------------------------------------
meeg         = cfg_choice; 
meeg.tag     = 'meeg';
meeg.name    = 'DCM for M/EEG';
meeg.help    = {'Dynamic Causal Modelling for M/EEG'};
meeg.values  = { estimate };

%==========================================================================
function out = spm_run_dcm_meeg_est(job)
%==========================================================================
for i=1:numel(job.dcmmat)
    spm_dcm_erp(job.dcmmat{i});
end
out = job.dcmmat;

%==========================================================================
function dep = vout_dcm_meeg(varargin)
%==========================================================================
dep(1)            = cfg_dep;
dep(1).sname      = 'DCM mat File(s)';
dep(1).src_output = substruct('.','dcmmat');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});