function mfx = spm_cfg_mfx
% SPM Configuration file for MFX
%______________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_cfg_mfx.m 4023 2010-07-28 18:41:36Z guillaume $

% ---------------------------------------------------------------------
% dir Directory
% ---------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.help    = {'Select a directory where the SPM.mat file containing the specified design matrix will be written.'};
dir.filter = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];

% ---------------------------------------------------------------------
% spmmat Select SPM.mat
% ---------------------------------------------------------------------
spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
spmmat.name    = 'Select SPM.mat files';
spmmat.help    = {...
    'Select the SPM.mat files that contains first-level designs.'
    'They must have the same number of parameters for each session.'
    ['These are assumed to represent session-specific realisations of ' ...
    '2nd-level effects.']};
spmmat.filter  = 'mat';
spmmat.ufilter = '^SPM\.mat$';
spmmat.num     = [1 Inf];

% ---------------------------------------------------------------------
% ffx Create first-level design
% ---------------------------------------------------------------------
ffx       = cfg_exbranch;
ffx.tag   = 'ffx';
ffx.name  = 'Create repeated-measure design';
ffx.val   = {dir spmmat};
ffx.help  = {'Create repeated-measure multi-session first-level design'};
ffx.prog  = @spm_local_ffx;
ffx.vout  = @vout_ffx;

% ---------------------------------------------------------------------
% spmmat Select SPM.mat
% ---------------------------------------------------------------------
spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
spmmat.name    = 'Select SPM.mat';
spmmat.help    = {...
    'Design and estimation structure after a 1st-level analysis'};
spmmat.filter  = 'mat';
spmmat.ufilter = '^SPM\.mat$';
spmmat.num     = [1 1];

% ---------------------------------------------------------------------
% spec MFX Specification
% ---------------------------------------------------------------------
spec       = cfg_exbranch;
spec.tag   = 'spec';
spec.name  = 'MFX Specification';
spec.val   = {spmmat};
spec.help  = {'MFX Specification'};
spec.prog  = @spm_local_mfx;
spec.vout  = @vout_mfx;

% ---------------------------------------------------------------------
% mfx Mixed-effects (MFX) analysis
% ---------------------------------------------------------------------
mfx         = cfg_choice;
mfx.tag     = 'mfx';
mfx.name    = 'Mixed-effects (MFX) analysis';
mfx.help    = {'Mixed-effects (MFX) analysis'};
mfx.values  = {ffx spec};

% =====================================================================
function out = spm_local_ffx(job)
spmmat = job.spmmat;
SPMS = cell(size(spmmat));
for i=1:numel(spmmat)
    load(spmmat{i},'SPM');
    SPMS{i} = SPM;
end

matlabbatch{1}.spm.stats.fmri_spec.dir            = cellstr(job.dir);
matlabbatch{1}.spm.stats.fmri_spec.timing.units   = SPMS{1}.xBF.UNITS;
matlabbatch{1}.spm.stats.fmri_spec.timing.RT      = SPMS{1}.xY.RT;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = SPMS{1}.xBF.T;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = SPMS{1}.xBF.T0;
switch SPMS{1}.xBF.name
    case 'hrf'
        matlabbatch{1}.spm.stats.fmri_spec.hrf.derivs = [0 0];
    case 'hrf (with time derivative)'
        matlabbatch{1}.spm.stats.fmri_spec.hrf.derivs = [1 0];
    case 'hrf (with time and dispersion derivatives)'
        matlabbatch{1}.spm.stats.fmri_spec.hrf.derivs = [1 1];
    case 'Fourier set'
        matlabbatch{1}.spm.stats.fmri_spec.bases.fourier.length = SPMS{1}.xBF.length;
        matlabbatch{1}.spm.stats.fmri_spec.bases.fourier.order  = SPMS{1}.xBF.order;
    case 'Fourier set (Hanning)'
        matlabbatch{1}.spm.stats.fmri_spec.bases.fourier_han.length = SPMS{1}.xBF.length;
        matlabbatch{1}.spm.stats.fmri_spec.bases.fourier_han.order  = SPMS{1}.xBF.order;
    case 'Gamma functions'
        matlabbatch{1}.spm.stats.fmri_spec.bases.gamma.length = SPMS{1}.xBF.length;
        matlabbatch{1}.spm.stats.fmri_spec.bases.gamma.order  = SPMS{1}.xBF.order;
    case 'Finite Impulse Response'
        matlabbatch{1}.spm.stats.fmri_spec.bases.fir.length = SPMS{1}.xBF.length;
        matlabbatch{1}.spm.stats.fmri_spec.bases.fir.order  = SPMS{1}.xBF.order;
end
matlabbatch{1}.spm.stats.fmri_spec.volt   = SPMS{1}.xBF.Volterra;
matlabbatch{1}.spm.stats.fmri_spec.global = SPMS{1}.xGX.iGXcalc;
if ~isempty(SPMS{1}.xM.VM)
    matlabbatch{1}.spm.stats.fmri_spec.mask = cellstr(SPMS{1}.xM.VM.fname); % can be intersection
end
if strncmp('AR',SPMS{1}.xVi.form,2)
    matlabbatch{1}.spm.stats.fmri_spec.cvi  = 'AR(1)';
else
    matlabbatch{1}.spm.stats.fmri_spec.cvi  = 'none';
end

k = 1;
for i=1:numel(SPMS)
    n = cumsum([1 SPMS{i}.nscan]);
    for j=1:numel(SPMS{i}.Sess)
        matlabbatch{1}.spm.stats.fmri_spec.sess(k).scans = cellstr(SPMS{i}.xY.P(n(j):n(j+1)-1,:));
        for l=1:numel(SPMS{i}.Sess(j).U)
            matlabbatch{1}.spm.stats.fmri_spec.sess(k).cond(l).name = SPMS{i}.Sess(j).U(l).name{1};
            matlabbatch{1}.spm.stats.fmri_spec.sess(k).cond(l).onset = SPMS{i}.Sess(j).U(l).ons;
            matlabbatch{1}.spm.stats.fmri_spec.sess(k).cond(l).duration = SPMS{i}.Sess(j).U(l).dur;
            o = 1;
            for m=1:numel(SPMS{i}.Sess(j).U(l).P)
                switch SPMS{i}.Sess(j).U(l).P(m).name
                    case 'time'
                        matlabbatch{1}.spm.stats.fmri_spec.sess(k).cond(l).tmod = SPMS{i}.Sess(j).U(l).P(m).h;
                    case 'none'
                    otherwise
                        matlabbatch{1}.spm.stats.fmri_spec.sess(k).cond(l).pmod(o).name  = SPMS{i}.Sess(j).U(l).P(m).name;
                        matlabbatch{1}.spm.stats.fmri_spec.sess(k).cond(l).pmod(o).param = SPMS{i}.Sess(j).U(l).P(m).P;
                        matlabbatch{1}.spm.stats.fmri_spec.sess(k).cond(l).pmod(o).poly  = SPMS{i}.Sess(j).U(l).P(m).h;
                        o = o + 1;
                end
            end
        end
        for l=1:numel(SPMS{i}.Sess(j).C.name)
            matlabbatch{1}.spm.stats.fmri_spec.sess(k).regress(l).name = SPMS{i}.Sess(j).C.name{l};
            matlabbatch{1}.spm.stats.fmri_spec.sess(k).regress(l).val  = SPMS{i}.Sess(j).C.C(:,l);
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(k).hpf = SPMS{i}.xX.K(j).HParam;
        k = k + 1;
    end
end
spm_jobman('run',matlabbatch);

out.spmmat{1} = fullfile(job.dir{1},'SPM.mat');

% =====================================================================
function dep = vout_ffx(job)
dep = cfg_dep;
dep.sname      = 'SPM.mat File';
dep.src_output = substruct('.','spmmat');
dep.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});

% =====================================================================
function out = spm_local_mfx(job)
load(job.spmmat{1},'SPM');
spm_mfx(SPM);
out.spmmat{1} = fullfile(fileparts(job.spmmat{1}),'mfx','SPM.mat');

% =====================================================================
function dep = vout_mfx(job)
dep = cfg_dep;
dep.sname      = 'SPM.mat File';
dep.src_output = substruct('.','spmmat');
dep.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
