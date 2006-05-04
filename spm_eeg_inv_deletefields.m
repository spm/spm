function D = spm_eeg_inv_deletefields(S,Cflags,Vcheck)

%=======================================================================
% Delete part of the whole of an inverse analysis
%
% FORMAT D = spm_eeg_inv_deletefields(S,Ian,Flock)
% Input:
%
% S		    - input data struct (optional)
% Ian       - vector of indices of the analysis to be deleted
%             (default: display the analysis and ask the user)
% Flock     - flag the enables/disables the deletion of associated files
%             1: unlocked
%             0: locked (default)
%
% Output:
% D			- same data struct including the new files and parameters
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_deletefields.m 507 2006-05-04 05:44:19Z Darren $

def_Ian   = [];
def_Flock = 0;

if nargin == 0
    D     = spm_select(1, '.mat', 'Select EEG/MEG mat file');
	D     = spm_eeg_ldata(D);
    Ian   = def_Ian;
    Flock = def_Flock;
elseif nargin == 1
    if isstruct(S);
        D = S;
        Ian   = def_Ian;
        Flock = def_Flock;
    else
        error(sprintf('Wrong input arguments\n'));
    end
elseif nargin == 2
    D = S;
    Flock = def_Flock;
elseif nargin == 3;
    D = S;
else
    error(sprintf('Wrong input arguments\n'));
end

if ~isfield(D,'inv')
    error(sprintf('No inverse structure for these data\n'));
end

val = length(D.inv);

if isempty(Ian)
    for i = 1:val
        disp(['Inverse analysis ' num2str(i) ': ' D.inv{i}.comment]);
    end
    Ian = spm_input('Analyses to be deleted','+1');
end

if max(Ian) > val
    error(sprintf('Wrong inverse analysis indices\n'));
end

Icomp = [1:val];
Icomp = setdiff(Icomp,Ian);

if Flock == 1   % Deleting files also
    
    tbd_priors = [];
    tbk_priors = [];
    for k = 1:val
        
        % Mesh
        files_sMRI(k)           = D.inv{k}.mesh.sMRI;
        files_nobias(k)         = D.inv{k}.mesh.nobias;
        files_def(k)            = D.inv{k}.mesh.def;
        files_invdef(k)         = D.inv{k}.mesh.invdef;
        files_msk_iskull(k)     = D.inv{k}.mesh.msk_iskull;
        files_msk_scalp(k)      = D.inv{k}.mesh.msk_scalp;
        files_tess_ctx(k)       = D.inv{k}.mesh.tess_ctx;
        files_tess_iskull(k)    = D.inv{k}.mesh.tess_iskull;
        files_tess_scalp(k)     = D.inv{k}.mesh.tess_scalp;
        files_CtxGeoDist(k)     = D.inv{k}.mesh.CtxGeoDist;
        
        % Datareg
        files_sens(k)           = D.inv{k}.datareg.sens;
        files_fid(k)            = D.inv{k}.datareg.fid;
        files_fidmri(k)         = D.inv{k}.datareg.fidmri;
        files_hsp(k)            = D.inv{k}.datareg.hsp;
        files_scalpvert(k)      = D.inv{k}.datareg.scalpvert;
        files_sens_coreg(k)     = D.inv{k}.datareg.sens_coreg;
        files_fid_coreg(k)      = D.inv{k}.datareg.fid_coreg;
        files_hsp_coreg(k)      = D.inv{k}.datareg.hsp_coreg;
        files_eeg2mri(k)        = D.inv{k}.datareg.eeg2mri;
        
        % Forward
        files_bst_options(k)    = D.inv{k}.forward.bst_options;
        files_bst_channel(k)    = D.inv{k}.forward.bst_channel;
        files_bst_tess(k)       = D.inv{k}.forward.bst_tess;
        files_bst_gainmat(k)    = D.inv{k}.forward.bst_gainmat;
        files_gainmat(k)        = D.inv{k}.forward.gainmat;
        files_pcagain(k)        = D.inv{k}.forward.pcagain;
        
        % Inverse
        files_resfile(k)        = D.inv{k}.inverse.resfile;
        if ~isempty(find(Ian == k))
            for j = 1:length(D.inv{k}.inverse.priors.level1)
                if ~isempty(D.inv{k}.inverse.priors.level1{j}) & ~strcmp(D.inv{k}.inverse.priors.level1{j},'none')
                    tbd_priors = [tbd_priors D.inv{k}.inverse.priors.level1{j}.filename];
                end
            end
            for j = 1:length(D.inv{k}.inverse.priors.level2)
                if ~isempty(D.inv{k}.inverse.priors.level2{j}) & ~strcmp(D.inv{k}.inverse.priors.level2{j},'none')
                    tbd_priors = [tbd_priors D.inv{k}.inverse.priors.level2{j}.filename];
                end
            end
        else
            for j = 1:length(D.inv(k).inverse.priors.level1)
                if ~isempty(D.inv{k}.inverse.priors.level1{j}) & ~strcmp(D.inv{k}.inverse.priors.level1{j},'none')
                    tbk_priors = [tbk_priors D.inv{k}.inverse.priors.level1{j}.filename];
                end
            end
            for j = 1:length(D.inv(k).inverse.priors.level2)
                if ~isempty(D.inv{k}.inverse.priors.level2{j}) & ~strcmp(D.inv{k}.inverse.priors.level2{j},'none')
                    tbk_priors = [tbk_priors D.inv{k}.inverse.priors.level2{j}.filename];
                end
            end
        end
        
    end
    
    % Mesh - deletion
    tbd_sMRI           = unique(files_sMRI(Ian));
    tbd_nobias         = unique(files_nobias(Ian));
    tbd_def            = unique(files_def(Ian));
    tbd_invdef         = unique(files_invdef(Ian));
    tbd_msk_iskull     = unique(files_msk_iskull(Ian));
    tbd_msk_scalp      = unique(files_msk_scalp(Ian));
    tbd_tess_ctx       = unique(files_tess_ctx(Ian));
    tbd_tess_iskull    = unique(files_tess_iskull(Ian));
    tbd_tess_scalp     = unique(files_tess_scalp(Ian));
    tbd_CtxGeoDist     = unique(files_CtxGeoDist(Ian));
    
    tbk_sMRI           = unique(files_sMRI(Icomp));
    tbk_nobias         = unique(files_nobias(Icomp));
    tbk_def            = unique(files_def(Icomp));
    tbk_invdef         = unique(files_invdef(Icomp));
    tbk_msk_iskull     = unique(files_msk_iskull(Icomp));
    tbk_msk_scalp      = unique(files_msk_scalp(Icomp));
    tbk_tess_ctx       = unique(files_tess_ctx(Icomp));
    tbk_tess_iskull    = unique(files_tess_iskull(Icomp));
    tbk_tess_scalp     = unique(files_tess_scalp(Icomp));
    tbk_CtxGeoDist     = unique(files_CtxGeoDist(Icomp));
    
    indd_sMRI = setdiff(tbd_sMRI,tbk_sMRI);
    for i = 1:length(indd_sMRI)
        delete(indd_sMRI(i));
    end
    indd_nobias = setdiff(tbd_nobias,tbk_nobias);
    for i = 1:length(indd_nobias)
        delete(indd_nobias(i));
    end
    indd_def = setdiff(tbd_def,tbk_def);
    for i = 1:length(indd_def)
        delete(indd_def(i));
    end
    indd_invdef = setdiff(tbd_invdef,tbk_invdef);
    for i = 1:length(indd_invdef)
        delete(indd_invdef(i));
    end
    indd_msk_iskull = setdiff(tbd_msk_iskull,tbk_msk_iskull);
    for i = 1:length(indd_msk_iskull)
        delete(indd_msk_iskull(i));
    end
    indd_msk_scalp = setdiff(tbd_msk_scalp,tbk_msk_scalp);
    for i = 1:length(indd_msk_scalp)
        delete(indd_msk_scalp(i));
    end
    indd_tess_ctx = setdiff(tbd_tess_ctx,tbk_tess_ctx);
    for i = 1:length(indd_tess_ctx)
        delete(indd_tess_ctx(i));
    end
    indd_tess_iskull = setdiff(tbd_tess_iskull,tbk_tess_iskull);
    for i = 1:length(indd_tess_iskull)
        delete(indd_tess_iskull(i));
    end
    indd_tess_scalp  = setdiff(tbd_tess_scalp,tbk_tess_scalp);
    for i = 1:length(indd_tess_scalp)
        delete(indd_tess_scalp(i));
    end
    indd_CtxGeoDist = setdiff(tbd_CtxGeoDist,tbk_CtxGeoDist);
    for i = 1:length(indd_CtxGeoDist)
        delete(indd_CtxGeoDist(i));
    end
    
    
    % Datareg
    tbd_sens           = unique(files_sens(Ian));
    tbd_fid            = unique(files_fid(Ian));
    tbd_fidmri         = unique(files_fidmri(Ian));
    tbd_hsp            = unique(files_hsp(Ian));
    tbd_scalpvert      = unique(files_scalpvert(Ian));
    tbd_sens_coreg     = unique(files_sens_coreg(Ian));
    tbd_fid_coreg      = unique(files_fid_coreg(Ian));
    tbd_hsp_coreg      = unique(files_hsp_coreg(Ian));
    tbd_eeg2mri        = unique(files_eeg2mri(Ian));
    
    tbk_sens           = unique(files_sens(Icomp));
    tbk_fid            = unique(files_fid(Icomp));
    tbk_fidmri         = unique(files_fidmri(Icomp));
    tbk_hsp            = unique(files_hsp(Icomp));
    tbk_scalpvert      = unique(files_scalpvert(Icomp));
    tbk_sens_coreg     = unique(files_sens_coreg(Icomp));
    tbk_fid_coreg      = unique(files_fid_coreg(Icomp));
    tbk_hsp_coreg      = unique(files_hsp_coreg(Icomp));
    tbk_eeg2mri        = unique(files_eeg2mri(Icomp));
    
    indd_sens = setdiff(tbd_sens,tbk_sens);
    for i = 1:length(indd_sens)
        delete(indd_sens(i));
    end
    indd_fid = setdiff(tbd_fid,tbk_fid);
    for i = 1:length(indd_fid)
        delete(indd_fid(i));
    end
    indd_fidmri = setdiff(tbd_fidmri,tbk_fidmri);
    for i = 1:length(indd_fidmri)
        delete(indd_fidmri(i));
    end
    indd_hsp = setdiff(tbd_hsp,tbk_hsp);
    for i = 1:length(indd_hsp)
        delete(indd_hsp(i));
    end
    indd_scalpvert = setdiff(tbd_scalpvert,tbk_scalpvert);
    for i = 1:length(indd_scalpvert)
        delete(indd_scalpvert(i));
    end
    indd_sens_coreg = setdiff(tbd_sens_coreg,tbk_sens_coreg);
    for i = 1:length(indd_sens_coreg)
        delete(indd_sens_coreg(i));
    end
    indd_fid_coreg = setdiff(tbd_sens,tbk_fid_coreg);
    for i = 1:length(indd_fid_coreg)
        delete(indd_fid_coreg(i));
    end
    indd_hsp_coreg = setdiff(tbd_hsp_coreg,tbk_hsp_coreg);
    for i = 1:length(indd_hsp_coreg)
        delete(indd_hsp_coreg(i));
    end
    indd_eeg2mri = setdiff(tbd_eeg2mri,tbk_eeg2mri);
    for i = 1:length(indd_eeg2mri)
        delete(indd_eeg2mri(i));
    end

    
    % Forward
    tbd_bst_options    = unique(files_bst_options(Ian));
    tbd_bst_channel    = unique(files_bst_channel(Ian));
    tbd_bst_tess       = unique(files_bst_tess(Ian));
    tbd_bst_gainmat    = unique(files_bst_gainmat(Ian));
    tbd_gainmat        = unique(files_gainmat(Ian));
    tbd_pcagain        = unique(files_pcagain(Ian));
    
    tbk_bst_options    = unique(files_bst_options(Icomp));
    tbk_bst_channel    = unique(files_bst_channel(Icomp));
    tbk_bst_tess       = unique(files_bst_tess(Icomp));
    tbk_bst_gainmat    = unique(files_bst_gainmat(Icomp));
    tbk_gainmat        = unique(files_gainmat(Icomp));
    tbk_pcagain        = unique(files_pcagain(Icomp));
    
    indd_bst_options   = setdiff(tbd_bst_options,tbk_bst_options);
    for i = 1:length(indd_bst_options)
        delete(indd_bst_options(i));
    end
    indd_bst_channel   = setdiff(tbd_bst_channel,tbk_bst_channel);
    for i = 1:length(indd_bst_channel)
        delete(indd_bst_channel(i));
    end
    indd_bst_tess      = setdiff(tbd_bst_tess,tbk_bst_tess);
    for i = 1:length(indd_bst_tess)
        delete(indd_bst_tess(i));
    end
    indd_bst_gainmat   = setdiff(tbd_bst_gainmat,tbk_bst_gainmat);
    for i = 1:length(indd_bst_gainmat)
        delete(indd_bst_gainmat(i));
    end
    indd_gainmat       = setdiff(tbd_gainmat,tbk_gainmat);
    for i = 1:length(indd_gainmat)
        delete(indd_gainmat(i));
    end
    indd_pcagain       = setdiff(tbd_pcagain,tbk_pcagain);
    for i = 1:length(indd_pcagain)
        delete(indd_pcagain(i));
    end

    
    % Inverse
    tbd_resfile        = unique(files_resfile(Ian));
    tbd_priors         = unique(tbd_priors);

    tbk_resfile        = unique(files_resfile(Icomp));
    tbk_priors         = unique(tbk_priors);
    
    indd_resfile       = setdiff(tbd_resfile,tbk_resfile);
    for i = 1:length(indd_resfile)
        delete(indd_resfile(i));
    end
    indd_priors        = setdiff(tbd_priors,tbk_priors);
    for i = 1:length(indd_priors)
        delete(indd_priors(i));
    end

end


% Deleting analyses
S = D.inv(Icomp);
D.inv = S;


if spm_matlab_version_chk('7.1') >= 0
	save(fullfile(D.path, D.fname), '-V6', 'D');
else
	save(fullfile(D.path, D.fname), 'D');
end
