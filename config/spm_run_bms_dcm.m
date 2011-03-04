function out = spm_run_bms_dcm (varargin)
% API to compare DCMs on the basis of their log-evidences. Four methods
% are available to identify the best among alternative models:
%
%  (1) single subject BMS using Bayes factors
%     (see Penny et al, NeuroImage, 2004)
%  (2) fixed effects group BMS using group Bayes factors
%     (see Stephan et al, NeuroImage, 2007)
%  (3) random effects group BMS using exceedance probabilities
%     (see Stephan et al, NeuroImage, 2009)
%  (4) comparing model families
%     (see Penny et al, PLOS-CB, submitted)
%
% Note: All functions use the negative free energy (F) as an approximation
% to the log model evidence.
% _________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% CC Chen & Maria Joao Rosa
% $Id: spm_run_bms_dcm.m 4229 2011-03-04 15:01:34Z maria $

% input
% -------------------------------------------------------------------------
job     = varargin{1};

fname   = 'BMS.mat';                     % Output filename
fname   = fullfile(job.dir{1},fname);    % Output filename (inc. directory)
ld_f    = ~isempty(job.load_f{1});       % Log-evidence file
ld_msp  = ~isempty(job.model_sp{1});     % Model space file
bma_do  = isfield(job.bma,'bma_yes');    % Compute BMA
data_se = ~isempty(job.sess_dcm);        % DCM files


% method
% -------------------------------------------------------------------------
if strcmp(job.method,'FFX');
    method = 'FFX';                      % Fixed-effects
else
    method = 'RFX';                      % Random-effects
end

% check DCM.mat files and BMA
% -------------------------------------------------------------------------
do_bma_famwin = 0;
do_bma_all    = 0;

if bma_do
    if ld_msp
        disp('Loading model space')
        load(job.model_sp{1});
        load(subj(1).sess(1).model(1).fname)
    else
        if data_se
            load(job.sess_dcm{1}(1).mod_dcm{1})
        else
            spm('alert','Please specify DCM.mat files or model_space.mat to do BMA!','Error');
            return
        end
    end
    
    nsamp       = 1e4;
    odds_ratio  = 1/20;
    
    if isfield(job.bma.bma_yes,'bma_famwin')
        do_bma_famwin  = 1;
    else
        if isfield(job.bma.bma_yes,'bma_all')
            do_bma_all = 1;
        else
            bma_fam    = double(job.bma.bma_yes.bma_part);
        end
    end
else
    if ld_msp
        disp('Loading model space')
        load(job.model_sp{1});
        if ~exist('subj','var')
            spm('alert','Incorrect model space file! File must contain ''subj'' structure.','Error');
            return
        end
    end
end

% initialise
% -------------------------------------------------------------------------
F      = [];                % Free energy
N      = {};                % Models

% prepare the data
% -------------------------------------------------------------------------
if  ld_f
    data      = job.load_f{1};
    load(data);
    nm        = size(F,2);                                % No of Models
    ns        = size(F,1);                                % No of Models
    N         = 1:nm;
    subj      = [];                                       % subj structure
    f_fname   = data;                                     % LogE file name
    fname_msp = [];                                       % Model space
    
else
    
    f_fname       = [];                                   % LogE file name
    
    % Verify dimensions of data and model space
    if ld_msp
        ns        = length(subj);                         % No of Subjects
        nsess     = length(subj(1).sess);                 % No of sessions
        nm        = length(subj(1).sess(1).model);        % No of Models
        fname_msp = [job.model_sp{1}];
    else
        ns        = size(job.sess_dcm,2);                 % No of Subjects
        nsess     = size(job.sess_dcm{1},2);              % No of sessions
        nm        = size(job.sess_dcm{1}(1).mod_dcm,1);   % No of Models
        fname_msp = fullfile(job.dir{1},'model_space.mat');
    end
    
    F       = zeros(ns,nm);
    N       = cell(nm);
    
    % Check if No of models > 2
    if nm < 2
        spm('alert','Please select more than one file','Error');
        return
    end
    
    for k=1:ns
        
        % Verify dimensions of data and model space for current session
        if ld_msp
            nsess_now       = length(subj(k).sess);
            nmodels         = length(subj(k).sess(1).model);
        else
            disp(sprintf('Loading DCMs for subject %d', k));
            nsess_now       = size(job.sess_dcm{k},2);
            nmodels         = size(job.sess_dcm{k}(1).mod_dcm,1);
        end
        
        if (nsess_now == nsess && nmodels== nm) % Check no of sess/mods
            
            ID = zeros(nm,nsess);
            
            for j=1:nm
                
                F_sess      = [];
                
                for h = 1:nsess_now
                    
                    if ~ld_msp
                        
                        clear DCM
                        
                        % Load DCM (model)
                        tmp = job.sess_dcm{k}(h).mod_dcm{j};
                        if ~job.verify_id
                            warning off MATLAB:load:variableNotFound
                            DCM.DCM = load(tmp, 'F', 'Ep', 'Cp');
                            warning on MATLAB:load:variableNotFound
                            if ~isfield(DCM.DCM, 'F')
                                DCM     = load(tmp);
                            end
                        else
                            DCM     = load(tmp);
                        end
                        
                        % Free energy for sessions
                        F_sess  = [F_sess,DCM.DCM.F];
                        
                        % Create model space
                        subj(k).sess(h).model(j).fname = tmp;
                        subj(k).sess(h).model(j).F     = DCM.DCM.F;
                        subj(k).sess(h).model(j).Ep    = DCM.DCM.Ep;
                        subj(k).sess(h).model(j).Cp    = DCM.DCM.Cp;
                        
                    else
                        
                        F_sess = [F_sess,subj(k).sess(h).model(j).F];
                        
                    end
                    
                    % Data ID verification. At least for now we'll
                    % re-compute the IDs rather than use the ones stored
                    % with the DCM.
                    if job.verify_id
                        
                        disp(sprintf('Identifying data: model %d', j));
                        
                        M   = DCM.DCM.M;
                        
                        if isfield(DCM.DCM, 'xY')
                            Y = DCM.DCM.xY;  % not fMRI
                        else
                            Y = DCM.DCM.Y;   % fMRI
                        end
                        
                        if isfield(M,'FS')
                            try
                                ID(j,h)  = spm_data_id(feval(M.FS,Y.y,M));
                            catch
                                ID(j,h)  = spm_data_id(feval(M.FS,Y.y));
                            end
                        else
                            ID(j,h) = spm_data_id(Y.y);
                        end
                        
                    end
                end
                
                % Sum over sessions
                F_mod       = sum(F_sess);
                F(k,j)      = F_mod;
                N{j}        = sprintf('model%d',j);
                
            end
            
            if job.verify_id
                failind = find(max(abs(diff(ID)),[],1) > eps);
                if ~isempty(failind)
                    out.files{1} = [];
                    str_error = ['The models for subject ' num2str(k) ...
                        ' session(s) ' num2str(failind) ' were not fitted to the same data.'];
                    spm('alert',str_error,'Error');
                    return
                end
            end
        else
            out.files{1} = [];
            spm('alert','The number of sessions/models should be the same for all subjects!','Error');
            return         
        end
        
    end
end

% bayesian model selection
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% free energy
sumF = sum(F,1);

% family or model level
% -------------------------------------------------------------------------
if isfield(job.family_level,'family_file')
    
    if ~isempty(job.family_level.family_file{1})
        
        load(job.family_level.family_file{1});
        do_family    = 1;
        family.infer = method;
        
        nfam    = size(family.names,2);
        npart   = length(unique(family.partition));
        maxpart = max(family.partition);
        m_indx  = 1:nm;
        
        if nfam ~= npart || npart == 1 || maxpart > npart
            spm('alert','Invalid family file!','Error');
            out.files{1} = [];
            return
        end
    else
        do_family = 0;
    end
    
else
    
    if isempty(job.family_level.family)
        do_family    = 0;
    else
        do_family    = 1;
        nfam         = size(job.family_level.family,2);
        
        names_fam  = {};
        models_fam = [];
        m_indx     = [];
        for f=1:nfam
            names_fam    = [names_fam,job.family_level.family(f).family_name];
            m_indx       = [m_indx,job.family_level.family(f).family_models'];
            models_fam(job.family_level.family(f).family_models) = f;
        end
        
        family.names     = names_fam;
        family.partition = models_fam;
        
        npart   = length(unique(family.partition));
        maxpart = max(family.partition);
        nmodfam = length(m_indx);
        
        if nfam ~= npart || npart == 1 || maxpart > npart || nmodfam > nm
            spm('alert','Invalid family!','Error');
            out.files{1} = [];
            return
        end
        
        family.infer     = method;
        
    end
    
end

% single subject BMS or 1st level( fixed effects) group BMS
% -------------------------------------------------------------------------
if strcmp(method,'FFX');
    
    disp('Computing FFX model/family posteriors ...');
    if ~do_family
        family         = [];
        model.post     = spm_api_bmc(sumF,N,[],[]);
    else
        Ffam           = F(:,sort(m_indx));
        [family,model] = spm_compare_families (Ffam,family);
        P              = spm_api_bmc(sumF,N,[],[],family);
    end
    
    if bma_do,
        if do_family
            if do_bma_famwin
                [fam_max,fam_max_i]  = max(family.post);
                indx  = find(family.partition==fam_max_i);
                post  = model.post(indx);
            else
                if do_bma_all
                    post = model.post;
                    indx = 1:nm;
                else
                    if bma_fam <= nfam && bma_fam > 0 && rem(bma_fam,1) == 0
                        indx  = find(family.partition==bma_fam);
                        post  = model.post(indx);
                    else
                        spm('alert','Incorrect family for BMA!','Error');
                        return
                    end
                end
            end
        else
            post  = model.post;
            indx = 1:nm;
        end
        
        disp('FFX Bayesian model averaging ...');
        
        % bayesian model averaging
        % -----------------------------------------------------------------
        bma = spm_dcm_bma(post,indx,subj,nsamp,odds_ratio);
        
    else
        
        bma = [];
        
    end
    
    if exist(fullfile(job.dir{1},'BMS.mat'),'file')
        load(fname);
        if  isfield(BMS,'DCM') && isfield(BMS.DCM,'ffx')
            spm('alert','Existing BMS.mat file has been over-written!',...
                'Warning');
        end
    end
    BMS.DCM.ffx.data    = fname_msp;
    BMS.DCM.ffx.F_fname = f_fname;
    BMS.DCM.ffx.F       = F;
    BMS.DCM.ffx.SF      = sumF;
    BMS.DCM.ffx.model   = model;
    BMS.DCM.ffx.family  = family;
    BMS.DCM.ffx.bma     = bma;
    
    disp('Saving BMS.mat file...')
    if spm_check_version('matlab','7') >= 0
        save(fname,'-V6','BMS');
    else
        save(fname,'BMS');
    end
    
    out.files{1} = fname;
    
    % 2nd-level (random effects) BMS
    % ---------------------------------------------------------------------
else
    
    disp('Computing RFX model/family posteriors...');
    if ~do_family
        [exp_r,xp,r_samp,g_post] = spm_BMS_gibbs(F);
        model.g_post             = g_post;
        alpha                    = [];
        model.alpha = alpha;
        model.exp_r = exp_r;
        model.xp    = xp;
        family      = [];
    else
        
        Ffam           = F(:,sort(m_indx));
        [family,model] = spm_compare_families(Ffam,family);
        
    end
    
    % display the result
    if do_family
        P = spm_api_bmc(sumF,N,model.exp_r,model.xp,family);
    else
        P = spm_api_bmc(sumF,N,model.exp_r,model.xp);
    end
    
    if bma_do,
        if do_family
            if do_bma_famwin
                [fam_max,fam_max_i]  = max(family.exp_r);
                indx = find(family.partition==fam_max_i);
                post  = model.g_post(:,indx);
            else
                if do_bma_all
                    post = model.g_post;
                    indx = 1:nm;
                else
                    if bma_fam <= nfam && bma_fam > 0 && rem(bma_fam,1) == 0
                        indx = find(family.partition==bma_fam);
                        post = model.g_post(:,indx);
                    else
                        spm('alert','Incorrect family for BMA!','Error');
                        return
                    end
                end
            end
        else
            post  = model.g_post;
            indx = 1:nm;
        end
        
        disp('RFX Bayesian model averaging ...');
        % bayesian model averaging
        % ------------------------------------------------------------------
        bma = spm_dcm_bma(post,indx,subj,nsamp,odds_ratio);
        
    else
        
        bma = [];
        
    end
    
    if exist(fullfile(job.dir{1},'BMS.mat'),'file')
        load(fname);
        if  isfield(BMS,'DCM') && isfield(BMS.DCM,'rfx')
            spm('alert','Existing BMS.mat file has been over-written!',...
                'Warning');
        end
    end
    BMS.DCM.rfx.data    = fname_msp;
    BMS.DCM.rfx.F_fname = f_fname;
    BMS.DCM.rfx.F       = F;
    BMS.DCM.rfx.SF      = sumF;
    BMS.DCM.rfx.model   = model;
    BMS.DCM.rfx.family  = family;
    BMS.DCM.rfx.bma     = bma;
    
    disp('Saving BMS.mat file...')
    if spm_check_version('matlab','7') >= 0
        save(fname,'-V6','BMS');
    else
        save(fname,'BMS');
    end
    out.files{1}= fname;
    
end

% Save model_space
% -------------------------------------------------------------------------
if ~ld_msp && data_se && ~ld_f
    disp('Saving model space...')
    if spm_check_version('matlab','7') >= 0
        save(fname_msp,'-V6','subj');
    else
        save(fname_msp,'subj');
    end
end

% Data verification
% -------------------------------------------------------------------------
if job.verify_id
    disp('Data identity has been verified');
else
    disp('Data identity has not been verified');
end

end

