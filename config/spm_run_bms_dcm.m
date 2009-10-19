function out = spm_run_bms_dcm (varargin)
% API to compare DCMs on the basis of their log-evidences. Four methods
% are available to identify the best among alternative models:
%
%  (1) single subject BMS using Bayes factors
%     (see Penny et al, NeuroImage, 2004)
%  (2) fixed effects group BMS using group Bayes factors
%     (see Stephan et al,NeuroImage, 2007)
%  (3) random effects group BMS using exceedance probabilities
%     (see Stephan et al,NeuroImage, 2009)
%  (4) comparing model families
%     (see Penny et al, PLOS-CB, submitted) 
%
% Note: All functions use the negative free energy (F) as an approximation
% to the log model evidence.
%
% __________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Chun-Chuan Chen
% $Id: spm_run_bms_dcm.m 3479 2009-10-19 10:10:55Z maria $

% input
% -------------------------------------------------------------------------
job     = varargin{1};
fname   = 'BMS.mat';                 % Output filename
fname   = fullfile(job.dir{1},fname);% Output filename (including directory)
priors  = job.priors;
ld_f    = ~isempty(job.load_f{1});
bma_do  = isfield(job.bma,'bma_yes');
data_se = ~isempty(job.sess_dcm);

do_bma_famwin = 0;
do_bma_all    = 0;

% check DCM.mat files and BMA 
% -------------------------------------------------------------------------
if bma_do
    if data_se
       load(job.sess_dcm{1}(1).mod_dcm{1})
       n  = size(DCM.a,2);
       m  = size(DCM.c,2); 
       mi = size(DCM.c,2); 
       bma.nsamp       = str2num(job.bma.bma_yes.bma_nsamp);
       bma.odds_ratio  = str2num(job.bma.bma_yes.bma_ratio);
       bma.a           = zeros(n,n,bma.nsamp);
       bma.b           = zeros(n,n,m,bma.nsamp);
       bma.c           = zeros(n,mi,bma.nsamp);
       
       if isfield(job.bma.bma_yes.bma_set,'bma_famwin')
           do_bma_famwin = 1;
       else
           if isfield(job.bma.bma_yes.bma_set,'bma_all')
               do_bma_all = 1;
           else
               bma_fam    = double(job.bma.bma_yes.bma_set.bma_part);
           end
       end
    else
        error('Plase specify DCM.mat files to do BMA!')
    end
end

F = [];
N = {};

% prepare the data
% -------------------------------------------------------------------------
if  ld_f
    data = job.load_f{1};
    load(data);
    nm      = size(F,2);                               % No of Models
    ns      = size(F,1);                               % No of Models
    N       = 1:nm;
    subj    = [];
    f_fname = data;
    
else

    f_fname = [];
    ns      = size(job.sess_dcm,2);                 % No of Subjects
    nsess   = size(job.sess_dcm{1},2);              % No of sessions
    nm      = size(job.sess_dcm{1}(1).mod_dcm,1);   % No of Models
    
    % Check if No of models > 2
    if nm < 2
        msgbox('Please select more than one file')
        return
    end
    
    for k=1:ns
        
        data{k}         = [job.sess_dcm{k}(:).mod_dcm];
        nsess_now       = size(job.sess_dcm{k},2);
        nmodels         = size(job.sess_dcm{k}(1).mod_dcm,1);
        
        if (nsess_now == nsess && nmodels== nm) % Check no of sess/mods
            
            ID = zeros(nsess, nm);
            
            for j=1:nm
                
                F_sess      = [];
                
                for h = 1:nsess_now
                    
                    tmp     = data{k}(j,h);
                    DCM_tmp = load(tmp{1});
                    
                    F_sess  = [F_sess,DCM_tmp.DCM.F];
                    subj(k).sess(h).model(j).fname = tmp{1};
                    
                    % Data ID verification. At least for now we'll
                    % re-compute the IDs rather than use the ones stored
                    % with the DCM.
                    if job.verify_id
                        M = DCM_tmp.DCM.M;
                        
                        if isfield(DCM_tmp.DCM, 'xY')
                            Y = DCM_tmp.DCM.xY;  %not fMRI
                        else
                            Y = DCM_tmp.DCM.Y;   % fMRI
                        end
                        
                        if isfield(M,'FS')
                            try
                                ID(h, j)  = spm_data_id(feval(M.FS,Y.y,M));
                            catch
                                ID(h, j)  = spm_data_id(feval(M.FS,Y.y));
                            end
                        else
                            ID(h, j) = spm_data_id(Y.y);
                        end
                        
                    end
                end
                
                F_mod       = sum(F_sess);
                F(k,j)      = F_mod;
                N{j}        = sprintf('model%d',j);
                
            end
            
            if job.verify_id
                failind = find(max(abs(diff(ID))) > eps);
                if ~isempty(failind)
                    out.files{1} = [];
                    msgbox(['Error: the models for subject ' num2str(k) ...
                        ' session(s) ' num2str(failind) ' were not fitted to the same data.']);
                    return
                end
            end
        else
            out.files{1} = [];
            msgbox('Error: the number of sessions/models should be the same for all subjects!')
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
        family.prior = priors;
        family.infer = job.method;
    
        nfam    = size(family.names,2);
        npart   = length(unique(family.partition));
        maxpart = max(family.partition);
        m_indx  = 1:nm;
    
        if nfam ~= npart || npart == 1 || maxpart > npart
            error('Invalid family file!')
            out.files{1} = [];
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
           error('Invalid family!')
           out.files{1} = [];
           return
        end
        
        family.prior     = priors;
        family.infer     = job.method;
      
    end
    
end

% single subject BMS or 1st level( fixed effects) group BMS
% -------------------------------------------------------------------------
if strcmp(job.method,'FFX'); 
    
    model.post     = spm_api_bmc(sumF,N);
    
    if ~do_family
        family     = [];
    else
        Ffam           = F(:,m_indx);
        [family,model] = spm_compare_families (Ffam,family);     
    end
    
    if bma_do,  
       if do_bma_famwin
           [fam_max,fam_max_i]  = max(family.post);
           post_indx = find(family.partition==fam_max_i);
           bma.post  = model.post(post_indx);
       else
           if do_bma_all  
              bma.post = model.post;
           else
               if bma_fam <= nfam && bma_fam > 0 && rem(bma_fam,1) == 0 
                    post_indx = find(family.partition==bma_fam);
                    bma.post  = model.post(post_indx);
               else
                    error('Incorrect family for BMA!');
               end
                   
           end
       end

       % bayesian model averaging
       % ------------------------------------------------------------------
       theta = spm_dcm_bma(bma.post,post_indx,subj,bma.nsamp,bma.odds_ratio);

       % reshape parameters
       % ------------------------------------------------------------------
       for i = 1:bma.nsamp,
           [A,B,C]     = spm_dcm_reshape(theta(:,i),m,n,1);
           bma.a(:,:,i)   = A(:,:);
           bma.b(:,:,:,i) = B(:,:,:);
           bma.c(:,:,i)   = C(:,:);
       end
       
    else
        
        bma = [];
        
    end

    if exist(fullfile(job.dir{1},'BMS.mat'),'file')
        load(fname);
        if  isfield(BMS,'DCM') && isfield(BMS.DCM,'ffx')
            str = { 'Warning: existing BMS.mat file has been over-written!'};
            msgbox(str)
        end
    end
    BMS.DCM.ffx.data    = subj;
    BMS.DCM.ffx.F_fname = f_fname;
    BMS.DCM.ffx.F       = F;
    BMS.DCM.ffx.SF      = sumF;
    BMS.DCM.ffx.model   = model;
    BMS.DCM.ffx.family  = family;
    BMS.DCM.ffx.bma     = bma;
     
    save(fname,'BMS')
    out.files{1} = fname;
    
% 2nd-level (random effects) BMS
% -------------------------------------------------------------------------
else   
    
    if ~do_family
       
        if nm <= ns
            [alpha,exp_r,xp] = spm_BMS(F, 1e6, 0);
        else
            [exp_r,xp] = spm_BMS_gibbs(F);
        end
        model.alpha = [];
        model.exp_r = exp_r;
        model.xp    = xp;
        family      = [];
    else
        
        Ffam           = F(:,m_indx);
        [family,model] = spm_compare_families (Ffam,family);
        
    end
    
    if bma_do,  
       if do_bma_famwin
           [fam_max,fam_max_i]  = max(family.exp_r);
           post_indx = find(family.partition==fam_max_i);
           bma.post  = model.exp_r(post_indx);
       else
           if do_bma_all  
              bma.post = model.exp_r;
           else
              if bma_fam <= nfam && bma_fam > 0 && rem(bma_fam,1) == 0 
                post_indx = find(family.partition==bma_fam);
                bma.post = model.exp_r(post_indx);
              else
                    error('Incorrect family for BMA!');
               end
           end
       end

       % bayesian model averaging
       % ------------------------------------------------------------------
       theta = spm_dcm_bma(bma.post,post_indx,subj,bma.nsamp,bma.odds_ratio);

       % reshape parameters
       % ------------------------------------------------------------------
       for i = 1:bma.nsamp,
           [A,B,C]     = spm_dcm_reshape(theta(:,i),m,n,1);
           bma.a(:,:,i)   = A(:,:);
           bma.b(:,:,:,i) = B(:,:,:);
           bma.c(:,:,i)   = C(:,:);
       end

    else
        
        bma = [];
        
    end

    if exist(fullfile(job.dir{1},'BMS.mat'),'file')
        load(fname);
        if  isfield(BMS,'DCM') && isfield(BMS.DCM,'rfx')
            str = { 'Warning:  existing BMS.mat file has been over-written!'};
            msgbox(str)
        end
    end
    BMS.DCM.rfx.data    = subj;
    BMS.DCM.rfx.F_fname = f_fname;  
    BMS.DCM.rfx.F       = F;
    BMS.DCM.rfx.SF      = sumF;
    BMS.DCM.rfx.model   = model;
    BMS.DCM.rfx.family  = family;
    BMS.DCM.rfx.bma     = bma;
    
    save(fname,'BMS')
    out.files{1}= fname;

    % display the result
    P = spm_api_bmc(sumF,N,model.exp_r,model.xp); 
    
end

% verify data
% -------------------------------------------------------------------------
spm_figure('GetWin','Graphics');
axes('position', [0.01, 0.01, 0.01, 0.01]); 
axis off
if job.verify_id
    text(0, 0, 'Data identity has been verified');
else
    text(0, 0, 'Data identity has not been verified');
end

end

