function out = spm_run_bms_dcm (varargin)

% API to compare DCMs/Log-evidences and allows one to indentify 
% the best model among models being tested 
% 
% This function can report the results from 
%  (1) the single subject BMC using bayes factors (see Penny et al,NeuroImage, 2004)
%  (2) the lst-level group BMC using ffx method (see Penny et al,NeuroImage, 2004)           
%  (3) the 2nd-level group BMC using rfx method (see Stephan et al,NeuroImage, 2009)
%              
% __________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging 

% Chun-Chuan Chen
% $Id: spm_run_bms_dcm.m 2625 2009-01-20 16:28:20Z maria $

job     = varargin{1};
ns      = size(job.sess_dcm,2);                 % No of Subjects
nm      = size(job.sess_dcm{1}(1).mod_dcm,1);   % No of Models

fname  ='BMS.mat';                  % Output filename
fname  = [job.dir{1},fname];        % Output filename (including directory)

% Check if No of models > 2
if nm < 2
    msgbox('Please select more than one file')
    return
end

if ns>1
    group_inx=1;  
else
    group_inx=0;
end

switch group_inx

    case (0) %% the single subject BMC

        F = [];
        N = {};

        data        = [job.sess_dcm{1}(:).mod_dcm];
        nsess       = size(job.sess_dcm{1},2);

        for j=1:nm

            F_sess      = [];

            for h = 1:nsess

                tmp     = data(j,h);
                DCM_tmp = load(tmp{1});
                F_sess  = [F_sess,DCM_tmp.DCM.F];

            end

            F_mod       = sum(F_sess);
            F(:,j)      = F_mod;
            N{j}        = sprintf('model%d',j);

        end

        P = spm_api_bmc(F,N);

        if exist(fullfile(job.dir{1},'BMS.mat'),'file')
            load(fname);
            if  isfield(BMS,'DCM') && isfield(BMS.DCM,'single')
                str = { 'Warning: existing BMS.mat file has been over-written!'};
                msgbox(str)
                BMS.DCM.single.F = F;
                BMS.DCM.single.P = P;
            else
                BMS.DCM.single.F = F;
                BMS.DCM.single.P = P;
            end
            save(fname,'BMS')
            out.files{1} = fname;
        else
            BMS.DCM.single.F = F;
            BMS.DCM.single.P = P;
            save(fname,'BMS')
            out.files{1} = fname;
        end


case (1) %% the group BMC

    F = [];
    N = {};

    nsess1 = size(job.sess_dcm{1},2);
    %%% prepare the data
    for k=1:ns

        data        = [job.sess_dcm{k}(:).mod_dcm];
        nsess       = size(job.sess_dcm{k},2);
        nmodels     = size(job.sess_dcm{k}(1).mod_dcm,1);

        if (nsess == nsess1 && nm == nmodels) % Check no of sess/mods

            for j=1:nm

                F_sess      = [];

                for h = 1:nsess

                    tmp     = data(j,h);
                    DCM_tmp = load(tmp{1});
                    F_sess  = [F_sess,DCM_tmp.DCM.F];

                end

                F_mod       = sum(F_sess);
                F(k,j)      = F_mod;
                N{j}        = sprintf('model%d',j);

            end
        else

            out.files{1} = [];
            msgbox('Error: the number of sessions/models should be the same for all subjects!')
            return

        end

    end
  %% make inference
    if strcmp(job.method,'FFX');                    %%% 1st-level group BMS

        P = spm_api_bmc(sum(F,1),N);
       
        if exist(fullfile(job.dir{1},'BMS.mat'),'file')
            load(fname);
            if  isfield(BMS,'DCM') && isfield(BMS.DCM.group,'ffx')
                str = { 'Warning: existing BMS.mat file has been over-written!'};
                 msgbox(str)
                BMS.DCM.group.ffx.F = F;
                BMS.DCM.group.ffx.P = P;
            else
                BMS.DCM.group.ffx.F = F;
                BMS.DCM.group.ffx.P = P;
            end
            save(fname,'BMS')
            out.files{1} = fname;
        else
            BMS.DCM.group.ffx.F = F;
            BMS.DCM.group.ffx.P = P;
            save(fname,'BMS')
            out.files{1} = fname;
        end
    else
        % 2nd-level
        if  nm==2
            [alpha,exp_r,xp] = spm_BMS(F, 1e6, 1);
        else
            [alpha,exp_r,xp] = spm_BMS(F, 1e6, 0);
        end

         if exist(fullfile(job.dir{1},'BMS.mat'),'file')
            load(fname);
            if  isfield(BMS,'DCM') && isfield(BMS.DCM.group,'rfx')
                str = { 'Warning: existing BMS.mat file has been over-written!'};
                msgbox(str)
                BMS.DCM.group.rfx.alpha = alpha;
                BMS.DCM.group.rfx.exp_r = exp_r;
                BMS.DCM.group.rfx.xp    = xp;
            else
                BMS.DCM.group.rfx.alpha = alpha;
                BMS.DCM.group.rfx.exp_r = exp_r;
                BMS.DCM.group.rfx.xp    = xp;
            end
            save(fname,'BMS')
            out.files{1}= fname;
         else
            BMS.DCM.group.rfx.alpha = alpha;
            BMS.DCM.group.rfx.exp_r = exp_r;
            BMS.DCM.group.rfx.xp    = xp;
            save(fname,'BMS')
            out.files{1}= fname;
         end

        P = spm_api_bmc(sum(F,1),N,alpha,exp_r,xp);  %%% display the result


    end
    end

