function out = spm_run_bms_dcm (varargin)

% API to compare DCMs on the basis of their log-evidences. Three methods
% are available to identify the best among alternative models:
%
%  (1) single subject BMS using Bayes factors
%     (see Penny et al, NeuroImage, 2004)
%  (2) fixed effects group BMS using group Bayes factors
%     (see Stephan et al,NeuroImage, 2007)
%  (3) random effects group BMS using exceedance probabilities
%     (see Stephan et al,NeuroImage, 2009)
%
% Note: All functions use the negative free energy (F) as an approximation
% to the log model evidence.
%
% __________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Chun-Chuan Chen
% $Id: spm_run_bms_dcm.m 3177 2009-06-03 08:47:41Z vladimir $

job     = varargin{1};
fname  ='BMS.mat';                  % Output filename
fname  = [job.dir{1},fname];        % Output filename (including directory)
F = [];
N = {};

% prepare the data
if  isempty(job.load_f{1})==0
    data=job.load_f{1};
    load(data);
    nm   = size(F,2);                               % No of Models
    N    = 1:nm;

else

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


% make inference and visualization
if strcmp(job.method,'FFX');                    %%% single subject BMS or 1st level( fixed effects) group BMS

    P = spm_api_bmc(sum(F,1),N);

    if exist(fullfile(job.dir{1},'BMS.mat'),'file')
        load(fname);
        if  isfield(BMS,'DCM') && isfield(BMS.DCM,'ffx')
            str = { 'Warning: existing BMS.mat file has been over-written!'};
            msgbox(str)
            BMS.DCM.ffx.F      = F;
            BMS.DCM.ffx.P      = P;
            BMS.DCM.ffx.SF     = sum(F,1);
            BMS.DCM.ffx.data   = data;
        else
            BMS.DCM.ffx.F     = F;
            BMS.DCM.ffx.P     = P;
            BMS.DCM.ffx.SF     = sum(F,1);
            BMS.DCM.ffx.data  = data;
        end
        save(fname,'BMS')
        out.files{1} = fname;
    else
        BMS.DCM.ffx.F     = F;
        BMS.DCM.ffx.P     = P;
        BMS.DCM.ffx.SF     = sum(F,1);
        BMS.DCM.ffx.data  = data;
        save(fname,'BMS')
        out.files{1} = fname;
    end
else
    % 2nd-level (random effects) BMS
    if  nm==2
        [alpha,exp_r,xp] = spm_BMS(F, 1e6, 1);
    else
        [alpha,exp_r,xp] = spm_BMS(F, 1e6, 0);
    end

    if exist(fullfile(job.dir{1},'BMS.mat'),'file')
        load(fname);
        if  isfield(BMS,'DCM') && isfield(BMS.DCM,'rfx')
            str = { 'Warning: existing BMS.mat file has been over-written!'};
            msgbox(str)
            BMS.DCM.rfx.F      = F;
            BMS.DCM.rfx.SF     = sum(F,1);
            BMS.DCM.rfx.alpha = alpha;
            BMS.DCM.rfx.exp_r = exp_r;
            BMS.DCM.rfx.xp    = xp;
            BMS.DCM.rfx.data  = data;
        else
            BMS.DCM.rfx.F      = F;
            BMS.DCM.rfx.SF     = sum(F,1);
            BMS.DCM.rfx.alpha  = alpha;
            BMS.DCM.rfx.exp_r  = exp_r;
            BMS.DCM.rfx.xp     = xp;
            BMS.DCM.rfx.data   = data;
        end
        save(fname,'BMS')
        out.files{1}= fname;
    else
        BMS.DCM.rfx.F      = F;
        BMS.DCM.rfx.SF     = sum(F,1);
        BMS.DCM.rfx.alpha = alpha;
        BMS.DCM.rfx.exp_r = exp_r;
        BMS.DCM.rfx.xp    = xp;
        BMS.DCM.rfx.data  = data;
        save(fname,'BMS')
        out.files{1}= fname;
    end

    P = spm_api_bmc(sum(F,1),N,alpha,exp_r,xp);  %%% display the result


end

spm_figure('GetWin','Graphics');
axes('position', [0.01, 0.01, 0.01, 0.01]); 
axis off
if job.verify_id
    text(0, 0, 'Data identity has been verified');
else
    text(0, 0, 'Data identity has not been verified');
end

end

