function spm_run_results(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_results.m 3095 2009-05-04 11:22:18Z volkmar $

cspec = job.conspec;
for k = 1:numel(cspec)
    job.conspec=cspec(k);
    if (numel(cspec(k).contrasts) == 1) && isinf(cspec(k).contrasts)
        tmp=load(job.spmmat{1});
        cspec1 = repmat(cspec(k),size(tmp.SPM.xCon));
        for l=1:numel(tmp.SPM.xCon)
            cspec1(l).contrasts = l;
        end;
        job1=job;
        job1.print = 1;
        job1.conspec = cspec1;
        spm_run_results(job1);
    else
        xSPM.swd        = spm_str_manip(job.spmmat{1},'H');
        xSPM.Ic         = job.conspec.contrasts;
        xSPM.u          = job.conspec.thresh;
        xSPM.Im         = [];
        if ~isempty(job.conspec.mask)
            xSPM.Im     = job.conspec.mask.contrasts;
            xSPM.pm     = job.conspec.mask.thresh;
            xSPM.Ex     = job.conspec.mask.mtype;
        end
        xSPM.thresDesc  = job.conspec.threshdesc;
        xSPM.title      = job.conspec.titlestr;
        xSPM.k          = job.conspec.extent;
        [hReg xSPM SPM] = spm_results_ui('Setup',xSPM);
        if job.print
            spm_list('List',xSPM,hReg);
            spm_figure('Print');
        end;
        assignin('base','hReg',hReg);
        assignin('base','xSPM',xSPM);
        assignin('base','SPM',SPM);
    end;
end;
