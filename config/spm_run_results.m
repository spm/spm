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

% $Id: spm_run_results.m 1185 2008-03-04 16:31:21Z volkmar $

cspec = job.conspec;
for k = 1:numel(cspec)
    job.conspec=cspec(k);
    if (numel(cspec(k).contrasts) == 1) && isinf(cspec(k).contrasts)
    tmp=load(job.spmmat{1});
    for l=1:numel(tmp.SPM.xCon)
        cspec1(l) = cspec(k);
        cspec1(l).contrasts = l;
    end;
    job1=job;
    job1.print = 1;
    job1.conspec = cspec1;
    spm_run_results(job1);
    else
    [hReg xSPM SPM]=spm_results_ui('Setup',job);
    if job.print
        spm_list('List',xSPM,hReg);
        spm_figure('Print');
    end;
    assignin('base','hReg',hReg);
    assignin('base','xSPM',xSPM);
    assignin('base','SPM',SPM);
    end;
end;
