function spm_run_results(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_results.m 4384 2011-07-06 17:00:20Z guillaume $

cspec = job.conspec;
for k = 1:numel(cspec)
    job.conspec=cspec(k);
    if (numel(cspec(k).contrasts) == 1) && isinf(cspec(k).contrasts)
        tmp    = load(job.spmmat{1});
        cspec1 = repmat(cspec(k),size(tmp.SPM.xCon));
        for l=1:numel(tmp.SPM.xCon)
            cspec1(l).contrasts = l;
        end
        job1           = job;
        job1.print     = 1;
        job1.conspec   = cspec1;
        spm_run_results(job1);
    else
        xSPM.swd       = spm_str_manip(job.spmmat{1},'H');
        xSPM.Ic        = job.conspec.contrasts;
        xSPM.u         = job.conspec.thresh;
        xSPM.Im        = [];
        if ~isempty(job.conspec.mask)
            xSPM.Im    = job.conspec.mask.contrasts;
            xSPM.pm    = job.conspec.mask.thresh;
            xSPM.Ex    = job.conspec.mask.mtype;
        end
        xSPM.thresDesc = job.conspec.threshdesc;
        xSPM.title     = job.conspec.titlestr;
        xSPM.k         = job.conspec.extent;
        switch job.units
            case 1
                xSPM.units = {'mm' 'mm' 'mm'};
            case 2
                xSPM.units = {'mm' 'mm' 'ms'};
            case 3
                xSPM.units = {'mm' 'mm' 'Hz'};
            case 4
                xSPM.units = {'Hz' 'ms' ''};
            case 5
                xSPM.units = {'Hz' 'Hz' ''};
            otherwise
                error('Unknown data type.');
        end
        [hReg xSPM SPM] = spm_results_ui('Setup',xSPM);
        TabDat = spm_list('List',xSPM,hReg);
        if job.print
            spm_figure('Print');
        end
        assignin('base','TabDat',TabDat);
        assignin('base', 'hReg', hReg);
        assignin('base', 'xSPM', xSPM);
        assignin('base', 'SPM',  SPM);
    end
end
