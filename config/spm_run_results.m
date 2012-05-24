function out = spm_run_results(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2008-2011 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_results.m 4747 2012-05-24 11:11:10Z guillaume $

cspec = job.conspec;
for k = 1:numel(cspec)
    job.conspec=cspec(k);
    if (numel(cspec(k).contrasts) == 1) && isinf(cspec(k).contrasts)
        tmp    = load(job.spmmat{1});
        cspec1 = repmat(cspec(k),size(tmp.SPM.xCon));
        for l = 1:numel(tmp.SPM.xCon)
            cspec1(l).contrasts = l;
        end
        job1           = job;
        job1.print     = 1;
        job1.conspec   = cspec1;
        out = spm_run_results(job1);
    else
        xSPM.swd       = spm_file(job.spmmat{1},'fpath');
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
        if ~spm('CmdLine')
            [hReg, xSPM, SPM] = spm_results_ui('Setup',xSPM);
            TabDat = spm_list('List',xSPM,hReg);
            if job.print
                spm_figure('Print');
            end
        else
            [SPM,xSPM] = spm_getSPM(xSPM);
            TabDat = spm_list('Table',xSPM);
            if job.print
                spm_list('TxtList',TabDat);
            end
            hReg = [];
        end
        assignin('base', 'TabDat', TabDat);
        assignin('base', 'hReg',   hReg);
        assignin('base', 'xSPM',   xSPM);
        assignin('base', 'SPM',    SPM);
        out.xSPMvar(k)   = xSPM;
        out.TabDatvar(k) = TabDat;
        
        fn = fieldnames(job.write);
        switch fn{1}
            case 'none'
            case {'tspm','binary','nary'}
                if numel(xSPM.Ic)>1, continue; end
                fname = spm_file(xSPM.Vspm.fname,...
                    'suffix',['_' job.write.(fn{1}).basename]);
                out.filtered{k} = fname;
                descrip = sprintf('SPM{%c}-filtered: u = %5.3f, k = %d',...
                    xSPM.STAT,xSPM.u,xSPM.k);
                switch fn{1} % see spm_results_ui.m
                    case 'tspm'
                        Z = xSPM.Z;
                    case 'binary'
                        Z = ones(size(xSPM.Z));
                    case 'nary'
                        Z       = spm_clusters(xSPM.XYZ);
                        num     = max(Z);
                        [n, ni] = sort(histc(Z,1:num), 2, 'descend');
                        n       = size(ni);
                        n(ni)   = 1:num;
                        Z       = n(Z);
                end
                spm_write_filtered(Z,xSPM.XYZ,xSPM.DIM,xSPM.M,...
                    descrip,fname);
            otherwise
                error('Unknown option.');
        end
    end
end
