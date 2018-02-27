function out = spm_run_results(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2008-2018 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_run_results.m 7268 2018-02-27 16:48:13Z guillaume $


cspec = job.conspec;
for k = 1:numel(cspec)
    job.conspec = cspec(k);
    
    %-Apply to all contrasts if Inf is entered
    %----------------------------------------------------------------------
    if (numel(cspec(k).contrasts) == 1) && isinf(cspec(k).contrasts)
        tmp    = load(job.spmmat{1});
        cspec1 = repmat(cspec(k),size(tmp.SPM.xCon));
        for l = 1:numel(tmp.SPM.xCon)
            cspec1(l).contrasts = l;
        end
        job1         = job;
        %job1.print  = spm_get_defaults('ui.print');
        job1.conspec = cspec1;
        out = spm_run_results(job1);
        continue;
    end
    
    %-Create xSPM variable
    %----------------------------------------------------------------------
    xSPM.swd       = spm_file(job.spmmat{1},'fpath');
    xSPM.Ic        = job.conspec.contrasts;
    xSPM.u         = job.conspec.thresh;
    xSPM.Im        = [];
    if ~isfield(job.conspec.mask,'none')
        if isfield(job.conspec.mask,'contrast')
            xSPM.Im    = job.conspec.mask.contrast.contrasts;
            xSPM.pm    = job.conspec.mask.contrast.thresh;
            xSPM.Ex    = job.conspec.mask.contrast.mtype;
        elseif isfield(job.conspec.mask,'image')
            xSPM.Im    = job.conspec.mask.image.name;
            xSPM.pm    = [];
            xSPM.Ex    = job.conspec.mask.image.mtype;
        end
    end
    xSPM.thresDesc = job.conspec.threshdesc;
    xSPM.title     = job.conspec.titlestr;
    xSPM.k         = job.conspec.extent;
    try
        xSPM.n     = job.conspec.conjunction;
    end
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
    
    %-Compute a specified and thresholded SPM
    %----------------------------------------------------------------------
    if ~spm('CmdLine')
        [hReg, xSPM, SPM] = spm_results_ui('Setup',xSPM);
        TabDat = spm_list('List',xSPM,hReg);
    else
        [SPM, xSPM] = spm_getSPM(xSPM);
        TabDat = spm_list('Table',xSPM);
        hReg = [];
    end
    
    assignin('base', 'TabDat', TabDat);
    assignin('base', 'hReg',   hReg);
    assignin('base', 'xSPM',   xSPM);
    assignin('base', 'SPM',    SPM);
    
    out.xSPMvar(k)   = xSPM;
    out.TabDatvar(k) = TabDat;
    out.filtered{k}  = {};
    
    %-Export
    %----------------------------------------------------------------------
    for i=1:numel(job.export)
        fn = char(fieldnames(job.export{i}));
        switch fn
            case {'tspm','binary','nary'}
                fname = spm_file(xSPM.Vspm(1).fname,...
                    'suffix',['_' job.export{i}.(fn).basename]);
                descrip = sprintf('SPM{%c}-filtered: u = %5.3f, k = %d',...
                    xSPM.STAT,xSPM.u,xSPM.k);
                switch fn % see spm_results_ui.m
                    case 'tspm'
                        Z = xSPM.Z;
                    case 'binary'
                        Z = ones(size(xSPM.Z));
                    case 'nary'
                        if ~isfield(xSPM,'G')
                            Z       = spm_clusters(xSPM.XYZ);
                            num     = max(Z);
                            [n, ni] = sort(histc(Z,1:num), 2, 'descend');
                            n       = size(ni);
                            n(ni)   = 1:num;
                            Z       = n(Z);
                        else
                            C       = NaN(1,size(xSPM.G.vertices,1));
                            C(xSPM.XYZ(1,:)) = ones(size(xSPM.Z));
                            C       = spm_mesh_clusters(xSPM.G,C);
                            Z       = C(xSPM.XYZ(1,:));
                        end
                end
                if isfield(xSPM,'G')
                    M     = gifti(xSPM.G);
                    C     = zeros(1,size(xSPM.G.vertices,1));
                    C(xSPM.XYZ(1,:)) = Z; % or use NODE_INDEX
                    M.cdata = C;
                    F     = spm_file(fname,'path',xSPM.swd);
                    save(M,F);
                    cmd   = 'spm_mesh_render(''Disp'',''%s'')';
                else
                    V = spm_write_filtered(Z,xSPM.XYZ,xSPM.DIM,xSPM.M,...
                        descrip,fname);
                    cmd = 'spm_image(''display'',''%s'')';
                    F = V.fname;
                end
                out.filtered{k} = F;
                fprintf('Written %s\n',spm_file(F,'link',cmd));
                
            case {'csv','xls'}
                ofile = spm_file(fullfile(xSPM.swd,...
                    ['spm_' datestr(now,'yyyymmmdd') '.' fn]),'unique');
                spm_list([upper(fn) 'List'],TabDat,ofile);
                if strcmp(fn,'csv'), cmd = 'open(''%s'')';
                else                 cmd = 'winopen(''%s'')'; end
                fprintf('Saving results to:\n  %s\n',spm_file(ofile,'link',cmd));
                
            case 'montage'
                % see myslover() in spm_results_ui.m
                so = slover;
                so.img.vol = spm_vol(char(job.export{i}.montage.background));
                so.img.prop = 1;
                so = add_spm(so, xSPM);
                so.transform = job.export{i}.montage.orientation;
                so = fill_defaults(so);
                so.slices = job.export{i}.montage.slices;
                so.figure = spm_figure('GetWin', 'SliceOverlay');
                so = paint(so);
                %spm_print('',so.figure);
                
            case 'nidm'
                opts = struct('mod',job.export{i}.nidm.modality, ...
                    'space',job.export{i}.nidm.refspace,...
                    'group',struct('N',[job.export{i}.nidm.group.nsubj],...
                        'name',{{job.export{i}.nidm.group.label}}));
                nidmfile = spm_results_nidm(SPM,xSPM,TabDat,opts);
                fprintf('Exporting results in:\n  %s\n',nidmfile);
                
            otherwise
                if ~spm('CmdLine')
                    spm_figure('Print','Graphics','',fn);
                else
                    spm_list('TxtList',TabDat);
                end
        end
    end
    
end
