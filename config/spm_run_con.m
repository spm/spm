function out = spm_run_con(job)
% SPM job execution function - Specify and estimate contrasts
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - struct containing contrast and SPM{.} images filename
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_con.m 7093 2017-06-05 16:34:04Z guillaume $


spm('FnBanner','spm_contrasts.m');

%-Change to the analysis directory
%--------------------------------------------------------------------------
cwd = pwd;
try
    swd = spm_file(job.spmmat{1},'fpath');
    cd(swd);
    fprintf('%-40s: %30s\n','Contrasts folder',spm_file(swd,'short30'));%-#
catch
    error('Failed to change directory %s.',swd)
end

%-Load SPM.mat file
%--------------------------------------------------------------------------
load(fullfile(swd,'SPM.mat'));
SPM.swd = swd;

try
    SPM.xVol.XYZ;
catch
    error('This model has not been estimated.');
end

bayes_con = isfield(SPM,'PPM');

%-Delete contrast (if requested)
%--------------------------------------------------------------------------
if job.delete && isfield(SPM,'xCon')
    fprintf('%-40s: ','Deleting contrasts');                            %-#
    for k=1:numel(SPM.xCon)
        if ~isempty(SPM.xCon(k).Vcon)
            f = SPM.xCon(k).Vcon.fname;
            switch spm_file(f,'ext')
                case 'img'
                    n = spm_file(f,'basename');
                    spm_unlink([n '.img'],[n '.hdr']);
                case 'gii'
                    n = spm_file(f,'basename');
                    spm_unlink([n '.gii'],[n '.dat']);
                otherwise
                    spm_unlink(f);
            end
        end
        if ~isempty(SPM.xCon(k).Vspm)
            f = SPM.xCon(k).Vspm.fname;
            switch spm_file(f,'ext')
                case 'img'
                    n = spm_file(f,'basename');
                    spm_unlink([n '.img'],[n '.hdr']);
                case 'gii'
                    n = spm_file(f,'basename');
                    spm_unlink([n '.gii'],[n '.dat']);
                otherwise
                    spm_unlink(f);
            end
        end
    end
    
    SPM.xCon = [];
    if bayes_con, SPM.PPM.xCon = []; end
    
    %-Save SPM if no new contrasts are specified
    if isempty(job.consess)
        save(fullfile(SPM.swd,'SPM.mat'), 'SPM', spm_get_defaults('mat.format'));
    end
    fprintf('%30s\n','...done');                                        %-#
end

%-Retrospectively label Bayesian contrasts as T's, if this info is missing
%--------------------------------------------------------------------------
if bayes_con
    if ~isfield(SPM.PPM,'xCon')
        SPM.PPM.xCon = [];
        for ii=1:length(SPM.xCon)
            SPM.PPM.xCon(ii).PSTAT = 'T';
        end
    end
end

%-Specify contrasts
%--------------------------------------------------------------------------
Ic = [];
for i = 1:length(job.consess)
    
    %-T-contrast
    %----------------------------------------------------------------------
    if isfield(job.consess{i},'tcon')
        name = job.consess{i}.tcon.name;
        if bayes_con
            STAT = 'P';
            SPM.PPM.xCon(end+1).PSTAT = 'T';
            SPM.xX.V = [];
        else
            STAT = 'T';
        end
        con = job.consess{i}.tcon.weights(:)';
        sessrep = job.consess{i}.tcon.sessrep;
        
    %-T-contrast (cond/sess based)
    %----------------------------------------------------------------------
    elseif isfield(job.consess{i},'tconsess')
        name = job.consess{i}.tconsess.name;
        if bayes_con
            STAT = 'P';
            SPM.PPM.xCon(end+1).PSTAT = 'T';
            SPM.xX.V = [];
        else
            STAT = 'T';
        end
        
        %-T-contrast for conditions
        if isfield(job.consess{i}.tconsess.coltype,'colconds')
            ccond = job.consess{i}.tconsess.coltype.colconds;
            con = zeros(1,size(SPM.xX.X,2)); % overall contrast
            for cs = job.consess{i}.tconsess.sessions
                for k=1:numel(ccond)
                    if SPM.xBF.order < ccond(k).colbf
                        error(['Session-based contrast %d:\n'...
                            'Basis function order (%d) in design less ' ...
                            'than specified basis function number (%d).'],...
                            i, SPM.xBF.order, ccond(k).colbf);
                    end
                    % Index into columns belonging to the specified
                    % condition
                    try
                        cind = ccond(k).colbf + ...
                            ccond(k).colmodord*SPM.xBF.order ...
                            *SPM.Sess(cs).U(ccond(k).colcond).P(ccond(k) ...
                            .colmod).i(ccond(k).colmodord+1);
                        con(SPM.Sess(cs).col(SPM.Sess(cs).Fc(ccond(k).colcond).i(cind))) ...
                            = ccond(k).conweight;
                    catch
                        error(['Session-based contrast %d:\n'...
                            'Column "Cond%d Mod%d Order%d" does not exist.'],...
                            i, ccond(k).colcond, ccond(k).colmod, ccond(k).colmodord);
                    end
                end
            end
            
        %-T-contrast for extra regressors
        else
            con = zeros(1,size(SPM.xX.X,2)); % overall contrast
            for cs = job.consess{i}.tconsess.sessions
                nC = size(SPM.Sess(cs).C.C,2);
                if nC < numel(job.consess{i}.tconsess.coltype.colreg)
                    error(['Session-based contrast %d:\n'...
                        'Contrast vector for extra regressors too long.'],...
                        i);
                end
                ccols = numel(SPM.Sess(cs).col)-(nC-1)+...
                    [0:numel(job.consess{i}.tconsess.coltype.colreg)-1];
                con(SPM.Sess(cs).col(ccols)) = job.consess{i}.tconsess.coltype.colreg;
            end
        end
        
        sessrep = 'none';
        
    %-F-contrast
    %----------------------------------------------------------------------
    else
        name = job.consess{i}.fcon.name;
        if bayes_con
            STAT = 'P';
            SPM.PPM.xCon(end+1).PSTAT = 'F';
            SPM.xX.V = [];
        else
            STAT = 'F';
        end
        con = job.consess{i}.fcon.weights;
        sessrep = job.consess{i}.fcon.sessrep;
    end

    %-Replicate contrast over sessions
    %----------------------------------------------------------------------
    if isfield(SPM,'Sess') && ~strcmp(sessrep,'none')
        nsessions = numel(SPM.Sess); % assume identical sessions
        switch sessrep
            case {'repl','replsc'}
                % within-session zero padding, replication over sessions
                cons = {zeros(size(con,1),size(SPM.xX.X,2))};
                for sess=1:nsessions
                    sfirst = SPM.Sess(sess).col(1);
                    cons{1}(:,sfirst:sfirst+size(con,2)-1) = con;
                end
                if strcmp(sessrep,'replsc')
                    cons{1} = cons{1} / nsessions;
                end
                names = {sprintf('%s - All Sessions', name)};
            case 'replna'
                % within-session zero padding, new rows per session
                cons = {zeros(nsessions*size(con,1),size(SPM.xX.X,2))};
                for sess=1:nsessions
                    sfirst = SPM.Sess(sess).col(1);
                    cons{1}((sess-1)*size(con,1)+(1:size(con,1)),sfirst-1+(1:size(con,2)))=con;
                end
                names = {sprintf('%s - All Sessions', name)};
            case 'sess'
                cons = cell(1,numel(SPM.Sess));
                names = cell(1,numel(SPM.Sess));
                for k=1:numel(SPM.Sess)
                    cons{k} = [zeros(size(con,1),SPM.Sess(k).col(1)-1) con];
                    names{k} = sprintf('%s - Session %d', name, k);
                end
            case {'both','bothsc'}
                cons = cell(1,numel(SPM.Sess));
                names = cell(1,numel(SPM.Sess));
                for k=1:numel(SPM.Sess)
                    cons{k} = [zeros(size(con,1),SPM.Sess(k).col(1)-1) con];
                    names{k} = sprintf('%s - Session %d', name, k);
                end
                if numel(SPM.Sess) > 1
                    % within-session zero padding, replication over sessions
                    cons{end+1} = zeros(size(con,1),size(SPM.xX.X,2));
                    for sess=1:nsessions
                        sfirst = SPM.Sess(sess).col(1);
                        cons{end}(:,sfirst:sfirst+size(con,2)-1)=con;
                    end
                    if strcmp(sessrep,'bothsc')
                        cons{end} = cons{end}/nsessions;
                    end
                    names{end+1} = sprintf('%s - All Sessions', name);
                end
        end
    else
        cons = {con};
        names = {name};
    end

    %-Store newly created contrasts in SPM.xCon
    %----------------------------------------------------------------------
    for k=1:numel(cons)

        %-Basic checking of contrast
        %------------------------------------------------------------------
        [c,I,emsg,imsg] = spm_conman('ParseCon',cons{k},SPM.xX.xKXs,STAT);
        if ~isempty(emsg)
            disp(emsg);
            error('Error in contrast specification');
        else
            %disp(imsg);
        end

        %-Fill-in the contrast structure
        %------------------------------------------------------------------
        if all(I)
            DxCon = spm_FcUtil('Set',names{k},STAT,'c',c,SPM.xX.xKXs);
        else
            DxCon = [];
        end

        %-Append to SPM.xCon
        %------------------------------------------------------------------
        if isempty(SPM.xCon)
            SPM.xCon = DxCon;
        elseif ~isempty(DxCon)
            SPM.xCon(end+1) = DxCon;
        end
        Ic = [Ic length(SPM.xCon)];
    end
end

%-Estimate newly created contrasts (and save SPM.mat)
%--------------------------------------------------------------------------
if ~isempty(Ic), SPM = spm_contrasts(SPM,Ic); end

fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#

%-Change back directory
%--------------------------------------------------------------------------
cd(cwd); 

%-Output structure
%--------------------------------------------------------------------------
out.spmmat = job.spmmat;
%out.spmvar = SPM;
if isfield(SPM, 'xCon') && ~isempty(SPM.xCon)
    Vcon = [SPM.xCon.Vcon]; % [SPM.xCon(Ic).Vcon] ?
    Vspm = [SPM.xCon.Vspm]; % [SPM.xCon(Ic).Vspm] ?
elseif isfield(SPM, 'PPM') && ~isempty(SPM.PPM)
    Vcon = cat(1,SPM.PPM.xCon.Vcon);
    Vspm = cat(1,SPM.PPM.xCon.Vspm);
else
    Vcon = ''; Vspm = '';
end
if ~isempty(Vcon) && ~isempty(Vspm)
    out.con = spm_file(cellstr(char(Vcon.fname)),'path',swd);
    out.spm =  spm_file(cellstr(char(Vspm.fname)),'path',swd);
else
    out.con = {}; out.spm = {};
end
