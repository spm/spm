function out = spm_run_con(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_con.m 3993 2010-07-13 11:59:32Z volkmar $


wd  = pwd;
job = varargin{1};

% Change to the analysis directory
%-----------------------------------------------------------------------
if ~isempty(job)
    try
        pth = fileparts(job.spmmat{:});
        cd(char(pth));
        fprintf('   Changing directory to: %s\n',char(pth));
    catch
        error('Failed to change directory. Aborting contrast setup.')
    end
end

% Load SPM.mat file
%-----------------------------------------------------------------------
load(job.spmmat{:},'SPM');

try
    SPM.xVol.XYZ;
catch
    error('This model has not been estimated.');
end

if ~strcmp(pth,SPM.swd)
    warning(['Path to SPM.mat: %s\n and SPM.swd: %s\n differ, using current ' ...
             'SPM.mat location as new working directory.'], pth, ...
            SPM.swd);
    SPM.swd = pth;
end

if job.delete && isfield(SPM,'xCon')
    for k=1:numel(SPM.xCon)
        if ~isempty(SPM.xCon(k).Vcon)
            [p n e v] = spm_fileparts(SPM.xCon(k).Vcon.fname);
            switch e,
                case '.img'
                    spm_unlink([n '.img'],[n '.hdr']);
                case '.nii'
                    spm_unlink(SPM.xCon(k).Vcon.fname);
            end
        end
        if ~isempty(SPM.xCon(k).Vspm)
            [p n e v] = spm_fileparts(SPM.xCon(k).Vspm.fname);
            switch e,
                case '.img'
                    spm_unlink([n '.img'],[n '.hdr']);
                case '.nii'
                    spm_unlink(SPM.xCon(k).Vspm.fname);
            end
        end
    end
    SPM.xCon = [];
end

bayes_con=isfield(SPM,'PPM');
if bayes_con
    if ~isfield(SPM.PPM,'xCon')
        % Retrospectively label Bayesian contrasts as T's, if this info is missing
        for ii=1:length(SPM.xCon)
            SPM.PPM.xCon(ii).PSTAT='T';
        end
    end
end

for i = 1:length(job.consess)
    if isfield(job.consess{i},'tcon')
        name = job.consess{i}.tcon.name;
        if bayes_con
            STAT = 'P';
            SPM.PPM.xCon(end+1).PSTAT = 'T';
            SPM.xX.V=[];
        else
            STAT = 'T';
        end
        con  = job.consess{i}.tcon.convec(:)';
        sessrep = job.consess{i}.tcon.sessrep;
    elseif isfield(job.consess{i},'tconsess')
        job.consess{i}.tconsess = job.consess{i}.tconsess; % save some typing
        name = job.consess{i}.tconsess.name;
        if bayes_con
            STAT = 'P';
            SPM.PPM.xCon(end+1).PSTAT = 'T';
            SPM.xX.V=[];
        else
            STAT = 'T';
        end
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
                    end;
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
                    end;
                end;
            end;
        else % convec on extra regressors
            con = zeros(1,size(SPM.xX.X,2)); % overall contrast
            for cs = job.consess{i}.tconsess.sessions
                nC = size(SPM.Sess(cs).C.C,2);
                if nC < numel(job.consess{i}.tconsess.coltype.colreg)
                    error(['Session-based contrast %d:\n'...
                        'Contrast vector for extra regressors too long.'],...
                        i);
                end;
                ccols = numel(SPM.Sess(cs).col)-(nC-1)+...
                    [0:numel(job.consess{i}.tconsess.coltype.colreg)-1];
                con(SPM.Sess(cs).col(ccols)) = job.consess{i}.tconsess.coltype.colreg;
            end;
        end;
        sessrep = 'none';
    else %fcon
        name = job.consess{i}.fcon.name;
        if bayes_con
            STAT = 'P';
            SPM.PPM.xCon(end+1).PSTAT = 'F';
            SPM.xX.V=[];
        else
            STAT = 'F';
        end
        try
            con  = cat(1,job.consess{i}.fcon.convec{:});
        catch
            error('Error concatenating F-contrast vectors. Sizes are:\n %s\n',... 
                   num2str(cellfun('length',job.consess{i}.fcon.convec)))
        end
        sessrep = job.consess{i}.fcon.sessrep;
    end

  
    if isfield(SPM,'Sess') && ~strcmp(sessrep,'none')
        % assume identical sessions, no check!
        nsessions=numel(SPM.Sess);
        switch sessrep
            case {'repl','replsc'}
                % within-session zero padding, replication over sessions
                cons = {zeros(size(con,1),size(SPM.xX.X,2))};
                for sess=1:nsessions
                    sfirst=SPM.Sess(sess).col(1);
                    cons{1}(:,sfirst:sfirst+size(con,2)-1)=con;
                end
                if strcmp(sessrep,'replsc')
                    cons{1} = cons{1}/nsessions;
                end
                names = {sprintf('%s - All Sessions', name)};
            case 'replna',
                % within-session zero padding, new rows per session
                cons= {zeros(nsessions*size(con,1),size(SPM.xX.X,2))};
                for sess=1:nsessions
                    sfirst=SPM.Sess(sess).col(1);
                    cons{1}((sess-1)*size(con,1)+(1:size(con,1)),sfirst-1+(1:size(con,2)))=con;
                end
                names = {sprintf('%s - All Sessions', name)};
            case 'sess',
                cons = cell(1,numel(SPM.Sess));
                names = cell(1,numel(SPM.Sess));
                for k=1:numel(SPM.Sess)
                    cons{k} = [zeros(size(con,1),SPM.Sess(k).col(1)-1) con];
                    names{k} = sprintf('%s - Session %d', name, k);
                end;
            case {'both','bothsc'}
                cons = cell(1,numel(SPM.Sess));
                names = cell(1,numel(SPM.Sess));
                for k=1:numel(SPM.Sess)
                    cons{k} = [zeros(size(con,1),SPM.Sess(k).col(1)-1) con];
                    names{k} = sprintf('%s - Session %d', name, k);
                end;
                if numel(SPM.Sess) > 1
                    % within-session zero padding, replication over sessions
                    cons{end+1}= zeros(size(con,1),size(SPM.xX.X,2));
                    for sess=1:nsessions
                        sfirst=SPM.Sess(sess).col(1);
                        cons{end}(:,sfirst:sfirst+size(con,2)-1)=con;
                    end
                    if strcmp(sessrep,'bothsc')
                        cons{end} = cons{end}/nsessions;
                    end
                    names{end+1} = sprintf('%s - All Sessions', name);
                end;
        end;
    else
        cons = {con};
        names = {name};
    end;

    % Loop over created contrasts
    %-------------------------------------------------------------------
    for k=1:numel(cons)

        % Basic checking of contrast
        %-------------------------------------------------------------------
        [c,I,emsg,imsg] = spm_conman('ParseCon',cons{k},SPM.xX.xKXs,STAT);
        if ~isempty(emsg)
            disp(emsg);
            error('Error in contrast specification');
        else
            disp(imsg);
        end;

        % Fill-in the contrast structure
        %-------------------------------------------------------------------
        if all(I)
            DxCon = spm_FcUtil('Set',names{k},STAT,'c',c,SPM.xX.xKXs);
        else
            DxCon = [];
        end

        % Append to SPM.xCon. SPM will automatically save any contrasts that
        % evaluate successfully.
        %-------------------------------------------------------------------
        if isempty(SPM.xCon)
            SPM.xCon = DxCon;
        elseif ~isempty(DxCon)
            SPM.xCon(end+1) = DxCon;
        end
        SPM = spm_contrasts(SPM,length(SPM.xCon));
    end
end;
% Change back directory
%-----------------------------------------------------------------------
fprintf('   Changing back to directory: %s\n', wd);
cd(wd); 
out.spmmat = job.spmmat;
%out.spmvar = SPM;
if isfield(SPM, 'xCon')
    Vcon = cat(1,SPM.xCon.Vcon);
    Vspm = cat(1,SPM.xCon.Vspm);
elseif isfield(SPM, 'PPM')
    Vcon = cat(1,SPM.PPM.xCon.Vcon);
    Vspm = cat(1,SPM.PPM.xCon.Vspm);
end;
out.con = cellstr(strvcat(Vcon.fname));
out.spm = cellstr(strvcat(Vspm.fname));
