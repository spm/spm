function out = spm_run_fmri_spec(job)
% Setting up the general linear model for fMRI time-series
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2005-2019 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_fmri_spec.m 7739 2019-12-02 14:00:18Z guillaume $


%-Check presence of previous analysis
%==========================================================================

%-Output directory
%--------------------------------------------------------------------------
cwd = pwd;
d   = spm_file(job.dir{1},'cpath');
if ~exist(d,'dir')
    sts = mkdir(d);
    if ~sts, error('Error creating output directory "%s".',d); end
end
cd(d);

%-Ask about overwriting files from previous analyses
%--------------------------------------------------------------------------
if exist(fullfile(pwd,'SPM.mat'),'file')
    str = {'Current directory contains existing SPM file:',...
           'Continuing will overwrite existing file!'};
    if spm_input(str,1,'bd','stop|continue',[1,0],1,mfilename)
        fprintf('%-40s: %30s\n\n',...
            'Abort...   (existing SPM file)',spm('time'));
        out = []; return
    end
end

%-Delete old analysis files
%--------------------------------------------------------------------------
files = {'^mask\..{3}$','^ResMS\..{3}$','^RPV\..{3}$',...
    '^beta_.{4}\..{3}$','^con_.{4}\..{3}$','^ResI_.{4}\..{3}$',...
    '^ess_.{4}\..{3}$', '^spm\w{1}_.{4}\..{3}$'};

for i=1:length(files)
    j = spm_select('List',pwd,files{i});
    for k=1:size(j,1)
        spm_unlink(deblank(j(k,:)));
    end
end


%-Timing parameters & Basis functions
%==========================================================================

%-Repeat time
%--------------------------------------------------------------------------
SPM.xY.RT = job.timing.RT;

%-Basis function parameters
%--------------------------------------------------------------------------
SPM.xBF.UNITS = job.timing.units;
SPM.xBF.T     = job.timing.fmri_t;
SPM.xBF.T0    = job.timing.fmri_t0;

%-Basis functions
%--------------------------------------------------------------------------
bf = char(fieldnames(job.bases));
if strcmp(bf,'none')
    SPM.xBF.name = 'NONE';
elseif strcmp(bf,'hrf')
    if all(job.bases.hrf.derivs == [0 0])
        SPM.xBF.name = 'hrf';
    elseif all(job.bases.hrf.derivs == [1 0])
        SPM.xBF.name = 'hrf (with time derivative)';
    elseif all(job.bases.hrf.derivs == [1 1])
        SPM.xBF.name = 'hrf (with time and dispersion derivatives)';
    else
        error('Unknown HRF derivative choices.');
    end
else
    switch bf
        case 'fourier'
            SPM.xBF.name = 'Fourier set';
        case 'fourier_han'
            SPM.xBF.name = 'Fourier set (Hanning)';
        case 'gamma'
            SPM.xBF.name = 'Gamma functions';
        case 'fir'
            SPM.xBF.name = 'Finite Impulse Response';
        otherwise
            error('Unknown basis functions.');
    end
    SPM.xBF.length = job.bases.(bf).length;
    SPM.xBF.order  = job.bases.(bf).order;
end

%-Model interactions (Volterra)
%--------------------------------------------------------------------------
SPM.xBF.Volterra = job.volt;


%-Data & Design
%==========================================================================

design_only = ~isfield(job,'mask');

if ~design_only
    SPM.xY.P = [];
end

for i = 1:numel(job.sess)
    
    sess = job.sess(i);

    %-Image filenames
    %----------------------------------------------------------------------
    if design_only
        SPM.nscan(i) = sess.nscan;
    else
        SPM.nscan(i) = numel(sess.scans);
        if SPM.nscan(i) == 1
            sess.scans   = cellstr(spm_select('Expand',sess.scans{1}));
            SPM.nscan(i) = numel(sess.scans);
        end
        SPM.xY.P     = strvcat(SPM.xY.P,sess.scans{:});
    end
    if SPM.nscan(i) == 1
        error('Not enough scans in session %d.',i);
    end

    %-Multiple conditions (structure from a MAT-file)
    %----------------------------------------------------------------------
    if ~isempty(sess.multi{1})
        
        %-Load MAT-file
        %------------------------------------------------------------------
        try
            multicond = load(sess.multi{1});
        catch
            error('Cannot load %s',sess.multi{1});
        end
        
        %-Check structure content
        %------------------------------------------------------------------
        if ~all(isfield(multicond, {'names','onsets','durations'})) || ...
           ~iscell(multicond.names) || ...
           ~iscell(multicond.onsets) || ...
           ~iscell(multicond.durations) || ...
           ~isequal(numel(multicond.names), numel(multicond.onsets), ...
                    numel(multicond.durations))
            error(['Multiple conditions MAT-file ''%s'' is invalid:\n',...
           'File must contain names, onsets, and durations '...
           'cell arrays of equal length.\n'],sess.multi{1});
        end
    
        for j=1:numel(multicond.onsets)
            
            %-Mutiple Conditions: names, onsets and durations
            %--------------------------------------------------------------
            cond.name     = multicond.names{j};
            if isempty(cond.name)
                error('MultiCond file: sess %d cond %d has no name.',i,j);
            end
            cond.onset    = multicond.onsets{j};
            if isempty(cond.onset)
                error('MultiCond file: sess %d cond %d has no onset.',i,j);
            end
            cond.duration = multicond.durations{j};
            if isempty(cond.onset)
                error('MultiCond file: sess %d cond %d has no duration.',i,j);
            end
            
            %-Mutiple Conditions: Time Modulation
            %--------------------------------------------------------------
            if ~isfield(multicond,'tmod')
                cond.tmod = 0;
            else
                try
                    cond.tmod = multicond.tmod{j};
                catch
                    error('Error specifying time modulation.');
                end
            end

            %-Mutiple Conditions: Parametric Modulation
            %--------------------------------------------------------------
            cond.pmod = [];
            if isfield(multicond,'pmod')
                try
                    %-Check if a PM is defined for that condition
                    if (j <= numel(multicond.pmod)) && ...
                            ~isempty(multicond.pmod(j).name)
                        for ii = 1:numel(multicond.pmod(j).name)
                            cond.pmod(ii).name  = multicond.pmod(j).name{ii};
                            cond.pmod(ii).param = multicond.pmod(j).param{ii};
                            cond.pmod(ii).poly  = multicond.pmod(j).poly{ii};
                        end
                    end
                catch
                    warning('Error specifying parametric modulation.');
                    rethrow(lasterror);
                end
            end
            
            %-Mutiple Conditions: Orthogonalisation of Modulations
            %--------------------------------------------------------------
            if isfield(multicond,'orth') && (j <= numel(multicond.orth))
                cond.orth    = multicond.orth{j};
            else
                cond.orth    = true;
            end
            
            %-Append to singly-specified conditions
            %--------------------------------------------------------------
            sess.cond(end+1) = cond;
        end
    end

    %-Conditions
    %----------------------------------------------------------------------
    U = [];
    
    for j = 1:numel(sess.cond)
        
        %-Name, Onsets, Durations
        %------------------------------------------------------------------
        cond      = sess.cond(j);
        U(j).name = {cond.name};
        U(j).ons  = cond.onset(:);
        U(j).dur  = cond.duration(:);
        U(j).orth = cond.orth;
        if isempty(U(j).orth), U(j).orth = true; end
        if length(U(j).dur) == 1
            U(j).dur = repmat(U(j).dur,size(U(j).ons));
        elseif numel(U(j).dur) ~= numel(U(j).ons)
            error('Mismatch between number of onset and number of durations.');
        end

        %-Modulations
        %------------------------------------------------------------------
        P  = [];
        q1 = 0;
        %-Time Modulation
        switch job.timing.units
            case 'secs'
                sf    = 1 / 60;
            case 'scans'
                sf    = job.timing.RT / 60;
            otherwise
                error('Unknown unit "%s".',job.timing.units);
        end
        if cond.tmod > 0
            P(1).name = 'time';
            P(1).P    = U(j).ons * sf;
            P(1).h    = cond.tmod;
            q1        = 1;
        end
        %-Parametric Modulations
        if ~isempty(cond.pmod)
            for q = 1:numel(cond.pmod)
                q1 = q1 + 1;
                P(q1).name = cond.pmod(q).name;
                P(q1).P    = cond.pmod(q).param(:);
                P(q1).h    = cond.pmod(q).poly;
            end
        end
        %-None
        if isempty(P)
            P.name = 'none';
            P.h    = 0;
        end
        
        U(j).P = P;

    end

    SPM.Sess(i).U = U;

    %-User-specified regressors
    %----------------------------------------------------------------------
    C     = [];
    Cname = cell(1,numel(sess.regress));
    for q = 1:numel(sess.regress)
        Cname{q} = sess.regress(q).name;
        if numel(sess.regress(q).val(:)) ~= SPM.nscan(i)
            error('Length of regressor is not commensurate with data points.');
        end
        C        = [C, sess.regress(q).val(:)];
    end
    
    %-Multiple regressors (from a TXT/MAT-file)
    %----------------------------------------------------------------------
    if ~isempty(sess.multi_reg{1})
        for q=1:numel(sess.multi_reg)
            tmp   = load(sess.multi_reg{q});
            names = {};
            if isstruct(tmp) % .mat
                if isfield(tmp,'R')
                    R = tmp.R;
                    if isfield(tmp,'names')
                        names = tmp.names;
                    end
                elseif isfield(tmp,'xY')
                    R = tmp.xY.u;
                    names = {tmp.xY.name};
                elseif isfield(tmp,'PPI')
                    R    = [tmp.PPI.ppi tmp.PPI.Y tmp.PPI.P];
                    names = {...
                        ['PPI Interaction: ' tmp.PPI.name],...
                        ['Main Effect: ' tmp.PPI.xY.name ' BOLD'],...
                        'Main Effect: "psych" condition'};
                else
                    error(['Variable ''R'' not found in multiple ' ...
                        'regressor file ''%s''.'], sess.multi_reg{q});
                end
            elseif isnumeric(tmp) % .txt
                R     = tmp;
            end
            
            if size(R,1) ~= SPM.nscan(i)
                error('Length of regressor is not commensurate with data points.');
            end
            if ~isempty(names) && numel(names) ~= size(R,2)
                warning('Mismatch between number of regressors and their names.');
                names = {};
            end
            C  = [C, R];
            if numel(sess.multi_reg) == 1
                nb_mult = '';
            else
                nb_mult = sprintf('_%d',q);
            end
            for j=1:size(R,2)
                if isempty(names)
                    Cname{end+1} = sprintf('R%d%s',j,nb_mult);
                else
                    Cname{end+1} = sprintf('%s%s',names{j},nb_mult);
                end
            end
        end
    end
    
    SPM.Sess(i).C.C    = C;
    SPM.Sess(i).C.name = Cname;

end

%-Factorial design
%--------------------------------------------------------------------------
if ~isempty(job.fact)
    for i=1:numel(job.fact)
        SPM.factor(i).name   = job.fact(i).name;
        SPM.factor(i).levels = job.fact(i).levels;
    end
    if prod([SPM.factor.levels]) ~= numel(SPM.Sess(1).U)
        error('Factors do not match conditions');
    end
else
    SPM.factor = [];
end

%-Globals
%--------------------------------------------------------------------------
SPM.xGX.iGXcalc = job.global;

%-Masking threshold
%--------------------------------------------------------------------------
SPM.xM.gMT = job.mthresh;

%-High Pass filter
%--------------------------------------------------------------------------
for i = 1:numel(job.sess)
    SPM.xX.K(i).HParam = job.sess(i).hpf;
end

%-Autocorrelation
%--------------------------------------------------------------------------
SPM.xVi.form = job.cvi;


%-Setting up the GLM
%==========================================================================
SPM = spm_fmri_spm_ui(SPM);

%-Explicit mask
%--------------------------------------------------------------------------
if ~design_only
    if ~isempty(job.mask{1})
        SPM.xM.VM         = spm_data_hdr_read(job.mask{1});
        SPM.xM.xs.Masking = [SPM.xM.xs.Masking, '+explicit mask'];
    end
end

%-Save SPM.mat
%--------------------------------------------------------------------------
fprintf('%-40s: ','Saving SPM configuration')                           %-#
fmt = spm_get_defaults('mat.format');
s = whos('SPM');
if s.bytes > 2147483647, fmt = '-v7.3'; end
save('SPM.mat','SPM', fmt);
fprintf('%30s\n','...SPM.mat saved')                                    %-#

fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#

out.spmmat{1} = fullfile(pwd, 'SPM.mat');

%-Change back directory
%--------------------------------------------------------------------------
cd(cwd);
