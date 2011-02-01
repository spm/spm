function out = spm_run_fmri_spec(job)
% Set up the design matrix and run a design.
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_fmri_spec.m 4185 2011-02-01 18:46:18Z guillaume $


original_dir = pwd;
my_cd(job.dir);

%-Ask about overwriting files from previous analyses
%--------------------------------------------------------------------------
if exist(fullfile(job.dir{1},'SPM.mat'),'file')
    str = { 'Current directory contains existing SPM file:',...
        'Continuing will overwrite existing file!'};
    if spm_input(str,1,'bd','stop|continue',[1,0],1,mfilename);
        fprintf('%-40s: %30s\n\n',...
            'Abort...   (existing SPM file)',spm('time'));
        return
    end
end

% If we've gotten to this point we're committed to overwriting files.
% Delete them so we don't get stuck in spm_spm
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

% Variables
%--------------------------------------------------------------------------
SPM.xY.RT = job.timing.RT;
SPM.xY.P = [];

% Slice timing
%--------------------------------------------------------------------------
% The following lines have the side effect of modifying the global
% defaults variable. This is necessary to pass job.timing.fmri_t to
% spm_hrf.m. The original values are saved here and restored at the end
% of this function, after the design has been specified. The original
% values may not be restored if this function crashes.
olddefs.stats.fmri.fmri_t  = spm_get_defaults('stats.fmri.fmri_t');
olddefs.stats.fmri.fmri_t0 = spm_get_defaults('stats.fmri.fmri_t0');
spm_get_defaults('stats.fmri.t',  job.timing.fmri_t);
spm_get_defaults('stats.fmri.t0', job.timing.fmri_t0);

% Basis function variables
%--------------------------------------------------------------------------
SPM.xBF.UNITS = job.timing.units;
SPM.xBF.dt    = job.timing.RT/job.timing.fmri_t;
SPM.xBF.T     = job.timing.fmri_t;
SPM.xBF.T0    = job.timing.fmri_t0;

% Basis functions
%--------------------------------------------------------------------------
if strcmp(fieldnames(job.bases),'hrf')
    if all(job.bases.hrf.derivs == [0 0])
        SPM.xBF.name = 'hrf';
    elseif all(job.bases.hrf.derivs == [1 0])
        SPM.xBF.name = 'hrf (with time derivative)';
    elseif all(job.bases.hrf.derivs == [1 1])
        SPM.xBF.name = 'hrf (with time and dispersion derivatives)';
    else
        error('Unrecognized hrf derivative choices.')
    end
else
    nambase = fieldnames(job.bases);
    if ischar(nambase)
        nam=nambase;
    else
        nam=nambase{1};
    end
    switch nam,
        case 'fourier',
            SPM.xBF.name = 'Fourier set';
        case 'fourier_han',
            SPM.xBF.name = 'Fourier set (Hanning)';
        case 'gamma',
            SPM.xBF.name = 'Gamma functions';
        case 'fir',
            SPM.xBF.name = 'Finite Impulse Response';
        otherwise
            error('Unrecognized hrf derivative choices.')
    end
    SPM.xBF.length = job.bases.(nam).length;
    SPM.xBF.order  = job.bases.(nam).order;
end
SPM.xBF          = spm_get_bf(SPM.xBF);
if isempty(job.sess),
    SPM.xBF.Volterra = false;
else
    SPM.xBF.Volterra = job.volt;
end

for i = 1:numel(job.sess),
    sess = job.sess(i);

    % Image filenames
    %----------------------------------------------------------------------
    SPM.nscan(i) = numel(sess.scans);
    SPM.xY.P     = strvcat(SPM.xY.P,sess.scans{:});
    U = [];

    % Augment the singly-specified conditions with the multiple
    % conditions specified in a .mat file provided by the user
    %----------------------------------------------------------------------
    if ~isempty(sess.multi{1})
        try
            multicond = load(sess.multi{1});
        catch
            error('Cannot load %s',sess.multi{1});
        end
        if ~(isfield(multicond,'names')&&isfield(multicond,'onsets')&&...
                isfield(multicond,'durations')) || ...
            ~all([numel(multicond.names),numel(multicond.onsets), ...
                numel(multicond.durations)]==numel(multicond.names))
            error(['Multiple conditions MAT-file ''%s'' is invalid.\n',...
           'File must contain names, onsets, and durations '...
           'cell arrays of equal length.\n'],sess.multi{1});
        end
    
        %-contains three cell arrays: names, onsets and durations
        for j=1:length(multicond.onsets)
            cond.name     = multicond.names{j};
            cond.onset    = multicond.onsets{j};
            cond.duration = multicond.durations{j};
            
            % Mutiple Conditions Time Modulation
            %--------------------------------------------------------------
            % initialise the variable.
            cond.tmod = 0;
            if isfield(multicond,'tmod');
                try
                    cond.tmod = multicond.tmod{j};
                catch
                    error('Error specifying time modulation.');
                end
            end

            % Mutiple Conditions Parametric Modulation
            %--------------------------------------------------------------
            % initialise the parametric modulation variable.
            cond.pmod = [];
            if isfield(multicond,'pmod')
                % only access existing modulators
                try
                    % check if there is a parametric modulator. this allows
                    % pmod structures with fewer entries than conditions.
                    % then check whether any cells are filled in.
                    if (j <= numel(multicond.pmod)) && ...
                            ~isempty(multicond.pmod(j).name)

                        % we assume that the number of cells in each
                        % field of pmod is the same (or should be).
                        for ii = 1:numel(multicond.pmod(j).name)
                            cond.pmod(ii).name  = multicond.pmod(j).name{ii};
                            cond.pmod(ii).param = multicond.pmod(j).param{ii};
                            cond.pmod(ii).poly  = multicond.pmod(j).poly{ii};
                        end
                    end;
                catch
                    error('Error specifying parametric modulation.');
                end
            end
            sess.cond(end+1) = cond;
        end
    end

    % Configure the input structure array
    %----------------------------------------------------------------------
    for j = 1:length(sess.cond),
        cond      = sess.cond(j);
        U(j).name = {cond.name};
        U(j).ons  = cond.onset(:);
        U(j).dur  = cond.duration(:);
        if length(U(j).dur) == 1
            U(j).dur    = U(j).dur*ones(size(U(j).ons));
        elseif length(U(j).dur) ~= length(U(j).ons)
            error('Mismatch between number of onset and number of durations.')
        end

        P  = [];
        q1 = 0;
        if cond.tmod>0,
            % time effects
            P(1).name = 'time';
            P(1).P    = U(j).ons*job.timing.RT/60;
            P(1).h    = cond.tmod;
            q1        = 1;
        end;
        if ~isempty(cond.pmod)
            for q = 1:numel(cond.pmod),
                % Parametric effects
                q1 = q1 + 1;
                P(q1).name = cond.pmod(q).name;
                P(q1).P    = cond.pmod(q).param(:);
                P(q1).h    = cond.pmod(q).poly;
            end;
        end
        if isempty(P)
            P.name = 'none';
            P.h    = 0;
        end
        U(j).P = P;

    end

    SPM.Sess(i).U = U;


    % User specified regressors
    %----------------------------------------------------------------------
    C     = [];
    Cname = cell(1,numel(sess.regress));
    for q = 1:numel(sess.regress),
        Cname{q} = sess.regress(q).name;
        C        = [C, sess.regress(q).val(:)];
    end

    % Augment the singly-specified regressors with the multiple regressors
    % specified in the regressors.mat file
    %----------------------------------------------------------------------
    if ~strcmp(sess.multi_reg,'')
        tmp = load(char(sess.multi_reg{:}));
        if isstruct(tmp) && isfield(tmp,'R')
            R = tmp.R;
        elseif isnumeric(tmp)
            % load from e.g. text file
            R = tmp;
        else
            warning('Can''t load user specified regressors in %s', ...
                char(sess.multi_reg{:}));
            R = [];
        end

        C  = [C, R];
        nr = size(R,2);
        nq = length(Cname);
        for inr=1:nr,
            Cname{inr+nq} = ['R',int2str(inr)];
        end
    end
    SPM.Sess(i).C.C    = C;
    SPM.Sess(i).C.name = Cname;

end

% Factorial design
%--------------------------------------------------------------------------
if isfield(job,'fact')
    if ~isempty(job.fact)
        NC=length(SPM.Sess(1).U); % Number of conditions
        CheckNC=1;
        for i=1:length(job.fact)
            SPM.factor(i).name=job.fact(i).name;
            SPM.factor(i).levels=job.fact(i).levels;
            CheckNC=CheckNC*SPM.factor(i).levels;
        end
        if ~(CheckNC==NC)
            disp('Error in fmri_spec job: factors do not match conditions');
            return
        end
    end
else
    SPM.factor=[];
end

% Globals
%--------------------------------------------------------------------------
SPM.xGX.iGXcalc = job.global;
SPM.xGX.sGXcalc = 'mean voxel value';
SPM.xGX.sGMsca  = 'session specific';

% High Pass filter
%--------------------------------------------------------------------------
for i = 1:numel(job.sess),
    SPM.xX.K(i).HParam = job.sess(i).hpf;
end

% Autocorrelation
%--------------------------------------------------------------------------
SPM.xVi.form = job.cvi;

% Let SPM configure the design
%--------------------------------------------------------------------------
SPM = spm_fmri_spm_ui(SPM);

if ~isempty(job.mask)&&~isempty(job.mask{1})
    SPM.xM.VM         = spm_vol(job.mask{:});
    SPM.xM.xs.Masking = [SPM.xM.xs.Masking, '+explicit mask'];
end

%-Save SPM.mat
%--------------------------------------------------------------------------
fprintf('%-40s: ','Saving SPM configuration')                           %-#
if spm_check_version('matlab','7') >= 0
    save('SPM.mat','-V6','SPM');
else
    save('SPM.mat','SPM');
end;

fprintf('%30s\n','...SPM.mat saved')                                    %-#

out.spmmat{1} = fullfile(pwd, 'SPM.mat');
my_cd(original_dir); % Change back dir
spm_get_defaults('stats.fmri.fmri_t',olddefs.stats.fmri.fmri_t); % Restore old timing
spm_get_defaults('stats.fmri.fmri_t0',olddefs.stats.fmri.fmri_t0); % parameters
fprintf('Done\n')
return

%==========================================================================
function my_cd(jobDir)
if ~isempty(jobDir)
    try
        cd(char(jobDir));
    catch
        error('Failed to change directory. Aborting run.')
    end
end
