function out = spm_run_factorial_design(job)
% SPM job execution function - factorial design specification
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - struct variable containing the path of the saved SPM.mat
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_run_factorial_design.m 7739 2019-12-02 14:00:18Z guillaume $

%--------------------------------------------------------------------------
% This function configures the design matrix (describing the general
% linear model), data specification, and other parameters necessary for
% the statistical analysis. These parameters are saved in a
% configuration file (SPM.mat) in the current directory, and are
% passed on to spm_spm.m (via the Estimate button) which estimates the
% design. Inference on these estimated parameters is then handled by the 
% SPM results section.
%
% This function, that sets up the necessary SPM structures, has been
% largely cannibalised from spm_spm_ui.m which we have retained for 
% developmental continuity.
%
% It has in common with spm_spm_ui.m its use of the I factor matrix,
% the H,C,B,G design matrix partitions, the sF,sCFI,CFIforms,sCC,CCforms,
% sGXcalc,sGloNorm,sGMsca option definition variables and use of the
% functions spm_DesMtx.m, spm_meanby.m and spm_get_vc.m.
%
% It departs from spm_spm_ui.m in that it does not use the design
% definition data structure D. Also, it uses the new SPM.factor field,
% which for the case of full factorial designs, is used to automatically
% generate contrasts testing for main effects and interactions.
%
% This function departs from spm_spm_ui.m in that it does not provide
% the same menu of design options (these were hardcoded in D). Instead it
% provides a number of options for simple designs (1) One-sample t-test,
% (2) Two-sample t-test, (3) Paired t-test and (4) Multiple regression.
% Two facilities are provided for specifying more complicated designs
% (5) Full-factorial and (6) Flexible-factorial. These should be able to
% specify all design options (and more) that were available in SPM2.
% For each of these design types one can additionally specify regressors
% using the `covariates' option.
%
% Options (5) and (6) differ in the efficiency (eg. number of key 
% strokes/button presses) with which a given design can be specified. For 
% example, one-way ANOVAs can be specified using either option, but (5) is 
% usually more efficient.
%
% Full-factorial designs
% ______________________
%
% This option is best used when you wish to test for all main effects and
% interactions in one-way, two-way or three-way ANOVAs.
%
% Design specification proceeds in 2 stages. Firstly, by creating new
% factors and specifying the number of levels and name for each.
% Nonsphericity, ANOVA-by-factor (for PET data) and scaling options (for
% PET data) can also be specified at this stage. Secondly, scans are 
% assigned separately to each cell. This accomodates unbalanced designs.
%
% For example, if you wish to test for a main effect in the population
% from which your subjects are drawn and have modelled that effect at the 
% first level using K basis functions (eg. K=3 informed basis functions) 
% you can use a one-way ANOVA with K-levels. Create a single factor with K
% levels and then assign the data to each cell eg. canonical, temporal
% derivative and dispersion derivative cells, where each cell is assigned
% scans from multiple subjects.
%
% SPM will automatically generate the contrasts necessary to test for all
% main effects and interactions.
%
% Flexible-factorial designs
% __________________________
%
% In this option the design matrix is created a block at a time. You can
% decide whether you wish each block to be a main effect or a (two-way)
% interaction.
%
% This option is best used for one-way, two-way or three-way ANOVAs but 
% where you do not wish to test for all possible main effects and
% interactions. This is perhaps most useful for PET where there is usually
% not enough data to test for all possible effects. Or for 3-way ANOVAs
% where you do not wish to test for all of the two-way interactions. A
% typical example here would be a group-by-drug-by-task analysis where,
% perhaps, only (i) group-by-drug or (ii) group-by-task interactions are
% of interest. In this case it is only necessary to have two-blocks in the
% design matrix - one for each interaction. The three-way interaction can
% then be tested for using a contrast that computes the difference between
% (i) and (ii).
%
% Design specification then proceeds in 3 stages. Firstly, factors are
% created and names specified for each. Nonsphericity, ANOVA-by-factor and
% scaling options can also be specified at this stage.
%
% Secondly, a list of scans is produced along with a factor matrix, I. 
% This is an nscan x 4 matrix of factor level indicators (see xX.I below).
% The first factor must be 'replication' but the other factors can be
% anything. Specification of I and the scan list can be achieved in one
% of two ways (a) the 'Specify All' option allows I to be typed in at the
% user interface or (more likely) loaded in from the matlab workspace. 
% All of the scans are then selected in one go. (b) the 'Subjects' option
% allows you to enter scans a subject at a time. The corresponding
% experimental conditions (ie. levels of factors) are entered at the same
% time. SPM will then create the factor matrix I. This style of interface
% is similar to that available in SPM2.
%
% Thirdly, the design matrix is built up a block at a time. Each block
% can be a main effect or a (two-way) interaction.
%
%--------------------------------------------------------------------------
%
% Variables saved in the SPM stucture:
%
% xY.VY         - nScan x 1 struct array of memory mapped images
%                 (see spm_vol for definition of the map structure)
% xX            - structure describing design matrix
% xX.I          - nScan x 4 matrix of factor level indicators
%                 I(n,i) is the level of factor i corresponding to image n
% xX.sF         - 1x4 cellstr containing the names of the four factors
%                 xX.sF{i} is the name of factor i
% xX.X          - design matrix
% xX.xVi        - correlation constraints for non-spericity correction
% xX.iH         - vector of H partition (condition effects) indices,
%                 identifying columns of X corresponding to H
% xX.iC         - vector of C partition (covariates of interest) indices
% xX.iB         - vector of B partition (block effects) indices
% xX.iG         - vector of G partition (nuisance variables) indices
% xX.name       - p x 1 cellstr of effect names corresponding to columns
%                 of the design matrix
%
% xC            - structure array of covariate details
% xC(i).rc      - raw (as entered) i-th covariate
% xC(i).rcname  - name of this covariate (string)
% xC(i).c       - covariate as appears in design matrix (after any scaling,
%                 centering of interactions)
% xC(i).cname   - cellstr containing names for effects corresponding to
%                 columns of xC(i).c
% xC(i).iCC     - covariate centering option
% xC(i).iCFI    - covariate by factor interaction option
% xC(i).type    - covariate type: 1=interest, 2=nuisance, 3=global
% xC(i).cols    - columns of design matrix corresponding to xC(i).c
% xC(i).descrip - cellstr containing a description of the covariate
%
% xGX           - structure describing global options and values
% xGX.iGXcalc   - global calculation option used
% xGX.sGXcalc   - string describing global calculation used
% xGX.rg        - raw globals (before scaling and such like)
% xGX.iGMsca    - grand mean scaling option
% xGX.sGMsca    - string describing grand mean scaling
% xGX.GM        - value for grand mean (/proportional) scaling
% xGX.gSF       - global scaling factor (applied to xGX.rg)
% xGX.iGC       - global covariate centering option
% xGX.sGC       - string describing global covariate centering option
% xGX.gc        - center for global covariate
% xGX.iGloNorm  - Global normalisation option
% xGX.sGloNorm  - string describing global normalisation option
%
% xM            - structure describing masking options
% xM.T          - Threshold masking value (-Inf=>None, real=>absolute,
%                 complex=>proportional (i.e. times global))
% xM.TH         - nScan x 1 vector of analysis thresholds, one per image
% xM.I          - Implicit masking (0=>none, 1=>implicit zero/NaN mask)
% xM.VM         - struct array of explicit mask images
%                 (empty if no explicit masks)
% xM.xs         - structure describing masking options
%                 (format is same as for xsDes described below)
%
% xsDes         - structure of strings describing the design:
%                 Fieldnames are essentially topic strings (use "_"'s for
%                 spaces), and the field values should be strings or cellstr's
%                 of information regarding that topic. spm_DesRep.m
%                 uses this structure to produce a printed description
%                 of the design, displaying the fieldnames (with "_"'s
%                 converted to spaces) in bold as topics, with
%                 the corresponding text to the right
%
%--------------------------------------------------------------------------

%-Output directory
%--------------------------------------------------------------------------
cwd = pwd;
d   = spm_file(job.dir{1},'cpath');
if ~exist(d,'dir')
    sts = mkdir(d);
    if ~sts, error('Error creating output directory "%s".',d); end
end
cd(d);

%-Ask about overwriting files from previous analyses...
%--------------------------------------------------------------------------
if exist(fullfile(job.dir{1},'SPM.mat'),'file')
    str = { 'Current directory contains existing SPM file:',...
        'Continuing will overwrite existing file!'};
    if spm_input(str,1,'bd','stop|continue',[1,0],1,mfilename)
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

%-Option definitions
%==========================================================================

%-Generic factor names
%--------------------------------------------------------------------------
sF = {'sF1','sF2','sF3','sF4'};

%-Covariate by factor interaction options
sCFI = {'<none>';...                                        %-1
    'with sF1';'with sF2';'with sF3';'with sF4';...         %-2:5
    'with sF2 (within sF4)';'with sF3 (within sF4)'};       %-6,7

%-DesMtx argument components for covariate by factor interaction options
% (Used for CFI's Covariate Centering (CC), GMscale & Global normalisation)
%--------------------------------------------------------------------------
CFIforms = {'[]',   'C',    '{}';...                        %-1
    'I(:,1)',       'FxC',  '{sF{1}}';...                   %-2
    'I(:,2)',       'FxC',  '{sF{2}}';...                   %-3
    'I(:,3)',       'FxC',  '{sF{3}}';...                   %-4
    'I(:,4)',       'FxC',  '{sF{4}}';...                   %-5
    'I(:,[4,2])',   'FxC',  '{sF{4},sF{2}}';...             %-6
    'I(:,[4,3])',   'FxC',  '{sF{4},sF{3}}' };              %-7

%-Centre (mean correction) options for covariates & globals            (CC)
% (options 9-12 are for centering of global when using AnCova GloNorm) (GC)
%--------------------------------------------------------------------------
sCC = {'around overall mean';...                            %-1
    'around sF1 means';...                                  %-2
    'around sF2 means';...                                  %-3
    'around sF3 means';...                                  %-4
    'around sF4 means';...                                  %-5
    'around sF2 (within sF4) means';...                     %-6
    'around sF3 (within sF4) means';...                     %-7
    '<no centering>';...                                    %-8
    'around user specified value';...                       %-9
    '(as implied by AnCova)';...                            %-10
    'GM';...                                                %-11
    '(redundant: not doing AnCova)'}';                      %-12

%-DesMtx I forms for covariate centering options
%--------------------------------------------------------------------------
CCforms = {'ones(nScan,1)',CFIforms{2:end,1},''}';

%-Global calculation options                                       (GXcalc)
%--------------------------------------------------------------------------
sGXcalc  = {'omit';...                                      %-1
    'user specified';...                                    %-2
    'mean voxel value (within per image fullmean/8 mask)'}; %-3

%-Global normalization options                                    (GloNorm)
%--------------------------------------------------------------------------
sGloNorm = {'AnCova';...                                    %-1
    'AnCova by sF1';...                                     %-2
    'AnCova by sF2';...                                     %-3
    'AnCova by sF3';...                                     %-4
    'AnCova by sF4';...                                     %-5
    'AnCova by sF2 (within sF4)';...                        %-6
    'AnCova by sF3 (within sF4)';...                        %-7
    'proportional scaling';...                              %-8
    '<no global normalisation>'};                           %-9


%-Grand mean scaling options                                        (GMsca)
% (NB: Grand mean scaling by subject is redundent for proportional scaling)
%--------------------------------------------------------------------------
sGMsca = {'scaling of overall grand mean';...               %-1
    'scaling of sF1 grand means';...                        %-2
    'scaling of sF2 grand means';...                        %-3
    'scaling of sF3 grand means';...                        %-4
    'scaling of sF4 grand means';...                        %-5
    'scaling of sF2 (within sF4) grand means';...           %-6
    'scaling of sF3 (within sF4) grand means';...           %-7
    '(implicit in PropSca global normalisation)';...        %-8
    '<no grand Mean scaling>'   };                          %-9

%-Conditions of no interest defaults
%--------------------------------------------------------------------------
B      = [];
Bnames = {};
factor = [];

switch char(fieldnames(job.des))
    
    %-One sample t-test
    %======================================================================
    case 't1',
        
        DesName = 'One sample t-test';

        P = job.des.t1.scans;
        n = length(P);
        I = (1:n)';
        I = [I,ones(n,3)];

        [H,Hnames] = spm_DesMtx(I(:,2),'-','mean');

        factor(1).name     = 'Group';
        factor(1).levels   = 1;
        factor(1).variance = 0;
        factor(1).dept     = 0;
        
    %-Two-sample t-test
    %======================================================================
    case 't2',
        
        DesName = 'Two-sample t-test';

        P  = job.des.t2.scans1;
        n1 = length(job.des.t2.scans1);
        P  = [P; job.des.t2.scans2];
        n2 = length(job.des.t2.scans2);

        I  = [(1:n1),(1:n2)]';
        I  = [I,[ones(n1,1);2*ones(n2,1)]];
        I  = [I,ones(n1+n2,2)];

        [H,Hnames] = spm_DesMtx(I(:,2),'-','Group');

        % Names and levels
        factor(1).name     = 'Group';
        factor(1).levels   = 2;

        % Ancova options
        factor(1).gmsca    = job.des.t2.gmsca;
        factor(1).ancova   = job.des.t2.ancova;

        % Nonsphericity options
        factor(1).variance = job.des.t2.variance;
        factor(1).dept     = job.des.t2.dept;

        if any([n1 n2]==1) && factor(1).variance == 1
            warning('Imposing equal variance between groups for 1 vs N comparison.');
            factor(1).variance = 0;
        end
        
    %-Paired t-test
    %======================================================================
    case 'pt',
        
        DesName = 'Paired t-test';

        Npairs  = length(job.des.pt.pair);
        P       = [];
        for p   = 1:Npairs
            P   = [P;job.des.pt.pair(p).scans];
        end

        I       = ones(Npairs*2,1);
        I(:,2)  = kron([1:Npairs]',ones(2,1));
        I(:,3)  = kron(ones(Npairs,1),[1 2]');
        I(:,4)  = I(:,1);

        [B,Bnames] = spm_DesMtx(I(:,2),'-','Subject');
        [H,Hnames] = spm_DesMtx(I(:,3),'-','Condition');

        % Names and levels
        factor(1).name     = 'Subject';
        factor(1).levels   = Npairs;
        factor(2).name     = 'Condition';
        factor(2).levels   = 2;

        % Ancova options
        factor(1).gmsca    = 0;
        factor(1).ancova   = 0;
        factor(2).gmsca    = job.des.pt.gmsca;
        factor(2).ancova   = job.des.pt.ancova;

        % Nonsphericity options
        factor(1).variance = 0;
        factor(1).dept     = 0;
        factor(2).variance = 0;
        factor(2).dept     = 0;

    %-Multiple regression
    %======================================================================
    case 'mreg',
        
        DesName = 'Multiple regression';

        P = job.des.mreg.scans;
        n = length(P);
        I = (1:n)';
        I = [I,ones(n,3)];

        % Names and levels
        factor(1).name     = '';
        factor(1).levels   = 1;

        % Nonsphericity options
        factor(1).variance = 0;
        factor(1).dept     = 0;

        if job.des.mreg.incint==0
            H = []; Hnames = '';
        else
            [H,Hnames] = spm_DesMtx(I(:,2),'-','mean');
        end

        for i=1:length(job.des.mreg.mcov)
            job.cov(end+1).c   = job.des.mreg.mcov(i).c;
            job.cov(end).cname = job.des.mreg.mcov(i).cname;
            job.cov(end).iCC   = job.des.mreg.mcov(i).iCC;
            job.cov(end).iCFI  = 1;
        end

    %-ANOVA
    %======================================================================
    case 'anova',
        
        DesName = 'ANOVA';
        
        job.des.anova.fact.name = 'Groups';
        
        % Automatically number cells 1 to levels, so user doesn't have to
        levels = length(job.des.anova.icell);
        job.des.anova.fact.levels = levels;
        for i=1:levels
            job.des.anova.icell(i).levels = i;
        end
        [I,P,H,Hnames] = spm_design_factorial(job.des.anova);

        
        % Names and levels
        factor(1).name     = 'Groups';
        factor(1).levels   = levels;
        
        % Ancova options
        factor(1).gmsca    = job.des.anova.gmsca;
        factor(1).ancova   = job.des.anova.ancova;
        
        % Nonsphericity options
        factor(1).variance = job.des.anova.variance;
        factor(1).dept     = job.des.anova.dept;
        
    %-ANOVA: within-subject
    %======================================================================
    case 'anovaw',
        
        DesName = 'ANOVA - within subject';
        
        anovaw  = job.des.anovaw;
        anovaw.fac(1).name      = 'Subject';
        anovaw.fac(1).dept      = 0;
        anovaw.fac(1).variance  = 0;
        anovaw.fac(1).gmsca     = 0;
        anovaw.fac(1).ancova    = 0;
        
        anovaw.fac(2).name      = 'Groups';
        anovaw.fac(2).dept      = job.des.anovaw.dept;
        anovaw.fac(2).variance  = job.des.anovaw.variance;
        anovaw.fac(2).gmsca     = job.des.anovaw.gmsca;
        anovaw.fac(2).ancova    = job.des.anovaw.ancova;
        
        anovaw.fsuball.fsubject = anovaw.fsubject;
        
        % Main effect of subject and group
        anovaw.maininters{1}.fmain.fnum = 1;
        anovaw.maininters{2}.fmain.fnum = 2;
        
        [I,P,job.cov] = spm_design_within_subject(anovaw,job.cov);
            
        [H,Hnames,B,Bnames] = spm_design_flexible(anovaw,I);
        
        factor           = anovaw.fac;
        factor(1).levels = length(job.des.anovaw.fsubject);
        
    %-Full Factorial Design
    %======================================================================
    case 'fd',
        
        DesName = 'Full factorial';

        [I,P,H,Hnames] = spm_design_factorial(job.des.fd);

        Nfactors = length(job.des.fd.fact);
        for i=1:Nfactors
            % Names and levels
            factor(i).name     = job.des.fd.fact(i).name;
            factor(i).levels   = job.des.fd.fact(i).levels;

            % Ancova options
            factor(i).gmsca    = job.des.fd.fact(i).gmsca;
            factor(i).ancova   = job.des.fd.fact(i).ancova;

            % Nonsphericity options
            factor(i).variance = job.des.fd.fact(i).variance;
            factor(i).dept     = job.des.fd.fact(i).dept;
        end
        
    %-Flexible factorial design
    %======================================================================
    case 'fblock',
        
        DesName = 'Flexible factorial';
        
        if isfield(job.des.fblock.fsuball,'fsubject')
            % Data has been entered subject by subject
            nf = length(job.des.fblock.fac);
            [I,P,job.cov] = spm_design_within_subject(job.des.fblock,job.cov);
        else
            % Specify all scans and factor matrix
            I = job.des.fblock.fsuball.specall.imatrix;
            [ns,nI] = size(I);
            % Pad out factorial matrix to cover the four canonical factors
            if nI < 4
                warning('Padding factor matrix to have four columns.');
                I = [I, ones(ns,4-nI)];
            end
            % Get number of factors
            nf = length(job.des.fblock.fac);
            P  = job.des.fblock.fsuball.specall.scans;
        end

        if isempty(job.des.fblock.maininters)
            warning('No main effects or interactions have been specified.');
        end
        [H,Hnames,B,Bnames] = spm_design_flexible(job.des.fblock,I);
        
        for i=1:nf
            % Names and levels
            factor(i).name     = job.des.fblock.fac(i).name;
            factor(i).levels   = length(unique(I(:,i+1)));

            % Ancova options
            factor(i).gmsca    = job.des.fblock.fac(i).gmsca;
            factor(i).ancova   = job.des.fblock.fac(i).ancova;

            % Nonsphericity options
            factor(i).variance = job.des.fblock.fac(i).variance;
            factor(i).dept     = job.des.fblock.fac(i).dept;
        end

end

nScan = size(I,1);


%-Covariate partition(s): interest (C) & nuisance (G) excluding global
%==========================================================================
dstr = {'covariate','nuisance variable'};
C  = []; Cnames = {};                 %-Covariate DesMtx partitions & names
G  = []; Gnames = {};

xC = [];                              %-Struct array to hold raw covariates

%-Multiple covariates
%--------------------------------------------------------------------------
for m=1:numel(job.multi_cov)
    for n=1:numel(job.multi_cov(m).files)
        tmp   = load(job.multi_cov(m).files{n});
        names = {};
        if isstruct(tmp) % .mat
            if isfield(tmp,'R')
                R = tmp.R;
                if isfield(tmp,'names')
                    names = tmp.names;
                end
            else
                error(['Variable ''R'' not found in multiple ' ...
                    'covariates file ''%s''.'], job.multi_cov(m).files{n});
            end
        elseif isnumeric(tmp) % .txt
            R     = tmp;
            % read names from first line if commented?
        end
        for j=1:size(R,2)
            job.cov(end+1).c   = R(:,j);
            if isempty(names)
                job.cov(end).cname = sprintf('R%d%s',j);
            else
                job.cov(end).cname = names{j};
            end
            job.cov(end).iCFI  = job.multi_cov(m).iCFI;
            job.cov(end).iCC   = job.multi_cov(m).iCC;
        end
    end
end


%-Covariates
%--------------------------------------------------------------------------
nc = length(job.cov);                 %-Number of covariates
for i=1:nc

    c      = job.cov(i).c;
    cname  = job.cov(i).cname;
    rc     = c;                       %-Save covariate value
    rcname = cname;                   %-Save covariate name
    if job.cov(i).iCFI==1
        iCFI=1;
    else
        % SPMs internal factor numbers are 1 higher than specified in user
        % interface as, internally, the first factor is always `replication'
        iCFI = job.cov(i).iCFI+1;
    end
    switch job.cov(i).iCC
        case 1
            iCC = 1;
        case {2,3,4}
            iCC = job.cov(i).iCC + 1;
        otherwise
            iCC = job.cov(i).iCC + 3;
    end

    %-Centre within factor levels as appropriate
    if any(iCC == (1:7))
        c = c - spm_meanby(c,eval(CCforms{iCC}));
    end

    %-Do any interaction (only for single covariate vectors)
    %----------------------------------------------------------------------
    if iCFI > 1                       %-(NB:iCFI=1 if size(c,2)>1)
        tI        = [eval(CFIforms{iCFI,1}),c];
        tConst    = CFIforms{iCFI,2};
        tFnames   = [eval(CFIforms{iCFI,3}),{cname}];
        [c,cname] = spm_DesMtx(tI,tConst,tFnames);
    elseif size(c,2)>1                %-Design matrix block
        [null,cname] = spm_DesMtx(c,'X',cname);
    else
        cname     = {cname};
    end

    %-Store raw covariate details in xC struct for reference
    %-Pack c into appropriate DesMtx partition
    %----------------------------------------------------------------------
    %-Construct description string for covariate
    str = {sprintf('%s',rcname)};
    if size(rc,2)>1, str = {sprintf('%s (block of %d covariates)',...
            str{:},size(rc,2))}; end
    if iCC < 8, str=[str;{['used centered ',sCC{iCC}]}]; end
    if iCFI> 1, str=[str;{['fitted as interaction ',sCFI{iCFI}]}]; end

    typ = 1;
    tmp = struct(...
        'rc',   rc,    'rcname', rcname,...
        'c',    c,     'cname',  {cname},...
        'iCC',  iCC,   'iCFI',   iCFI,...
        'type', typ,...
        'cols', [1:size(c,2)] + size([H,C],2) + size([B,G],2)*min(typ-1,1),...
        'descrip', {str});
    if isempty(xC), xC = tmp; else xC = [xC,tmp]; end
    C     = [C,c];
    Cnames = [Cnames; cname];

end
clear c tI tConst tFnames


%==========================================================================
% - C O N F I G U R E   D E S I G N -
%==========================================================================

%-Images & image info: Map Y image files and check consistency of
% dimensions and orientation / voxel size
%==========================================================================
fprintf('%-40s: ','Mapping files')                                      %-#
VY    = spm_data_hdr_read(char(P));

%-Check compatibility of images
%--------------------------------------------------------------------------
spm_check_orientations(VY);

fprintf('%30s\n','...done')                                             %-#


%-Global values, scaling and global normalisation
%==========================================================================
%-Compute global values
%--------------------------------------------------------------------------
switch char(fieldnames(job.globalc))
    case 'g_omit',
        iGXcalc = 1;
    case 'g_user',
        iGXcalc = 2;
    case 'g_mean',
        iGXcalc = 3;
end

switch job.globalm.glonorm
    case 1,
        iGloNorm = 9;
    case 2,
        iGloNorm = 8;
    case 3,
        iGloNorm = 1;
end
if factor(1).levels > 1
    % Override if factor-specific ANCOVA has been specified
    for i=1:length(factor)
        if factor(i).ancova
            iGloNorm=i+2;
        end
    end
end

%-Analysis threshold mask
%--------------------------------------------------------------------------
%-Work out available options:
% -Inf=>None, real=>absolute, complex=>proportional, (i.e. times global)
M_T = -Inf;
switch char(fieldnames(job.masking.tm)),
    case 'tma',
        % Absolute
        M_T = job.masking.tm.tma.athresh;
    case 'tmr',
        % Relative
        M_T = job.masking.tm.tmr.rthresh*sqrt(-1);
        % Need to force calculation of globals
        if iGXcalc~=2, iGXcalc=3; end
    case 'tm_none'
        % None
        M_T = -Inf;
end

if iGXcalc==1 && (any(iGloNorm == [1:5 8]) || ...
        (factor(1).levels > 1 && any([factor.gmsca])))
    % Over-ride omission of global calculation if we need it
    disp(' ');
    disp('SPM needs estimates of global activity.');
    disp('But you have specified to omit this computation.');
    disp('SPM has overridden this omission and will automatically compute ');
    disp('globals as the mean value of within brain voxels.');
    disp(' ');
    iGXcalc = 3;
end
sGXcalc = sGXcalc{iGXcalc};

switch iGXcalc,
    case 1
        %-Don't compute => no GMsca (iGMsca==9) or GloNorm (iGloNorm==9)
        g = [];
    case 2
        %-User specified globals
        g = job.globalc.g_user.global_uval;
    case 3
        %-Compute as mean voxel value (within per image fullmean/8 mask)
        g = zeros(nScan,1);
        fprintf('%-40s: %30s','Calculating globals',' ')                %-#
        for i = 1:nScan
            str = sprintf('%3d/%-3d',i,nScan);
            fprintf('%s%30s',repmat(sprintf('\b'),1,30),str)            %-#
            g(i) = spm_global(VY(i)); % FIXME % for meshes
        end
        fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')        %-#
    otherwise
        error('illegal iGXcalc')
end
rg = g;

fprintf('%-40s: ','Design configuration')                               %-#

%-Grand mean scaling options                                        (GMsca)
%--------------------------------------------------------------------------
if iGloNorm==8
    iGMsca=8;   %-grand mean scaling implicit in PropSca GloNorm
else
    switch char(fieldnames(job.globalm.gmsca))
        case 'gmsca_yes',
            iGMsca=1;
        case 'gmsca_no',
            iGMsca=9;
    end
    if factor(1).levels > 1
        % Over-ride if factor-specific scaling has been specified
        for i=1:numel(factor)
            if factor(i).gmsca
                iGMsca=i+2;
            end
        end
    end
end

%-Value for PropSca / GMsca                                            (GM)
%--------------------------------------------------------------------------
switch iGMsca,
    case 9                                %-Not scaling (GMsca or PropSca)
        GM = 0;                           %-Set GM to zero when not scaling
    case 1                                %-Ask user value of GM
        GM = job.globalm.gmsca.gmsca_yes.gmscv;
    otherwise
        if iGloNorm==8
            switch char(fieldnames(job.globalm.gmsca))
                case 'gmsca_yes',
                    % Proportionally scale to this value
                    GM = job.globalm.gmsca.gmsca_yes.gmscv;
                case 'gmsca_no',
                    GM = 50;
            end
        else
            % Grand mean scaling by factor eg. scans are scaled so that the
            % mean global value over each level of the factor is set to GM
            GM=50;
        end
end

%-If GM is zero then don't GMsca! or PropSca GloNorm
if GM==0,
    iGMsca=9;
    if iGloNorm==8,
        iGloNorm=9;
    end
end

%-Sort out description strings for GloNorm and GMsca
%--------------------------------------------------------------------------
sGloNorm = sGloNorm{iGloNorm};
sGMsca   = sGMsca{iGMsca};
if iGloNorm==8
    sGloNorm = sprintf('%s to %-4g',sGloNorm,GM);
elseif iGMsca<8
    sGMsca   = sprintf('%s to %-4g',sGMsca,GM);
end

%-Scaling: compute global scaling factors gSF required to implement
% proportional scaling global normalisation (PropSca) or grand mean
% scaling (GMsca), as specified by iGMsca (& iGloNorm)
%--------------------------------------------------------------------------
switch iGMsca,
    case 8
        %-Proportional scaling global normalisation
        if iGloNorm~=8, error('iGloNorm-iGMsca(8) mismatch for PropSca'), end
        gSF    = GM./g;
        g      = GM*ones(nScan,1);
    case {1,2,3,4,5,6,7}
        %-Grand mean scaling according to iGMsca
        if iGXcalc==1, error('Global calculation option is not appropriate.'), end
        gSF    = GM./spm_meanby(g,eval(CCforms{iGMsca}));
        g      = g.*gSF;
    case 9
        %-No grand mean scaling
        gSF    = ones(nScan,1);
    otherwise
        error('illegal iGMsca')
end

%-Apply gSF to memory-mapped scalefactors to implement scaling
%--------------------------------------------------------------------------
for i = 1:nScan
    VY(i).pinfo(1:2,:) = VY(i).pinfo(1:2,:)*gSF(i); % FIXME % for meshes
end

%-Global centering (for AnCova GloNorm)                                (GC)
%-If not doing AnCova then GC is irrelevant
%--------------------------------------------------------------------------
if ~any(iGloNorm == [1:7])
    iGC = 12;
    gc  = [];
else
    iGC = 10;
    gc = 0;
end

%-AnCova: Construct global nuisance covariates partition (if AnCova)
%--------------------------------------------------------------------------
if any(iGloNorm == [1:7])

    %-Centre global covariate as requested
    %----------------------------------------------------------------------
    switch iGC, case {1,2,3,4,5,6,7}    %-Standard sCC options
        gc = spm_meanby(g,eval(CCforms{iGC}));
        case 8                  %-No centering
            gc = 0;
        case 9                  %-User specified centre
            %-gc set above
        case 10                 %-As implied by AnCova option
            gc = spm_meanby(g,eval(CCforms{iGloNorm}));
        case 11                 %-Around GM
            gc = GM;
        otherwise               %-unknown iGC
            error('unexpected iGC value')
    end

    %-AnCova - add scaled centred global to DesMtx `G' partition
    %----------------------------------------------------------------------
    rcname     = 'global';
    tI         = [eval(CFIforms{iGloNorm,1}),g - gc];
    tConst     = CFIforms{iGloNorm,2};
    tFnames    = [eval(CFIforms{iGloNorm,3}),{rcname}];
    [f,gnames]  = spm_DesMtx(tI,tConst,tFnames);
    clear tI tConst tFnames

    %-Save GX info in xC struct for reference
    %----------------------------------------------------------------------
    str     = {sprintf('%s: %s',dstr{2},rcname)};
    if any(iGMsca==[1:7]), str=[str;{['(after ',sGMsca,')']}]; end
    if iGC ~= 8, str=[str;{['used centered ',sCC{iGC}]}]; end
    if iGloNorm > 1
        str=[str;{['fitted as interaction ',sCFI{iGloNorm}]}];
    end
    tmp  = struct(  'rc',rg.*gSF,       'rcname',rcname,...
        'c',f,          'cname' ,{gnames},...
        'iCC',iGC,      'iCFI'  ,iGloNorm,...
        'type',         3,...
        'cols',[1:size(f,2)] + size([H C B G],2),...
        'descrip',      {str}       );

    G = [G,f]; Gnames = [Gnames; gnames];
    if isempty(xC), xC = tmp; else xC = [xC,tmp]; end

elseif iGloNorm==8 || iGXcalc>1

    %-Globals calculated, but not AnCova: Make a note of globals
    %----------------------------------------------------------------------
    if iGloNorm==8
        str = { 'global values: (used for proportional scaling)';...
            '("raw" unscaled globals shown)'};
    elseif isfinite(M_T) && ~isreal(M_T)
        str = { 'global values: (used to compute analysis threshold)'};
    else
        str = { 'global values: (computed but not used)'};
    end

    rcname ='global';
    tmp     = struct('rc',rg,    'rcname',rcname,...
        'c',{[]},   'cname' ,{{}},...
        'iCC',0,    'iCFI'  ,0,...
        'type',     3,...
        'cols',     {[]},...
        'descrip',  {str}           );

    if isempty(xC), xC = tmp; else xC = [xC,tmp]; end
end

%-Save info on global calculation in xGX structure
%--------------------------------------------------------------------------
xGX = struct(...
    'iGXcalc', iGXcalc,  'sGXcalc', sGXcalc,  'rg',rg,...
    'iGMsca',  iGMsca,   'sGMsca',  sGMsca,   'GM',GM,    'gSF',gSF,...
    'iGC',     iGC,      'sGC',     sCC{iGC}, 'gc',gc,...
    'iGloNorm',iGloNorm, 'sGloNorm',sGloNorm);

%-Make a description string
%--------------------------------------------------------------------------
if isinf(M_T)
    xsM.Analysis_threshold = 'None (-Inf)';
elseif isreal(M_T)
    xsM.Analysis_threshold = sprintf('images thresholded at %6g',M_T);
else
    xsM.Analysis_threshold = sprintf(['images thresholded at %6g ',...
        'times global'],imag(M_T));
end

%-Construct masking information structure and compute actual analysis
% threshold using scaled globals (rg.*gSF)
%--------------------------------------------------------------------------
if isreal(M_T),
    M_TH = M_T  * ones(nScan,1);    %-NB: -Inf is real
else
    M_TH = imag(M_T) * (rg.*gSF);
end

%-Implicit masking: Ignore zero voxels in low data-types?
%--------------------------------------------------------------------------
% (Implicit mask is NaN in higher data-types.)
if ~spm_type(VY(1).dt(1),'nanrep')
    M_I = job.masking.im;  % Implicit mask ?
    if M_I
        xsM.Implicit_masking = 'Yes: zero''s treated as missing';
    else
        xsM.Implicit_masking = 'No';
    end
else
    M_I = 1;
    xsM.Implicit_masking = 'Yes: NaN''s treated as missing';
end

%-Explicit masking
%--------------------------------------------------------------------------
if isempty(job.masking.em{:})
    VM = [];
    xsM.Explicit_masking = 'No';
else
    VM = spm_data_hdr_read(char(job.masking.em));
    xsM.Explicit_masking = 'Yes';
end

xM     = struct('T',M_T, 'TH',M_TH, 'I',M_I, 'VM',{VM}, 'xs',xsM);


%-Construct full design matrix (X), parameter names and structure (xX)
%==========================================================================
X      = [H C B G];
tmp    = cumsum([size(H,2), size(C,2), size(B,2), size(G,2)]);
xX     = struct(...
    'X',        X,...
    'iH',       [1:size(H,2)],...
    'iC',       [1:size(C,2)] + tmp(1),...
    'iB',       [1:size(B,2)] + tmp(2),...
    'iG',       [1:size(G,2)] + tmp(3),...
    'name',     {[Hnames; Cnames; Bnames; Gnames]},...
    'I',        I,...
    'sF',       {sF});


%-Design description (an nx2 cellstr) - for saving and display
%==========================================================================
tmp = {sprintf('%d condition, +%d covariate, +%d block, +%d nuisance',...
    size(H,2),size(C,2),size(B,2),size(G,2));...
    sprintf('%d total, having %d degrees of freedom',...
    size(X,2),rank(X));...
    sprintf('leaving %d degrees of freedom from %d images',...
    size(X,1)-rank(X),size(X,1))};
xsDes = struct('Design',    {DesName},...
    'Global_calculation',   {sGXcalc},...
    'Grand_mean_scaling',   {sGMsca},...
    'Global_normalisation', {sGloNorm},...
    'Parameters',           {tmp});

fprintf('%30s\n','...done')                                             %-#


%-Generate error covariance components (non-sphericity)
%==========================================================================
Vi          = spm_get_vc(I, factor);


%-Assemble SPM structure
%==========================================================================
SPM.xY.P    = P;            % filenames
SPM.xY.VY   = VY;           % mapped data
SPM.nscan   = size(xX.X,1); % scan number
SPM.xX      = xX;           % design structure
SPM.xC      = xC;           % covariate structure
SPM.xGX     = xGX;          % global structure
SPM.xM      = xM;           % mask structure
SPM.xsDes   = xsDes;        % description
SPM.xVi.I   = I;            % factor matrix
if numel(Vi) == 1
    SPM.xVi.V  = Vi{1};     % non-sphericity matrix
else
    SPM.xVi.Vi = Vi;        % non-sphericity variance components
end

%-Automatic contrast generation for 'Full factorial'
%--------------------------------------------------------------------------
if strcmp(char(fieldnames(job.des)),'fd') && job.des.fd.contrasts
    SPM.factor = factor;
end

%-Save SPM.mat and set output argument
%--------------------------------------------------------------------------
fprintf('%-40s: ','Saving SPM configuration')                           %-#
fmt = spm_get_defaults('mat.format');
s = whos('SPM');
if s.bytes > 2147483647, fmt = '-v7.3'; end
save('SPM.mat', 'SPM', fmt);
fprintf('%30s\n','...SPM.mat saved')                                    %-#

out.spmmat{1} = fullfile(pwd, 'SPM.mat');


%-Display Design report
%==========================================================================
if ~spm('CmdLine') && ~isempty(spm_figure('FindWin','Graphics'))
    fprintf('%-40s: ','Design reporting')                               %-#
    fname     = cat(1,{SPM.xY.VY.fname}');
    spm_DesRep('DesMtx',SPM.xX,fname,SPM.xsDes)
    fprintf('%30s\n','...done')                                         %-#
end

fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#

%-Change back directory
%--------------------------------------------------------------------------
cd(cwd);
