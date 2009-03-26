function out = spm_run_factorial_design(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_factorial_design.m 2964 2009-03-26 16:18:28Z guillaume $


global defaults
if isempty(defaults)
    spm_defaults;
end

original_dir = pwd;
cd(job.dir{1});

%-Ask about overwriting files from previous analyses...
%-------------------------------------------------------------------
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
%------------------------------------------------------------------------
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
%-------------------------------------------------------------------
%-Generic factor names
sF = {'sF1','sF2','sF3','sF4'};

%-Covariate by factor interaction options
sCFI = {'<none>';...                            %-1
    'with sF1';'with sF2';'with sF3';'with sF4';...         %-2:5
    'with sF2 (within sF4)';'with sF3 (within sF4)'};       %-6,7

%-DesMtx argument components for covariate by factor interaction options
% (Used for CFI's Covariate Centering (CC), GMscale & Global normalisation)
CFIforms = {    '[]',       'C',    '{}';...            %-1
    'I(:,1)',       'FxC',  '{sF{1}}';...           %-2
    'I(:,2)',       'FxC',  '{sF{2}}';...           %-3
    'I(:,3)',       'FxC',  '{sF{3}}';...           %-4
    'I(:,4)',       'FxC',  '{sF{4}}';...           %-5
    'I(:,[4,2])',   'FxC',  '{sF{4},sF{2}}';... %-6
    'I(:,[4,3])',   'FxC',  '{sF{4},sF{3}}' };  %-7

%-Centre (mean correction) options for covariates & globals            (CC)
% (options 9-12 are for centering of global when using AnCova GloNorm) (GC)
sCC = {     'around overall mean';...               %-1
    'around sF1 means';...                  %-2
    'around sF2 means';...                  %-3
    'around sF3 means';...                  %-4
    'around sF4 means';...                  %-5
    'around sF2 (within sF4) means';...         %-6
    'around sF3 (within sF4) means';...         %-7
    '<no centering>';...                    %-8
    'around user specified value';...           %-9
    '(as implied by AnCova)';...                %-10
    'GM';...                        %-11
    '(redundant: not doing AnCova)'}';          %-12
%-DesMtx I forms for covariate centering options
CCforms = {'ones(nScan,1)',CFIforms{2:end,1},''}';

%-Global calculation options                                       (GXcalc)
sGXcalc  = {    'omit';...                      %-1
    'user specified';...                    %-2
    'mean voxel value (within per image fullmean/8 mask)'}; %-3


%-Global normalization options  (GloNorm)
sGloNorm = {    'AnCova';...                        %-1
    'AnCova by sF1';...                 %-2
    'AnCova by sF2';...                 %-3
    'AnCova by sF3';...                 %-4
    'AnCova by sF4';...                 %-5
    'AnCova by sF2 (within sF4)';...            %-6
    'AnCova by sF3 (within sF4)';...            %-7
    'proportional scaling';...              %-8
    '<no global normalisation>'};               %-9


%-Grand mean scaling options                                        (GMsca)
sGMsca = {  'scaling of overall grand mean';...         %-1
    'scaling of sF1 grand means';...            %-2
    'scaling of sF2 grand means';...            %-3
    'scaling of sF3 grand means';...            %-4
    'scaling of sF4 grand means';...            %-5
    'scaling of sF2 (within sF4) grand means';...       %-6
    'scaling of sF3 (within sF4) grand means';...       %-7
    '(implicit in PropSca global normalisation)';...    %-8
    '<no grand Mean scaling>'   };          %-9
%-NB: Grand mean scaling by subject is redundent for proportional scaling

% Conditions of no interest defaults
B=[];
Bnames={};

switch strvcat(fieldnames(job.des)),
    case 't1',
        % One sample t-test
        DesName='One sample t-test';

        P=job.des.t1.scans;
        n=length(P);
        I=[1:n]';
        I=[I,ones(n,3)];

        [H,Hnames]=spm_DesMtx(I(:,2),'-','mean');

        SPM.factor(1).name='Group';
        SPM.factor(1).levels=1;
        SPM.factor(1).variance=0;
        SPM.factor(1).dept=0;
    case 't2',
        % Two-sample t-test
        DesName='Two-sample t-test';

        P=job.des.t2.scans1;
        n1=length(job.des.t2.scans1);
        P=[P;job.des.t2.scans2];
        n2=length(job.des.t2.scans2);

        I=[];
        I=[1:n1]';
        I=[I;[1:n2]'];
        I=[I,[ones(n1,1);2*ones(n2,1)]];
        I=[I,ones(n1+n2,2)];

        [H,Hnames]=spm_DesMtx(I(:,2),'-','Group');

        % Names and levels
        SPM.factor(1).name='Group';
        SPM.factor(1).levels=2;

        % Ancova options
        SPM.factor(1).gmsca=job.des.t2.gmsca;
        SPM.factor(1).ancova=job.des.t2.ancova;

        % Nonsphericity options
        SPM.factor(1).variance=job.des.t2.variance;
        SPM.factor(1).dept=job.des.t2.dept;

    case 'pt',
        % Paired t-test
        DesName='Paired t-test';

        Npairs=length(job.des.pt.pair);
        P=[];
        for p=1:Npairs,
            P=[P;job.des.pt.pair(p).scans];
        end

        I=ones(Npairs*2,1);
        I(:,2)=kron([1:Npairs]',ones(2,1));
        I(:,3)=kron(ones(Npairs,1),[1 2]');
        I(:,4)=I(:,1);

        [H,Hnames]=spm_DesMtx(I(:,2),'-','Subject');
        [B,Bnames]=spm_DesMtx(I(:,3),'-','Condition');

        % Names and levels
        SPM.factor(1).name='Subject';
        SPM.factor(1).levels=Npairs;
        SPM.factor(2).name='Condition';
        SPM.factor(2).levels=2;

        % Ancova options
        SPM.factor(1).gmsca=0;
        SPM.factor(1).ancova=0;
        SPM.factor(2).gmsca=job.des.pt.gmsca;
        SPM.factor(2).ancova=job.des.pt.ancova;

        % Nonsphericity options
        SPM.factor(1).variance=0;
        SPM.factor(1).dept=0;
        SPM.factor(2).variance=job.des.pt.variance;
        SPM.factor(2).dept=job.des.pt.dept;

    case 'mreg',
        % Multiple regression
        DesName='Multiple regression';

        P=job.des.mreg.scans;
        n=length(P);
        I=[1:n]';
        I=[I,ones(n,3)];

        % Names and levels
        SPM.factor(1).name='';
        SPM.factor(1).levels=1;

        % Nonsphericity options
        SPM.factor(1).variance=0;
        SPM.factor(1).dept=0;

        H=[];Hnames=[];
        if job.des.mreg.incint==0
            B = []; Bnames = '';
        else
            [B,Bnames] = spm_DesMtx(I(:,2),'-','mean');
        end

        for i=1:length(job.des.mreg.mcov)
            job.cov(end+1).c   = job.des.mreg.mcov(i).c;
            job.cov(end).cname = job.des.mreg.mcov(i).cname;
            job.cov(end).iCC   = job.des.mreg.mcov(i).iCC;
            job.cov(end).iCFI  = 1;
        end

    case 'fd',
        % Full Factorial Design
        DesName='Full factorial';

        [I,P,H,Hnames] = spm_set_factorial_design (job);

        Nfactors=length(job.des.fd.fact);
        for i=1:Nfactors,
            % Names and levels
            SPM.factor(i).name=job.des.fd.fact(i).name;
            SPM.factor(i).levels=job.des.fd.fact(i).levels;

            % Ancova options
            SPM.factor(i).gmsca=job.des.fd.fact(i).gmsca;
            SPM.factor(i).ancova=job.des.fd.fact(i).ancova;

            % Nonsphericity options
            SPM.factor(i).variance=job.des.fd.fact(i).variance;
            SPM.factor(i).dept=job.des.fd.fact(i).dept;

        end

    case 'fblock',
        % Flexible factorial design
        DesName='Flexible factorial';

        if isfield(job.des.fblock.fsuball,'fsubject')
            nsub=length(job.des.fblock.fsuball.fsubject);
            % Specify design subject-by-subject
            P=[];I=[];
            subj=[];
            for s=1:nsub,
                P  = [P; job.des.fblock.fsuball.fsubject(s).scans];
                ns = length(job.des.fblock.fsuball.fsubject(s).scans);
                cc = job.des.fblock.fsuball.fsubject(s).conds;

                [ccr,ccc] = size(cc);
                if ~(ccr==ns) && ~(ccc==ns)
                    disp(sprintf('Error for subject %d: conditions not specified for each scan',s));
                    return
                elseif ~(ccr==ccc) && (ccc==ns)
                    %warning('spm:transposingConditions',['Condition matrix ',...
                    %    'appears to be transposed. Transposing back to fix.\n',...
                    %    'Alert developers if it is not actually transposed.'])
                    cc=cc';
                end
                subj=[subj;s*ones(ns,1)];
                % get real replications within each subject cell
                [unused cci  ccj] = unique(cc,'rows');
                repl = zeros(ns, 1);
                for k=1:max(ccj)
                    repl(ccj==k) = 1:sum(ccj==k);
                end;
                I = [I; [repl cc]];
            end

            nf=length(job.des.fblock.fac);
            subject_factor=0;
            for i=1:nf,
                if strcmpi(job.des.fblock.fac(i).name,'repl')
                    % Copy `replications' column to create explicit `replications' factor
                    nI=I(:,1:i);
                    nI=[nI,I(:,1)];
                    nI=[nI,I(:,i+1:end)];
                    I=nI;
                end
                if strcmpi(job.des.fblock.fac(i).name,'subject')
                    % Create explicit `subject' factor
                    nI=I(:,1:i);
                    nI=[nI,subj];
                    nI=[nI,I(:,i+1:end)];
                    I=nI;
                    subject_factor=1;
                end

            end

            % Re-order scans conditions and covariates into standard format
            % This is to ensure compatibility with how variance components are created
            if subject_factor
                U=unique(I(:,2:nf+1),'rows');
                Un=length(U);
                Uc=zeros(Un,1);
                r=1;rj=[];
                for k=1:Un,
                    for j=1:size(I,1),
                        match=sum(I(j,2:nf+1)==U(k,:))==nf;
                        if match
                            Uc(k)=Uc(k)+1;
                            Ir(r,:)=[Uc(k),I(j,2:end)];
                            r=r+1;
                            rj=[rj;j];
                        end
                    end
                end
                P=P(rj); % -scans
                I=Ir;    % -conditions
                for k=1:numel(job.cov) % -covariates
                    job.cov(k).c = job.cov(k).c(rj);
                end;
            end

        else % specify all scans and factor matrix
            [ns,nc]=size(job.des.fblock.fsuball.specall.imatrix);
            if ~(nc==4)
                disp('Error: factor matrix must have four columns');
                return
            end
            I=job.des.fblock.fsuball.specall.imatrix;

            % Get number of factors
            nf=length(job.des.fblock.fac);
            %        nf=0;
            %        for i=1:4,
            %            if length(unique(I(:,i)))>1
            %                nf=nf+1;
            %            end
            %        end

            P=job.des.fblock.fsuball.specall.scans;
        end

        % Pad out factorial matrix to cover the four canonical factors
        [ns,nI]=size(I);
        if nI < 4
            I = [I, ones(ns,4-nI)];
        end

        % Sort main effects and interactions
        fmain = struct('fnum',{});
        inter = struct('fnums',{});
        for k=1:numel(job.des.fblock.maininters)
            if isfield(job.des.fblock.maininters{k},'fmain')
                fmain(end+1)=job.des.fblock.maininters{k}.fmain;
            elseif isfield(job.des.fblock.maininters{k},'inter')
                inter(end+1)=job.des.fblock.maininters{k}.inter;
            end;
        end;

        % Create main effects
        H=[];Hnames=[];
        nmain=length(fmain);
        for f=1:nmain,
            fcol=fmain(f).fnum;
            fname=job.des.fblock.fac(fcol).name;

            % Augment H partition - explicit factor numbers are 1 lower than in I matrix
            [Hf,Hfnames]=spm_DesMtx(I(:,fcol+1),'-',fname);
            H=[H,Hf];
            Hnames=[Hnames;Hfnames];
        end

        % Create interactions
        ni=length(inter);
        for i=1:ni,
            % Get the two factors for this interaction
            fnums=inter(i).fnums;
            f1=fnums(1);f2=fnums(2);

            % Names
            iname{1}=job.des.fblock.fac(f1).name;
            iname{2}=job.des.fblock.fac(f2).name;

            % Augment H partition - explicit factor numbers are 1 lower than in I matrix
            Isub=[I(:,f1+1),I(:,f2+1)];
            [Hf,Hfnames]=spm_DesMtx(Isub,'-',iname);
            H=[H,Hf];
            Hnames=[Hnames;Hfnames];

        end

        if nmain==0 && ni==0
            disp('Error in design specification: You have not specified any main effects or interactions');
            return
        end

        for i=1:nf,
            % Names and levels
            SPM.factor(i).name=job.des.fblock.fac(i).name;
            SPM.factor(i).levels=length(unique(I(:,i+1)));

            % Ancova options
            SPM.factor(i).gmsca=job.des.fblock.fac(i).gmsca;
            SPM.factor(i).ancova=job.des.fblock.fac(i).ancova;

            % Nonsphericity options
            SPM.factor(i).variance=job.des.fblock.fac(i).variance;
            SPM.factor(i).dept=job.des.fblock.fac(i).dept;

        end


end
nScan=size(I,1); %-#obs

% Set up data structures for non-sphericity routine
SPM.xVi.I=I;
SPM = spm_get_vc(SPM);

%-Covariate partition(s): interest (C) & nuisance (G) excluding global
%===================================================================
dstr   = {'covariate','nuisance variable'};
C  = []; Cnames = [];  %-Covariate DesMtx partitions & names
G  = []; Gnames = [];

xC = [];                         %-Struct array to hold raw covariates

% Covariate options:
nc=length(job.cov); % number of covariates
for i=1:nc,

    c      = job.cov(i).c;
    cname  = job.cov(i).cname;
    rc     = c;                         %-Save covariate value
    rcname = cname;                     %-Save covariate name
    if job.cov(i).iCFI==1,
        iCFI=1;
    else
        % SPMs internal factor numbers are 1 higher than specified in user
        % interface as, internally, the first factor is always `replication'
        iCFI=job.cov(i).iCFI+1;
    end
    switch job.cov(i).iCC,
        case 1
            iCC=1;
        case {2,3,4}
            iCC=job.cov(i).iCC+1;
        otherwise
            iCC=job.cov(i).iCC+3;
    end

    %-Centre within factor levels as appropriate
    if any(iCC == [1:7]),
        c = c - spm_meanby(c,eval(CCforms{iCC}));
    end

    %-Do any interaction (only for single covariate vectors)
    %-----------------------------------------------------------
    if iCFI > 1             %-(NB:iCFI=1 if size(c,2)>1)
        tI        = [eval(CFIforms{iCFI,1}),c];
        tConst    = CFIforms{iCFI,2};
        tFnames   = [eval(CFIforms{iCFI,3}),{cname}];
        [c,cname] = spm_DesMtx(tI,tConst,tFnames);
    elseif size(c,2)>1          %-Design matrix block
        [null,cname] = spm_DesMtx(c,'X',cname);
    else
        cname = {cname};
    end

    %-Store raw covariate details in xC struct for reference
    %-Pack c into appropriate DesMtx partition
    %-----------------------------------------------------------
    %-Construct description string for covariate
    str = {sprintf('%s',rcname)};
    if size(rc,2)>1, str = {sprintf('%s (block of %d covariates)',...
            str{:},size(rc,2))}; end
    if iCC < 8, str=[str;{['used centered ',sCC{iCC}]}]; end
    if iCFI> 1, str=[str;{['fitted as interaction ',sCFI{iCFI}]}]; end

    typ = 1;
    tmp       = struct( 'rc',rc,    'rcname',rcname,...
        'c',c,      'cname',{cname},...
        'iCC',iCC,  'iCFI',iCFI,...
        'type',typ,...
        'cols',[1:size(c,2)] + ...
        size([H,C],2) +  ...
        size([B,G],2)*min(typ-1,1),...
        'descrip',{str}             );
    if isempty(xC), xC = tmp; else, xC = [xC,tmp]; end
    C     = [C,c];
    Cnames = [Cnames; cname];

end
clear c tI tConst tFnames

xGX=[];
xM=[];



%===================================================================
% - C O N F I G U R E   D E S I G N -
%===================================================================

%-Images & image info: Map Y image files and check consistency of
% dimensions and orientation / voxel size
%===================================================================
fprintf('%-40s: ','Mapping files')                               %-#
VY    = spm_vol(char(P));

%-Check compatability of images (Bombs for single image)
%-------------------------------------------------------------------
spm_check_orientations(VY);

fprintf('%30s\n','...done')                                      %-#

%-Global values, scaling and global normalisation
%===================================================================
%-Compute global values
%-------------------------------------------------------------------
switch strvcat(fieldnames(job.globalc))
    case 'g_omit',
        iGXcalc=1;
    case 'g_user',
        iGXcalc=2;
    case 'g_mean',
        iGXcalc=3;
end

switch job.globalm.glonorm
    case 1,
        iGloNorm=9;
    case 2,
        iGloNorm=8;
    case 3,
        iGloNorm=1;
end
if SPM.factor(1).levels > 1
    % Over-ride if factor-specific ANCOVA has been specified
    for i=1:length(SPM.factor),
        if SPM.factor(i).ancova
            iGloNorm=i+2;
        end
    end
end

%-Analysis threshold mask
%-------------------------------------------------------------------
%-Work out available options:
% -Inf=>None, real=>absolute, complex=>proportional, (i.e. times global)
M_T = -Inf;
switch strvcat(fieldnames(job.masking.tm)),
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

if (any(iGloNorm == [1:5]) || iGloNorm==8) && iGXcalc==1
    % Over-ride omission of global calculation if we need it
    disp(' ');
    disp(sprintf('For %s, SPM needs estimates of global activity.',sGloNorm{iGloNorm}));
    disp('But you have specified to omit this computation.');
    disp('SPM has overridden this omission and will automatically compute ');
    disp('globals as the mean value of within brain voxels');
    disp(' ');
    iGXcalc=3;
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
        g = zeros(nScan,1 );
        fprintf('%-40s: %30s','Calculating globals',' ')             %-#
        for i = 1:nScan
            str = sprintf('%3d/%-3d',i,nScan);
            fprintf('%s%30s',repmat(sprintf('\b'),1,30),str)%-#
            g(i) = spm_global(VY(i));
        end
        fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')       %-#
    otherwise
        error('illegal iGXcalc')
end
rg = g;

fprintf('%-40s: ','Design configuration')                        %-#

%-Grand mean scaling options                                 (GMsca)
%-------------------------------------------------------------------
if iGloNorm==8
    iGMsca=8;   %-grand mean scaling implicit in PropSca GloNorm
else
    switch strvcat(fieldnames(job.globalm.gmsca))
        case 'gmsca_yes',
            iGMsca=1;
        case 'gmsca_no',
            iGMsca=9;
    end
    if SPM.factor(1).levels > 1
        % Over-ride if factor-specific scaling has been specified
        for i=1:length(SPM.factor),
            if SPM.factor(i).gmsca
                iGMsca=i+2;
            end
        end
    end
end

%-Value for PropSca / GMsca                                     (GM)
%-------------------------------------------------------------------
switch iGMsca,
    case 9                      %-Not scaling (GMsca or PropSca)
        GM = 0;                         %-Set GM to zero when not scaling
    case 1                                %-Ask user value of GM
        GM = job.globalm.gmsca.gmsca_yes.gmscv;
    otherwise
        if iGloNorm==8
            switch strvcat(fieldnames(job.globalm.gmsca))
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
%-------------------------------------------------------------------
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
%-------------------------------------------------------------------
switch iGMsca,
    case 8
        %-Proportional scaling global normalisation
        if iGloNorm~=8, error('iGloNorm-iGMsca(8) mismatch for PropSca'), end
        gSF    = GM./g;
        g      = GM*ones(nScan,1);
    case {1,2,3,4,5,6,7}
        %-Grand mean scaling according to iGMsca
        gSF    = GM./spm_meanby(g,eval(CCforms{iGMsca}));
        g      = g.*gSF;
    case 9
        %-No grand mean scaling
        gSF    = ones(nScan,1);
    otherwise
        error('illegal iGMsca')
end


%-Apply gSF to memory-mapped scalefactors to implement scaling
%-------------------------------------------------------------------
for i = 1:nScan
    VY(i).pinfo(1:2,:) = VY(i).pinfo(1:2,:)*gSF(i);
end

%-Global centering (for AnCova GloNorm)                         (GC)
%-If not doing AnCova then GC is irrelevant
if ~any(iGloNorm == [1:7])
    iGC = 12;
    gc  = [];
else
    iGC = 10;
    gc = 0;
end

%-AnCova: Construct global nuisance covariates partition (if AnCova)
%-------------------------------------------------------------------
if any(iGloNorm == [1:7])

    %-Centre global covariate as requested
    %---------------------------------------------------------------
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
    %---------------------------------------------------------------
    rcname     = 'global';
    tI         = [eval(CFIforms{iGloNorm,1}),g - gc];
    tConst     = CFIforms{iGloNorm,2};
    tFnames    = [eval(CFIforms{iGloNorm,3}),{rcname}];
    [f,gnames]  = spm_DesMtx(tI,tConst,tFnames);
    clear tI tConst tFnames

    %-Save GX info in xC struct for reference
    %---------------------------------------------------------------
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
    if isempty(xC), xC = tmp; else, xC = [xC,tmp]; end


elseif iGloNorm==8 || iGXcalc>1

    %-Globals calculated, but not AnCova: Make a note of globals
    %---------------------------------------------------------------
    if iGloNorm==8
        str = { 'global values: (used for proportional scaling)';...
            '("raw" unscaled globals shown)'};
    elseif isfinite(M_T) && ~isreal(M_T)
        str = { 'global values: (used to compute analysis threshold)'};
    else
        str = { 'global values: (computed but not used)'};
    end

    rcname ='global';
    tmp     = struct(   'rc',rg,    'rcname',rcname,...
        'c',{[]},   'cname' ,{{}},...
        'iCC',0,    'iCFI'  ,0,...
        'type',     3,...
        'cols',     {[]},...
        'descrip',  {str}           );

    if isempty(xC), xC = tmp; else, xC = [xC,tmp]; end
end


%-Save info on global calculation in xGX structure
%-------------------------------------------------------------------
xGX = struct(...
    'iGXcalc',iGXcalc,  'sGXcalc',sGXcalc,  'rg',rg,...
    'iGMsca',iGMsca,    'sGMsca',sGMsca,    'GM',GM,'gSF',gSF,...
    'iGC',  iGC,        'sGC',  sCC{iGC},   'gc',   gc,...
    'iGloNorm',iGloNorm,    'sGloNorm',sGloNorm);

%-Make a description string
%-------------------------------------------------------------------
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
%-------------------------------------------------------------------
if isreal(M_T),
    M_TH = M_T  * ones(nScan,1);    %-NB: -Inf is real
else
    M_TH = imag(M_T) * (rg.*gSF);
end

%-Implicit masking: Ignore zero voxels in low data-types?
%-------------------------------------------------------------------
% (Implicit mask is NaN in higher data-types.)
type = getfield(spm_vol(P{1,1}),'dt')*[1,0]';
if ~spm_type(type,'nanrep')
    M_I = job.masking.im;  % Implicit mask ?
    if M_I,
        xsM.Implicit_masking = 'Yes: zero''s treated as missing';
    else,
        xsM.Implicit_masking = 'No';
    end
else
    M_I = 1;
    xsM.Implicit_masking = 'Yes: NaN''s treated as missing';
end

%-Explicit masking
%-------------------------------------------------------------------
if isempty(job.masking.em{:})
    VM = [];
    xsM.Explicit_masking = 'No';
else
    VM = spm_vol(char(job.masking.em));
    xsM.Explicit_masking = 'Yes';
end

xM     = struct('T',M_T, 'TH',M_TH, 'I',M_I, 'VM',{VM}, 'xs',xsM);

%-Construct full design matrix (X), parameter names and structure (xX)
%===================================================================
X      = [H C B G];
tmp    = cumsum([size(H,2), size(C,2), size(B,2), size(G,2)]);
xX     = struct(    'X',        X,...
    'iH',       [1:size(H,2)],...
    'iC',       [1:size(C,2)] + tmp(1),...
    'iB',       [1:size(B,2)] + tmp(2),...
    'iG',       [1:size(G,2)] + tmp(3),...
    'name',     {[Hnames; Cnames; Bnames; Gnames]},...
    'I',        I,...
    'sF',       {sF});


%-Design description (an nx2 cellstr) - for saving and display
%===================================================================
tmp = { sprintf('%d condition, +%d covariate, +%d block, +%d nuisance',...
    size(H,2),size(C,2),size(B,2),size(G,2));...
    sprintf('%d total, having %d degrees of freedom',...
    size(X,2),rank(X));...
    sprintf('leaving %d degrees of freedom from %d images',...
    size(X,1)-rank(X),size(X,1))                };
xsDes = struct( 'Design',           {DesName},...
    'Global_calculation',       {sGXcalc},...
    'Grand_mean_scaling',       {sGMsca},...
    'Global_normalisation',     {sGloNorm},...
    'Parameters',           {tmp}           );


fprintf('%30s\n','...done')                                      %-#

%-Assemble SPM structure
%===================================================================
SPM.xY.P    = P;            % filenames
SPM.xY.VY   = VY;           % mapped data
SPM.nscan   = size(xX.X,1); % scan number
SPM.xX      = xX;           % design structure
SPM.xC      = xC;           % covariate structure
SPM.xGX     = xGX;          % global structure
SPM.xM      = xM;           % mask structure
SPM.xsDes   = xsDes;        % description

% Automatic contrast generation only works for 'Full factorials'
if ~strcmp(DesName,'Full factorial')
    % Remove the .factor field to prevent attempted automatic contrast generation
    SPM=rmfield(SPM,'factor');
end

%-Save SPM.mat and set output argument
%-------------------------------------------------------------------
fprintf('%-40s: ','Saving SPM configuration')                    %-#

if spm_matlab_version_chk('7') >= 0
    save('SPM', 'SPM', '-V6');
else
    save('SPM', 'SPM');
end;
fprintf('%30s\n','...SPM.mat saved')                             %-#
varargout = {SPM};

%-Display Design report
%===================================================================
fprintf('%-40s: ','Design reporting')                            %-#
fname     = cat(1,{SPM.xY.VY.fname}');
spm_DesRep('DesMtx',SPM.xX,fname,SPM.xsDes)
fprintf('%30s\n','...done')

out.spmmat{1} = fullfile(pwd, 'SPM.mat');
cd(original_dir); % Change back dir
fprintf('Done\n')
