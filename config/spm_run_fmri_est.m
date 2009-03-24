function out = spm_run_fmri_est(job)
% Set up the design matrix and run a design.
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_fmri_est.m 2928 2009-03-24 08:54:32Z lee $


global defaults
if isempty(defaults)
    spm_defaults;
end;
if ~isfield(defaults,'modality')
    defaults.modality = 'FMRI';
end;

%-Load SPM.mat file
%-----------------------------------------------------------------------
SPM = [];
load(job.spmmat{:});
% pre-set output
out.spmmat = job.spmmat;

original_dir = pwd;

%-Move to the directory where the SPM.mat file is
%-----------------------------------------------------------------------
cd(fileparts(job.spmmat{:}));

% COMMENTED OUT BY DRG. THIS SHOULD BE TAKEN CARE OF WITHIN SPM_SPM AND
% SPM_SPM_BAYES. REMOVE THIS SECTION ONCE THIS HAS BEEN VERIFIED.
%-If we've gotten to this point we're committed to overwriting files.
% Delete them so we don't get stuck in spm_spm
%-----------------------------------------------------------------------
% files = {'^mask\..{3}$','^ResMS\..{3}$','^RPV\..{3}$',...
%          '^beta_.{4}\..{3}$','^con_.{4}\..{3}$','^ResI_.{4}\..{3}$',...
%          '^ess_.{4}\..{3}$', '^spm\w{1}_.{4}\..{3}$'};
% 
% for i=1:length(files)
%     j = spm_select('List',pwd,files{i});
%     for k=1:size(j,1)
%         spm_unlink(deblank(j(k,:)));
%     end
% end
% END COMMENTED OUT BY DRG

%=======================================================================
% B A Y E S I A N   2nd   L E V E L   E S T I M A T I O N
%=======================================================================
if isfield(job.method,'Bayesian2')
    %out.spmvar = spm_spm_Bayes(SPM);
    cd(original_dir); % Change back
    fprintf('Done\n');
    return
end

%=======================================================================
% R E M L   E S T I M A T I O N
%=======================================================================
if isfield(job.method,'Classical'),
    
    SPM = spm_spm(SPM);
    
    %-Automatically set up contrasts for factorial designs
    %-------------------------------------------------------------------
    if isfield(SPM,'factor')
        if SPM.factor(1).levels > 1
        % don't both if you've only got 1 level and 1 factor
            cons = spm_design_contrasts(SPM);
        
            %-Create F-contrasts
            %-----------------------------------------------------------
            for i=1:length(cons)
                con  = cons(i).c;
                name = cons(i).name;
                STAT = 'F';
                [c,I,emsg,imsg] = spm_conman('ParseCon',con,SPM.xX.xKXs,STAT);
                if all(I)
                    DxCon = spm_FcUtil('Set',name,STAT,'c',c,SPM.xX.xKXs);
                else
                    DxCon = [];
                end
                if isempty(SPM.xCon),
                    SPM.xCon = DxCon;
                else
                    SPM.xCon(end+1) = DxCon;
                end
               SPM = spm_contrasts(SPM,length(SPM.xCon));
            end
        
            %-Create t-contrasts
            %-----------------------------------------------------------
            for i=1:length(cons)
                % Create a t-contrast for each row of each F-contrast
                % The resulting contrast image can be used in a 2nd-level analysis
                Fcon  = cons(i).c;
                nrows = size(Fcon,1);
                STAT  = 'T';
                for r=1:nrows,
                    con = Fcon(r,:); 
                    str = cons(i).name;
                    if ~isempty(strmatch('Interaction',str))
                        name = ['Positive ',str,'_',int2str(r)];
                    else
                        sp1  = min(find(isspace(str))); 
                        name = ['Positive',str(sp1:end),'_',int2str(r)];
                    end
                
                    [c,I,emsg,imsg] = spm_conman('ParseCon',con,SPM.xX.xKXs,STAT);
                    if all(I)
                        DxCon = spm_FcUtil('Set',name,STAT,'c',c,SPM.xX.xKXs);
                    else
                        DxCon = [];
                    end
                    if isempty(SPM.xCon),
                        SPM.xCon = DxCon;
                    else
                        SPM.xCon(end+1) = DxCon;
                    end
                   SPM = spm_contrasts(SPM,length(SPM.xCon));
                end
            end
        end % if SPM.factor(1).levels > 1
    end % if isfield(SPM,'factor')
    
    %out.spmvar = SPM;
    out.beta = cellfun(@(fn)fullfile(SPM.swd,fn), cellstr(char(SPM.Vbeta(:).fname)),'UniformOutput',false);
    out.mask = {fullfile(SPM.swd,SPM.VM.fname)};
    out.resms = {fullfile(SPM.swd,SPM.VResMS.fname)};
    cd(original_dir); % Change back
    fprintf('Done\n');
    return
end

%=======================================================================
% B A Y E S I A N   1st   L E V E L   E S T I M A T I O N
%=======================================================================

%-Analyse specific slices or whole volume
%-----------------------------------------------------------------------
switch char(fieldnames(job.method.Bayesian.space))
  case 'volume'
      SPM.PPM.space_type = 'volume';
      SPM.PPM.block_type = lower(job.method.Bayesian.space.volume.block_type);
  case 'slices'
      SPM.PPM.space_type = 'slices';
      SPM.PPM.AN_slices  = job.method.Bayesian.space.slices.numbers;
      SPM.PPM.block_type = lower(job.method.Bayesian.space.slices.block_type);
  case 'clusters'
      SPM.PPM.space_type = 'clusters';
      SPM.PPM.clustermask  = job.method.Bayesian.space.clusters.mask;
      SPM.PPM.block_type = lower(job.method.Bayesian.space.clusters.block_type);
  otherwise
      SPM.PPM.space_type = 'volume';
      SPM.PPM.block_type = 'slices';
end
%-Regression coefficient priors
%-----------------------------------------------------------------------
switch job.method.Bayesian.signal
    case 'UGL',
        SPM.PPM.priors.W = 'Spatial - UGL';
    case 'GMRF',
        SPM.PPM.priors.W = 'Spatial - GMRF';
    case 'LORETA',
        SPM.PPM.priors.W = 'Spatial - LORETA';
    case 'WGL',
        SPM.PPM.priors.W = 'Spatial - WGL';
    case 'Global',
        SPM.PPM.priors.W = 'Voxel - Shrinkage';
    case 'Uninformative',
        SPM.PPM.priors.W = 'Voxel - Uninformative';
    otherwise
        error('Unkown prior for W in spm_config_fmri_est');
end

%-Number of AR coefficients
%-----------------------------------------------------------------------
SPM.PPM.AR_P = job.method.Bayesian.ARP;

%-AR coefficient priors
%-----------------------------------------------------------------------
if isfield(job.method.Bayesian.noise,'UGL')
    SPM.PPM.priors.A  = 'Spatial - UGL';
elseif isfield(job.method.Bayesian.noise,'GMRF')
    SPM.PPM.priors.A  = 'Spatial - GMRF';
elseif isfield(job.method.Bayesian.noise,'LORETA')
    SPM.PPM.priors.A  = 'Spatial - LORETA';
elseif isfield(job.method.Bayesian.noise,'tissue_type')
    SPM.PPM.priors.A  = 'Discrete';
    SPM.PPM.priors.SY = char(job.method.Bayesian.noise.tissue_type);
elseif isfield(job.method.Bayesian.noise,'Robust')
    SPM.PPM.priors.A  = 'Robust';
    SPM.PPM.AR_P=0;
    SPM.PPM.update_F=1;
end

%-Define an empty contrast
%-----------------------------------------------------------------------
NullCon      = spm_FcUtil('Set','','P','c',[],1);
NullCon.X0   = [];
NullCon.iX0  = [];
NullCon.X1o  = [];
NullCon.eidf = 1;

SPM.xCon = [];
%-Set up contrasts for 2nd-level ANOVA
%-----------------------------------------------------------------------
if strcmp(job.method.Bayesian.anova.second,'Yes')
    if isfield(SPM,'factor')
        cons=spm_design_contrasts(SPM);
        for i=1:length(cons),
            % Create a simple contrast for each row of each F-contrast
            % The resulting contrast image can be used in a 2nd-level analysis
            Fcon=cons(i).c;
            nrows=size(Fcon,1);
            STAT='P';
            for r=1:nrows,
                con=Fcon(r,:);     
                
                % Normalise contrast st. sum of positive elements is 1
                % and sum of negative elements  is 1
                s1=length(find(con==1));
                con=con./s1;
                
                % Change name 
                str=cons(i).name;
                sp1=min(find(str==' '));
                if strcmp(str(1:11),'Interaction')
                    name=['Positive ',str,'_',int2str(r)];
                else
                    name=['Positive',str(sp1:end),'_',int2str(r)];
                end
                
                DxCon=NullCon;
                DxCon.name=name;
                DxCon.c=con';
                
                if isempty(SPM.xCon),
                    SPM.xCon = DxCon;
                else
                    SPM.xCon(end+1) = DxCon;
                end
            end
        end
    end
end

%-Compute F 
%-----------------------------------------------------------------------
if strcmp(job.method.Bayesian.LogEv,'Yes')
    SPM.PPM.update_F      = 1;
    SPM.PPM.compute_det_D = 1;
end

%-Set up user-specified simple contrasts
%-----------------------------------------------------------------------
ncon = length(job.method.Bayesian.gcon);
K    = size(SPM.xX.X,2);
for c = 1:ncon
    DxCon = NullCon;
    DxCon.name = job.method.Bayesian.gcon(c).name;
    convec = job.method.Bayesian.gcon(c).convec(:);
    if length(convec) == K
        DxCon.c = convec;
    else
        str = ['Error in contrast specification:' ...
        sprintf('\n    contrast has %d entries ', length(convec)) ...
        sprintf('but there are %d regressors !\n', K)];
        error(str);
    end
    
    if isempty(SPM.xCon),
        SPM.xCon = DxCon;
    else
        SPM.xCon(end+1) = DxCon;
    end
end
    
%-1st level Bayesian ANOVA ?
%-----------------------------------------------------------------------
bayes_anova = 0;
if strcmp(job.method.Bayesian.anova.first,'Yes')
    bayes_anova           = 1;
    SPM.PPM.update_F      = 1; % Compute evidence for each model
    SPM.PPM.compute_det_D = 1; 
end

%-Variational Bayes estimation
%-----------------------------------------------------------------------
SPM = spm_spm_vb(SPM);

%-Bayesian ANOVA using model comparison
%-----------------------------------------------------------------------
if bayes_anova
    % We don't want to estimate contrasts for each different model
    SPM.xCon = [];
    spm_vb_ppm_anova(SPM);
end

%out.spmvar = SPM;
cd(original_dir); % Change back
fprintf('Done\n')
return
