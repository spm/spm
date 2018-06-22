function out = spm_run_fmri_est(job)
% Estimate parameters of a specified model
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_fmri_est.m 7354 2018-06-22 10:44:22Z guillaume $


%-Load SPM.mat file
%--------------------------------------------------------------------------
load(job.spmmat{1},'SPM');
if ~exist('SPM','var')
    error('The MAT-file does not contain an SPM variable.');
end
out.spmmat = job.spmmat;

%-Move to the directory where the SPM.mat file is
%--------------------------------------------------------------------------
original_dir = pwd;
cd(fileparts(job.spmmat{:}));

%==========================================================================
%                       R E M L   E S T I M A T I O N
%==========================================================================
if isfield(job.method,'Classical')
    
    %-ReML estimation of the model
    %----------------------------------------------------------------------
    SPM = spm_spm(SPM);
    
    %-Automatically set up contrasts for factorial designs
    %----------------------------------------------------------------------
    if isfield(SPM,'factor') && ~isempty(SPM.factor) && SPM.factor(1).levels > 1

        Ic = [];
        
        %-Generate contrasts
        %------------------------------------------------------------------
        cons = spm_design_contrasts(SPM);

        %-Create F-contrasts
        %------------------------------------------------------------------
        for i=1:length(cons)
            con       = cons(i).c;
            name      = cons(i).name;
            STAT      = 'F';
            [c,I]     = spm_conman('ParseCon',con,SPM.xX.xKXs,STAT);
            if all(I)
                DxCon = spm_FcUtil('Set',name,STAT,'c',c,SPM.xX.xKXs);
                if isempty(SPM.xCon)
                    SPM.xCon = DxCon;
                else
                    SPM.xCon(end+1) = DxCon;
                end
                Ic = [Ic length(SPM.xCon)];
            end
        end

        %-Create t-contrasts
        %------------------------------------------------------------------
        for i=1:length(cons)
            % Create a t-contrast for each row of each F-contrast
            % The resulting contrast image can be used in a 2nd-level analysis
            Fcon          = cons(i).c;
            STAT          = 'T';
            for r=1:size(Fcon,1)
                con       = Fcon(r,:);
                str       = cons(i).name;
                if strncmp('Interaction',str,11)
                    name  = ['Positive ',str,'_',int2str(r)];
                else
                    sp1   = find(isspace(str), 1);
                    name  = ['Positive',str(sp1:end),'_',int2str(r)];
                end
                
                [c,I]     = spm_conman('ParseCon',con,SPM.xX.xKXs,STAT);
                if all(I)
                    DxCon = spm_FcUtil('Set',name,STAT,'c',c,SPM.xX.xKXs);
                    if isempty(SPM.xCon)
                        SPM.xCon = DxCon;
                    else
                        SPM.xCon(end+1) = DxCon;
                    end
                    Ic = [Ic length(SPM.xCon)];
                end
            end
        end
        
        %-Estimate constrasts
        %------------------------------------------------------------------
        if ~isempty(Ic)
            spm('FnBanner','spm_contrasts.m');
            SPM = spm_contrasts(SPM,Ic);
        end
        
    end
    
    %-Residuals
    %----------------------------------------------------------------------
    if job.write_residuals
        VRes = spm_write_residuals(SPM,NaN);
    end
    
    %-Computation results
    %----------------------------------------------------------------------
    %out.spmvar = SPM;
    out.beta  = spm_file({SPM.Vbeta(:).fname}','path',SPM.swd);
    out.mask  = {fullfile(SPM.swd,SPM.VM.fname)};
    out.resms = {fullfile(SPM.swd,SPM.VResMS.fname)};
    if job.write_residuals
        out.res = spm_file({VRes(:).fname}','path',SPM.swd);
    end
    cd(original_dir);
    fprintf('Done\n');
    return
end

if job.write_residuals
    warning('Save residuals option is only implemented for classical inference.');
    job.write_residuals = false;
end

    
%==========================================================================
%        B A Y E S I A N   2nd   L E V E L   E S T I M A T I O N
%==========================================================================
if isfield(job.method,'Bayesian2')
    
    SPM = spm_spm_Bayes(SPM);
    
    %out.spmvar = SPM;
    cd(original_dir);
    fprintf('Done\n');
    return
end


%==========================================================================
%        B A Y E S I A N   1st   L E V E L   E S T I M A T I O N
%==========================================================================

%-Analyse specific slices or whole volume
%--------------------------------------------------------------------------
switch char(fieldnames(job.method.Bayesian.space))
  case 'volume'
      SPM.PPM.space_type  = 'volume';
      SPM.PPM.block_type  = lower(job.method.Bayesian.space.volume.block_type);
  case 'slices'
      SPM.PPM.space_type  = 'slices';
      SPM.PPM.AN_slices   = job.method.Bayesian.space.slices.numbers;
      SPM.PPM.block_type  = lower(job.method.Bayesian.space.slices.block_type);
  case 'clusters'
      SPM.PPM.space_type  = 'clusters';
      SPM.PPM.clustermask = job.method.Bayesian.space.clusters.mask;
      SPM.PPM.block_type  = lower(job.method.Bayesian.space.clusters.block_type);
  otherwise
      SPM.PPM.space_type  = 'volume';
      SPM.PPM.block_type  = 'slices';
end

%-Regression coefficient priors
%--------------------------------------------------------------------------
switch job.method.Bayesian.signal
    case 'UGL'
        SPM.PPM.priors.W  = 'Spatial - UGL';
    case 'GMRF'
        SPM.PPM.priors.W  = 'Spatial - GMRF';
    case 'LORETA'
        SPM.PPM.priors.W  = 'Spatial - LORETA';
    case 'WGL'
        SPM.PPM.priors.W  = 'Spatial - WGL';
    case 'Global'
        SPM.PPM.priors.W  = 'Voxel - Shrinkage';
    case 'Uninformative'
        SPM.PPM.priors.W  = 'Voxel - Uninformative';
    otherwise
        error('Unknown prior for W in Bayesian 1st level estimation.');
end

%-Number of AR coefficients
%--------------------------------------------------------------------------
SPM.PPM.AR_P              = job.method.Bayesian.ARP;

%-AR coefficient priors
%--------------------------------------------------------------------------
if isfield(job.method.Bayesian.noise,'UGL')
    SPM.PPM.priors.A      = 'Spatial - UGL';
elseif isfield(job.method.Bayesian.noise,'GMRF')
    SPM.PPM.priors.A      = 'Spatial - GMRF';
elseif isfield(job.method.Bayesian.noise,'LORETA')
    SPM.PPM.priors.A      = 'Spatial - LORETA';
elseif isfield(job.method.Bayesian.noise,'tissue_type')
    SPM.PPM.priors.A      = 'Discrete';
    SPM.PPM.priors.SY     = char(job.method.Bayesian.noise.tissue_type);
elseif isfield(job.method.Bayesian.noise,'Robust')
    SPM.PPM.priors.A      = 'Robust';
    SPM.PPM.AR_P          = 0;
    SPM.PPM.update_F      = 1;
end

%-Define an empty contrast
%--------------------------------------------------------------------------
NullCon      = spm_FcUtil('Set','','P','c',[],1);
NullCon.X0   = [];
NullCon.iX0  = [];
NullCon.X1o  = [];
NullCon.eidf = 1;

SPM.xCon     = [];
%-Set up contrasts for 2nd-level ANOVA
%--------------------------------------------------------------------------
if strcmp(job.method.Bayesian.anova.second,'Yes') && isfield(SPM,'factor')
    cons = spm_design_contrasts(SPM);
    for i=1:length(cons)
        % Create a simple contrast for each row of each F-contrast
        % The resulting contrast image can be used in a 2nd-level analysis
        Fcon = cons(i).c;
        for r=1:size(Fcon,1)
            % Normalise contrast st. sum of positive elements is 1
            % and sum of negative elements is 1
            con = Fcon(r,:);
            con = con / length(find(con==1));

            % Change name
            str = cons(i).name;
            sp1 = find(str==' ', 1);
            if strcmp(str(1:11),'Interaction')
                name = ['Positive ',str,'_',int2str(r)];
            else
                name = ['Positive',str(sp1:end),'_',int2str(r)];
            end

            DxCon = NullCon;
            DxCon.name = name;
            DxCon.c = con';

            if isempty(SPM.xCon)
                SPM.xCon = DxCon;
            else
                SPM.xCon(end+1) = DxCon;
            end
        end
    end
end

%-Set up user-specified simple contrasts
%--------------------------------------------------------------------------
K               = size(SPM.xX.X,2);
for c = 1:length(job.method.Bayesian.gcon)
    DxCon       = NullCon;
    DxCon.name  = job.method.Bayesian.gcon(c).name;
    convec      = job.method.Bayesian.gcon(c).convec(:);
    if length(convec) == K
        DxCon.c = convec;
    elseif length(convec) < K
        fprintf('Zero padding.');                                       %-#
        DxCon.c = [convec; zeros(K-length(convec),1)];
    else
        warning(['User-specified contrast nb %d has %d entries '...
            'but there are %d regressors - ignored.'], c,length(convec),K);
        DxCon   = [];
    end
    
    if isempty(SPM.xCon)
        SPM.xCon = DxCon;
    else
        SPM.xCon(end+1) = DxCon;
    end
end

%-Compute log evidence maps
%--------------------------------------------------------------------------
if strcmp(job.method.Bayesian.LogEv,'Yes')
    SPM.PPM.update_F      = 1;
    SPM.PPM.compute_det_D = 1;
end

%-1st level Bayesian ANOVA ?
%--------------------------------------------------------------------------
bayes_anova = 0;
if strcmp(job.method.Bayesian.anova.first,'Yes')
    bayes_anova           = 1;
    SPM.PPM.update_F      = 1; % Compute evidence for each model
    SPM.PPM.compute_det_D = 1; 
end

%-Variational Bayes estimation
%--------------------------------------------------------------------------
SPM = spm_spm_vb(SPM);

%-Bayesian ANOVA using model comparison
%--------------------------------------------------------------------------
if bayes_anova
    % We don't want to estimate contrasts for each different model
    SPM.xCon = [];
    spm_vb_ppm_anova(SPM);
end

%out.spmvar = SPM;
cd(original_dir);
fprintf('Done\n')
