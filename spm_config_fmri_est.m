function conf = spm_config_fmri_est
% Configuration file for estimation of fMRI model
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Darren Gitelman and Will Penny
% $Id: spm_config_fmri_est.m 832 2007-06-22 11:33:31Z will $


% Define inline types.
%-----------------------------------------------------------------------

entry = inline(['struct(''type'',''entry'',''name'',name,'...
    '''tag'',tag,''strtype'',strtype,''num'',num,''help'',hlp)'],...
    'name','tag','strtype','num','hlp');

files = inline(['struct(''type'',''files'',''name'',name,'...
    '''tag'',tag,''filter'',fltr,''num'',num,''help'',hlp)'],...
    'name','tag','fltr','num','hlp');

mnu = inline(['struct(''type'',''menu'',''name'',name,'...
    '''tag'',tag,''labels'',{labels},''values'',{values},''help'',hlp)'],...
    'name','tag','labels','values','hlp');

branch = inline(['struct(''type'',''branch'',''name'',name,'...
    '''tag'',tag,''val'',{val},''help'',hlp)'],...
    'name','tag','val','hlp');

repeat = inline(['struct(''type'',''repeat'',''name'',name,'...
    '''tag'',tag,''values'',{values},''help'',hlp)'],...
    'name','tag','values','hlp');

choice = inline(['struct(''type'',''choice'',''name'',name,'...
    '''tag'',tag,''values'',{values},''help'',hlp)'],...
    'name','tag','values','hlp');

%-----------------------------------------------------------------------

sp_text = ['                                                      ',...
           '                                                      '];

%-------------------------------------------------------------------------

spm.type = 'files';
spm.name = 'Select SPM.mat';
spm.tag  = 'spmmat';
spm.num  = [1 1];
spm.filter  = 'mat';
spm.ufilter = '^SPM\.mat$';
spm.help   = {'Select the SPM.mat file that contains the design specification. ',...
        'The directory containing this file is known as the input directory.'};

 
% Bayesian estimation over slices or whole volume ?

slices  = entry('Slices','Slices','e',[Inf 1],'');
p1=['Enter Slice Numbers. This can be a single slice or multiple slices. ',...
    'If you select a single slice or only a few slices you must be aware of the ',...
    'interpolation options when, after estimation, displaying the estimated ',...
    'images eg. images of contrasts or AR maps. ',...
    'The default interpolation option may need to be changed to nearest neighbour (NN) ',...
    '(see bottom right hand of graphics window) for you slice maps to be visible.'];
slices.help={p1};

volume  = struct('type','const','name','Volume','tag','Volume','val',{{1}});
p1=['You have selected the Volume option. SPM will analyse fMRI ',...
    'time series in all slices of each volume.'];
volume.help={p1};

space   = choice('Analysis Space','space',{volume,slices},'');
p1=['Because estimation can be time consuming an option is provided to analyse ',...
   'selected slices rather than the whole volume.'];
space.help={p1};
space.val={volume};

% Regression coefficient  priors for Bayesian estimation

w_prior  = mnu('Signal priors','signal',{'GMRF','LORETA','Global','Uninformative'},...
    {'GMRF','LORETA','Global','Uninformative'},{'Signal priors'});
w_prior.val={'GMRF'};
p1=['[GMRF] Gaussian Markov Random Field. This spatial prior is the recommended option. ',...
        'Regression coefficients at a given voxel are (softly) constrained to be similar to those at ',...
        'nearby voxels. The strength of this constraint is determined by a spatial precision ',...
        'parameter that is estimated from the data. Different regression coefficients have ',...
        'different spatial precisions allowing each putative experimental effect to have its ',...
        'own spatial regularity. '];
p2=['[LORETA] Low resolution Tomography Prior. This spatial prior is very similar to the GMRF ',...
        'prior and is a standatd choice for EEG source localisation algorithms. It does, however, ',...
        'have undesirable edge effects.'];
p3=['[Global] Global Shrinkage prior. This is not a spatial prior in the sense that ',...
        'regression coefficients are constrained to be similar to neighboring voxels. ',...
        'Instead, the average effect over all voxels (global effect) is assumed to be zero ',...
        'and all regression coefficients are shrunk towards this value in proporation to the ',...
        'prior precision. This is the same prior that is used for Bayesian estimation at the ',...
        'second level models, except that here the prior precision is estimated separaetly ',...
        'for each slice. '];
p4=['[Uninformative] A flat prior. Essentially, no prior information is used. ',...
        'If you select this option then VB reduces to Maximum Likelihood (ML)',...
        'estimation. This option is useful if, for example, you do not wish to ',...
        'use a spatial prior but wish to take advantage of the voxel-wise AR(P) ',...
        'modelling of noise processes. In this case, you would apply the algorithm ',...
        'to images that have been spatially smoothed. For P=0, ML estimation ',...
        'in turn reduces to Ordinary Least Squares (OLS) estimates, and for P>0 ',...
        'ML estimation is equivalent to a weighted least squares (WLS) but where the ',...
        'weights are different at each voxel (reflecting the different noise correlation ',...
        'at each voxel). '];
w_prior.help={p1,sp_text,p2,sp_text,p3,sp_text,p4};

% AR model order for Bayesian estimation

arp   = entry('AR model order','ARP','e',[Inf 1],'Enter AR model order');
arp.val={3};
p1=['An AR model order of 3 is the default. Cardiac and respiratory artifacts ',...
        'are periodic in nature and therefore require an AR order of at least 2. In ',...
        'previous work, voxel-wise selection of the optimal model order showed that a ',...
        'value of 3 was the highest order required. '];
p2=['Higher model orders have little ',...
        'effect on the estimation time. If you select a model order of zero this ',...
        'corresponds to the assumption that the errors are IID. This AR specification ',...
        'overrides any choices that were made in the model specification stage.'];
p3=['Voxel-wise AR models are fitted separately for each session of data. For each ',...
        'session this therefore produces maps of AR(1), AR(2) etc coefficients in the ',...
        'output directory. '];
arp.help={p1,sp_text,p2,sp_text,p3};

% AR coefficient  priors for Bayesian estimation

a_gmrf = struct('type','const','name','GMRF','tag','GMRF','val',{{1}});
p1=['[GMRF] Gaussian Markov Random Field. This is the default option. ',...
        'This spatial prior is the same as that used for the ',...
        'regression coefficients. Spatial precisions are estimated separately for each ',...
        'AR coefficient eg. the AR(1) coefficient over space, AR(2) over space etc. '];
a_gmrf.help={p1};

a_loreta = struct('type','const','name','LORETA','tag','LORETA','val',{{1}});
p1=['[LORETA] Low resolution Tomography Prior. See comments on LORETA priors ',...
        'for regresion coefficients.'];
a_loreta.help={p1};

a_tissue_type = files('Tissue-type','tissue_type','image',[1 Inf],'Select tissue-type images');
p1=['[Tissue-type] AR estimates at each voxel are biased towards typical ',...
    'values for that tissue type (eg. gray, white, CSF). If you select this ',...
    'option you will need to then ',...
    'select files that contain tissue type maps (see below). ',...
    'These are typically chosen to be ',...
    'Grey Matter, White Matter and CSF images derived from segmentation of ',...
    'registered structural scans.'];
p2=['Previous work has shown that ',...
    'there is significant variation in AR values with tissue type. However, GMRF priors ',...
    'have previously been favoured by Bayesian model comparison.'];
a_tissue_type.help={p1,sp_text,p2};

a_robust = struct('type','const','name','Robust','tag','Robust','val',{{1}});
p1=['Robust GLM. Uses Mixture of Gaussians noise model.'];
a_robust.help={p1};

a_prior = choice('Noise priors','noise',{a_gmrf,a_loreta,a_tissue_type,a_robust},'Noise priors');
a_prior.val={a_gmrf};
a_prior.help={'There are four noise prior options here (1) GMRF, (2) LORETA ',...
        '(3) Tissue-type and (4) Robust'};

% ANOVA options

first  = mnu('First level','first',...
    {'No','Yes'},{'No','Yes'},{''});
first.val={'No'};
p1=['This is implemented using Bayesian model comparison. For example, to test for ',...
        'the main effect of a factor two models are compared, one where the levels are ',...
        'represented using different regressors and one using the same regressor. ',...
    'This therefore requires explicit fitting of several models at each voxel and is ',...
    'computationally demanding (requiring several hours of computation). ',...
    'The recommended option is therefore NO.'];
p2=['To use this option you must have already specified your factorial design ',...
        'during the model specification stage. '];
first.help={p1,sp_text,p2};

second  = mnu('Second level','second',...
    {'No','Yes'},{'No','Yes'},{''});
second.val={'Yes'};
p1=['This option tells SPM to automatically generate ',...
    'the simple contrasts that are necessary to produce the contrast images ',...
    'for a second-level (between-subject) ANOVA. Naturally, these contrasts can ',...
    'also be used to characterise simple effects for each subject. '];
p2= ['With the Bayesian estimation ',...
    'option it is recommended that contrasts are computed during the parameter ',...
    'estimation stage (see ''simple contrasts'' below). ',...
    'The recommended option here is therefore YES.'];
p3=['To use this option you must have already specified your factorial design ',...
        'during the model specification stage. '];
p4=['If you wish to use these contrast images for a second-level analysis then ',...
    'you will need to spatially smooth them to take into account between-subject ',...
    'differences in functional anatomy ie. the fact that one persons V5 may be in ',...
    'a different position than anothers. '];
second.help={p1,sp_text,p2,sp_text,p3,sp_text,p4};

anova.type   = 'branch';
anova.name   = 'ANOVA';
anova.tag    = 'anova';
anova.val    = {first,second};
anova.help = {'Perform 1st or 2nd level Analysis of Variance.'};

% Contrasts to be computed during Bayesian estimation

name.type    = 'entry';
name.name    = 'Name';
name.tag     = 'name';
name.strtype = 's';
name.num     = [1 1];
name.help    = {'Name of contrast eg. ''Positive Effect'''};

gconvec   = entry('Contrast vector','convec','e',[Inf 1],'');
p1=['These contrasts are used to generate PPMs which ',...
        'characterise effect ',...
        'sizes at each voxel. This is in contrast to SPMs in which eg. maps of t-statistics ',...
        'show the ratio of the effect size to effect variability (standard deviation). ',...
        'SPMs are therefore a-dimensional. This is not the case for PPMs as the size of the ',...
        'effect is of primary interest. Some care is therefore needed about the scaling of ',...
        'contrast vectors. For example, if you are interested in the differential effect ',...
        'size averaged over conditions then the contrast 0.5 0.5 -0.5 -0.5 would be more ',...
        'suitable than the 1 1 -1 -1 contrast which looks at the differential effect size ',...
        'summed over conditions. '];
gconvec.help    = {p1};
            
gcon.type   = 'branch';
gcon.name   = 'Simple contrast';
gcon.tag    = 'gcon';
gcon.val    = {name,gconvec};
gcon.help = {''};
            
contrast.type = 'repeat';
contrast.name = 'Simple contrasts';
contrast.tag  = 'contrasts';
contrast.values = {gcon};
p1=['''Simple'' contrasts refers to a contrast that spans one-dimension ie. ',...
    'to assess an effect that is increasing or decreasing.'];
p2 =['If you have a factoral design then the contrasts needed to generate ',...
     'the contrast images for a 2nd-level ANOVA (or to assess these simple effects ',...
     'within-subject) can be specified automatically ',...
     'using the ANOVA->Second level option.'];
p3 =['When using the Bayesian estimation option it is computationally more ',...
     'efficient to compute the contrasts when the parameters are estimated. ',...
     'This is because estimated parameter vectors have potentially different ',...
     'posterior covariance matrices at different voxels ',...
     'and these matrices are not stored. If you compute contrasts ',...
     'post-hoc these matrices must be recomputed (an approximate reconstruction ',...
     'based on a Taylor series expansion is used). ',...
     'It is therefore recommended to specify as many contrasts as possible ',...
     'prior to parameter estimation.'];
p4=['If you wish to use these contrast images for a second-level analysis then ',...
    'you will need to spatially smooth them to take into account between-subject ',...
    'differences in functional anatomy ie. the fact that one persons V5 may be in ',...
    'a different position than anothers. '];
contrast.help={p1,sp_text,p2,sp_text,p3,sp_text,p4};

% Bayesian estimation

est_bayes2  = struct('type','const','name','Bayesian 2nd-level','tag','Bayesian2','val',{{1}});
p1=['Bayesian estimation of 2nd level models. This option uses the Empirical Bayes ',...
        'algorithm with global shrinkage priors that was previously implemented in SPM2. ',...
        'Use of the global shrinkage prior embodies a prior belief that, on average over all ',...
        'voxels, there is no net experimental effect. Some voxels will respond negatively ',...
        'and some positively with a variability determined by the prior precision. ',...
        'This prior precision can be estimated from the data using Empirical Bayes. '];
est_bayes2.help={p1};

est_bayes1 = branch('Bayesian 1st-level','Bayesian',{space,w_prior,arp,a_prior,anova,contrast},'Bayesian Estimation');
bayes_1 = ['Model parameters are estimated using Variational Bayes (VB). ',...
     'This allows you to specify spatial priors for regression coefficients ',...
     'and regularised voxel-wise AR(P) models for fMRI noise processes. ',...
     'The algorithm does not require functional images to be spatially smoothed. ',...
     'Estimation will take about 5 times longer than with the classical approach. ',...
     'This is why VB is not the default estimation option. '];
p2=['Model estimation using this option is only efficient ',...
     'if MATLAB can load a whole slice of data into physical memory. With modern PCs ',...
     'this is usually achieved if the within-plane voxel sizes are 3 by 3 mm. This is therefore ',...
     'the minimum recommended voxel size that your spatial normalisation process should ',...
     'produce. Within-plane voxel sizes of 2 by 2 mm usually result in too many voxels per slice and result ',...
     'in estimation times lasting several hours or days. Such a high resolution is therefore to be ',...
     'avoided. '];
bayes_2 = ['After estimation, contrasts are used to find regions with effects larger ',...
      'than a user-specified size eg. 1 per cent of the global mean signal. ',...
      'These effects are assessed statistically using a Posterior Probability Map (PPM).'];
est_bayes1.help={bayes_1,sp_text,p2,sp_text,bayes_2};

% Classical (ReML) estimation

est_class = struct('type','const','name','Classical','tag','Classical','val',{{1}});
classical_1 =['Model parameters are estimated using Restricted Maximum ',...
     'Likelihood (ReML). This assumes the error correlation structure is the ',...
     'same at each voxel. This correlation can be specified using either an ',...
     'AR(1) or an Independent and Identically Distributed (IID) error ',...
     'model. These options are chosen at the model specification stage. ',...
     'ReML estimation should be applied to spatially ',...
     'smoothed functional images.'];
classical_2 = ['After estimation, specific profiles of parameters are tested using a linear ',...
      'compound or contrast with the T or F statistic. The resulting statistical map ',...
      'constitutes an SPM. The SPM{T}/{F} is then characterised in terms of ',...
      'focal or regional differences by assuming that (under the null hypothesis) ',...
      'the components of the SPM (ie. residual fields) behave as smooth stationary ',...
      'Gaussian fields.'];
est_class.help={classical_1,sp_text,classical_2};

% Select method of estimation - Bayesian or classical

meth   = choice('Method','method',{est_class,est_bayes1,est_bayes2},{'Type of estimation procedure'});
meth.val={est_class};
p1=['There are three possible estimation procedures for fMRI models (1) classical (ReML) estimation of ',...
     'first or second level models, (2) Bayesian estimation of first level models and ',...
     '(3) Bayesian estimation of second level models. ',...
     'Option (2) uses a Variational Bayes (VB) algorithm that is new to SPM5. ',...
     'Option (3) uses the Empirical ',...
     'Bayes algorithm with global shrinkage priors that was also in SPM2. '];
p2=['To use option (3) you must have already estimated the model using option (1). ',...
        'That is, for second-level models you must run a ReML estimation before ',...
        'running a Bayesian estimation. This is not necessary for option (2). Bayesian ',...
        'estimation of 1st-level models using VB does not require a prior ReML estimation.'];
meth.help={p1,sp_text,p2};


%-------------------------------------------------------------------------

conf = branch('Model estimation','fmri_est',...
    {spm,meth},'Model estimation');
conf.prog   = @run_est;
conf.vfiles = @vfiles_stats;
conf.modality = {'FMRI','PET'};
p1 = [...
  'Model parameters can be estimated using ',...
  'classical (ReML - Restricted Maximum Likelihood) or Bayesian algorithms. ',...
  'After parameter estimation, the RESULTS button can be used to specify ',...
  'contrasts that will produce Statistical Parametric Maps (SPMs) or ',...
  'Posterior Probability Maps (PPMs) and tables of statistics.'];

conf.help = {p1};

return;
%=======================================================================

%=======================================================================
function run_est(job)
% Set up the design matrix and run a design.

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
    SPM = spm_spm_Bayes(SPM);
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
    
    cd(original_dir); % Change back
    fprintf('Done\n');
    return
end

%=======================================================================
% B A Y E S I A N   1st   L E V E L   E S T I M A T I O N
%=======================================================================

%-Analyse specific slices or whole volume
%-----------------------------------------------------------------------
if isfield(job.method.Bayesian.space,'Slices')
    SPM.PPM.space_type = 'Slices';
    SPM.PPM.AN_slices  = job.method.Bayesian.space.Slices;
else
    SPM.PPM.space_type = 'Volume';
end

%-Regression coefficient priors
%-----------------------------------------------------------------------
switch job.method.Bayesian.signal
    case 'GMRF',
        SPM.PPM.priors.W = 'Spatial - GMRF';
    case 'LORETA',
        SPM.PPM.priors.W = 'Spatial - LORETA';
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
if isfield(job.method.Bayesian.noise,'GMRF')
    SPM.PPM.priors.A  = 'Spatial - GMRF';
elseif isfield(job.method.Bayesian.noise,'LORETA')
    SPM.PPM.priors.A  = 'Spatial - LORETA';
elseif isfield(job.method.Bayesian.noise,'tissue_type')
    SPM.PPM.priors.A  = 'Discrete';
    SPM.PPM.priors.SY = job.method.Bayesian.noise.tissue_type;
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

cd(original_dir); % Change back
fprintf('Done\n')
return


%=======================================================================
function vf = vfiles_stats(job)
vf = {job.spmmat{:}};

% Should really create a few vfiles for beta images etc here as well.
