function conf = spm_config_fmri_est
% Configuration file for estimation of fMRI model
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Darren Gitelman and Will Penny
% $Id$


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
spm.filter = '^SPM\.mat$';
spm.help   = {'Select SPM.mat file for contrasts'};

%-------------------------------------------------------------------------

cvi   = mnu('Serial correlations','cvi',{'none','AR(1)'},{'none','AR(1)'},...
    {'Correct for serial correlations'});
cvi.val={'AR(1)'};
p1 = [...
    'Serial correlations in fast fMRI time-series are dealt with as ',...
    'described in spm_spm.  At this stage you need to specify the filtering ',...
    'that will be applied to the data (and design matrix) to give a ',...
    'generalized least squares (GLS) estimate of the parameters required. ',...
    'This filtering is important to ensure that the GLS estimate is ',...
    'efficient and that the error variance is estimated in an unbiased way.'];
p2 = ['                                                      ',...
      '                                                      '];
p3 = [...
    'The serial correlations will be estimated with a ReML (restricted ',...
    'maximum likelihood) algorithm using an autoregressive AR(1) plus ',...
    ' white noise model during parameter estimation.  This estimate assumes the same ',...
    'correlation structure for each voxel, within each session.  The ReML ',...
    'estimates are then used to correct for non-sphericity during inference ',...
    'by adjusting the statistics and degrees of freedom appropriately.  The ',...
    'discrepancy between estimated and actual intrinsic (i.e. prior to ',...
    'filtering) correlations are greatest at low frequencies.  Therefore ',...
    'specification of the high-pass filter is particularly important.'];
cvi.help = {p1,p2,p3};
 
% Bayesian estimation over slices or whole volume ?

slices  = entry('Slices','Slices','e',[Inf 1],'Enter Slice Numbers');

volume  = struct('type','const','name','Volume','tag','Volume','val',{{1}});
p1=['You have selected the Volume option. SPM will analyse fMRI ',...
    'time series in all slices of each volume.'];
volume.help={p1};

space   = choice('Analysis Space','space',{volume,slices},'Analyse whole volume or selected slices only');
space.val={volume};

% Regression coefficient  priors for Bayesian estimation

w_prior  = mnu('Signal priors','signal',{'GMRF','LORETA','Global','Uninformative'},...
    {'GMRF','LORETA','Global','Uninformative'},{'Signal priors'});
w_prior.val={'GMRF'};
p1=['[GMRF] = Gaussian Markov Random Field. This spatial prior is the recommended option. '];
p2=['[LORETA] = Low resolution Tomography Prior. This spatial prior is popular in the EEG world. '];
p3=['[Global] = Global Shrinkage prior. This is not a spatial prior.'];
p4=['[Uninformative] = A flat prior. Essentially, no prior information is used. '];
w_prior.help={p1,p2,p3,p4};

% AR model order for Bayesian estimation

arp   = entry('AR model order','ARP','e',[Inf 1],'Enter AR model order');
arp.val={3};
p1=['An AR model order of 3 is recommended'];
arp.help={p1};

% AR coefficient  priors for Bayesian estimation

a_gmrf = struct('type','const','name','GMRF','tag','GMRF','val',{{1}});
a_loreta = struct('type','const','name','LORETA','tag','LORETA','val',{{1}});
a_tissue_type = files('Tissue-type','tissue_type','image',[1 Inf],'Select tissue-type images');
p1=['Select files that specify tissue types. These are typically chosen to be ',...
    'Grey Matter, White Matter and CSF images derived from segmentation of ',...
    'registered structural scans.'];
a_tissue_type.help={p1};
a_prior = choice('Noise priors','noise',{a_gmrf,a_loreta,a_tissue_type},'Noise priors');
a_prior.val={a_gmrf};
p1=['[GMRF] = Gaussian Markov Random Field. This spatial prior is the recommended option. '];
p2=['[LORETA] = Low resolution Tomography Prior. This spatial prior is popular in the EEG world. '];
p3=['[Tissue-type] = AR estimates at each voxel are biased towards typical ',...
    'values for that tissue type. '];
a_prior.help={p1,p2,p3};

% ANOVA options

first  = mnu('First level','first',...
    {'No','Yes'},{'No','Yes'},{''});
first.val={'No'};
p1=['[First level ANOVA ?] '];
p2=['This is implemented using Bayesian model comparison. ',...
    'This requires explicit fitting of several models at each voxel and is ',...
    'computationally demanding (requiring several hours of computation). ',...
    'The recommended option is therefore NO.'];
p3=['To use this option you must also specify your Factorial design (see options ',...
    'under FMRI Stats).'];
first.help={p1,sp_text,p2,sp_text,p3};

second  = mnu('Second level','second',...
    {'No','Yes'},{'No','Yes'},{''});
second.val={'Yes'};
p1=['[Second level ANOVA ?] '];
p2=['This option tells SPM to automatically generate ',...
    'the simple contrasts that are necessary to produce the contrast images ',...
    'for a second-level (between-subject) ANOVA. With the Bayesian estimation ',...
    'option it is recommended that contrasts are computed during the parameter ',...
    'estimation stage (see HELP for Simple contrasts). ',...
    'The recommended option here is therefore YES.'];
p3=['To use this option you must also specify your Factorial design (see options ',...
    'under FMRI Stats).'];
second.help={p1,sp_text,p2,sp_text,p3};

anova.type   = 'branch';
anova.name   = 'ANOVA';
anova.tag    = 'anova';
anova.val    = {first,second};
anova.help = {'Perform 1st or 2nd level ANOVAs'};

% Contrasts to be computed during Bayesian estimation

name.type    = 'entry';
name.name    = 'Name';
name.tag     = 'name';
name.strtype = 's';
name.num     = [1 1];
name.help    = {'Name of contrast'};

gconvec.type    = 'entry';
gconvec.name    = 'Contrast vector';
gconvec.tag     = 'convec';
gconvec.strtype = 's';
gconvec.num     = [1 1];
gconvec.help    = {''};
            
gcon.type   = 'branch';
gcon.name   = 'Simple contrast';
gcon.tag    = 'gcon';
gcon.val    = {name,gconvec};
gcon.help = {''};
            
contrast.type = 'repeat';
contrast.name = 'Simple contrasts';
contrast.tag  = 'contrasts';
contrast.values = {gcon};
p1 =['Specify simple one-dimensional contrasts'];
p2 =['If you have a factoral design then the contrasts needed to generate ',...
     'the contrast images for a 2nd-level ANOVA can be specified automatically ',...
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
contrast.help={p1,sp_text,p2,sp_text,p3};

% Bayesian estimation

est_bayes2  = struct('type','const','name','Bayesian 2nd-level','tag','Bayesian2','val',{{1}});
p1=['Bayesian analysis of 2nd level data'];
est_bayes2.help={p1};

est_bayes1 = branch('Bayesian 1st-level','Bayesian',{space,w_prior,arp,a_prior,anova,contrast},'Bayesian Estimation');
bayes_1 = ['[Bayesian 1st-level] - model parameters are estimated using Variational Bayes. ',...
     'This allows you to specify spatial priors for regression coefficients ',...
     'and regularised voxel-wise AR(P) models for fMRI noise processes. ',...
     'The algorithm does not require functional images to be spatially smoothed. ',...
     'Estimation will take about 5 times longer than with the classical approach.' ];
bayes_2 = ['After estimation, contrasts are used to find regions with effects larger ',...
      'than a user-specified size eg. 1 per cent of the global mean signal. ',...
      'These effects are assessed statistically using a Posterior Probability Map (PPM).'];
est_bayes1.help={bayes_1,sp_text,bayes_2};

% Classical (ReML) estimation

est_class = branch('Classical','Classical',{cvi},{'Classical Estimation'});
classical_1 =['[Classical] - model parameters are estimated using Restricted Maximum ',...
     'Likelihood (ReML). This assumes the error correlation structure is the ',...
     'same at each voxel. This correlation can be specified using an ',...
     'AR(1) plus white noise model. The algorithm should be applied to spatially ',...
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
meth.help={classical_1,sp_text,classical_2,sp_text,bayes_1,sp_text,bayes_2};

%-------------------------------------------------------------------------

cdir = files('Directory','dir','dir',1,'');
cdir.help = {[...
'Select an analysis directory to change to. ',...
'This is where the results of the estimation will be written.']};

%-------------------------------------------------------------------------

conf = branch('fMRI model estimation','fmri_est',...
    {spm,cdir,meth},'fMRI model estimation');
conf.prog   = @run_est;
conf.vfiles = @vfiles_stats;
conf.check  = @check_dir;
conf.modality = {'FMRI'};
p1 = [...
  'This configures the design matrix, data specification and ',...
  'filtering that specify the ensuing statistical analysis. These ',...
  'arguments are passed to spm_spm that then performs the actual parameter ',...
  'estimation.'];
p2 = [...
  'The design matrix defines the experimental design and the nature of ',...
  'hypothesis testing to be implemented.  The design matrix has one row ',...
  'for each scan and one column for each effect or explanatory variable. ',...
  '(e.g. regressor or stimulus function).  The parameters are estimated using ',...
  'Bayesian or Restricted Maximum Likelihood algorithms. Specific profiles ',...
  'within these parameters are tested using a linear compound or contrast ',...
  'with the T or F statistic. '];
p3 = [...
  'You (i) specify a statistical model in terms ',...
  'of a design matrix, (ii) associate some data with a pre-specified design ',...
  '[or (iii) specify both the data and design] and then proceed to estimate ',...
  'the parameters of the model. ',...
  'Inferences can be made about the ensuing parameter estimates (at a first ',...
  'or fixed-effect level) in the results section, or they can be re-entered ',...
  'into a second (random-effect) level analysis by treating the session or ',...
  'subject-specific [contrasts of] parameter estimates as new summary data. ',...
  'Inferences at any level obtain by specifying appropriate T or F contrasts ',...
  'in the results section to produce SPMs/PPMs and tables of statistics.'];
p4 = [...
  'You can build design matrices with separable ',...
  'session-specific partitions.  Each partition may be the same (in which ',...
  'case it is only necessary to specify it once) or different.  Responses ',...
  'can be either event- or epoch related, The only distinction is the duration ',...
  'of the underlying input or stimulus function. Mathematically they are both ',...
  'modeled by convolving a series of delta (stick) or box functions (u), ',...
  'indicating the onset of an event or epoch with a set of basis ',...
  'functions.  These basis functions model the hemodynamic convolution, ',...
  'applied by the brain, to the inputs.  This convolution can be first-order ',...
  'or a generalized convolution modeled to second order (if you specify the ',...
  'Volterra option). [The same inputs are used by the hemodynamic model or ',...
  'or dynamic causal models which model the convolution explicitly in terms of ',...
  'hidden state variables (see spm_hdm_ui and spm_dcm_ui).] ',...
  'Basis functions can be used to plot estimated responses to single events ',...
  'once the parameters (i.e. basis function coefficients) have ',...
  'been estimated.  The importance of basis functions is that they provide ',...
  'a graceful transition between simple fixed response models (like the ',...
  'box-car) and finite impulse response (FIR) models, where there is one ',...
  'basis function for each scan following an event or epoch onset.  The ',...
  'nice thing about basis functions, compared to FIR models, is that data ',...
  'sampling and stimulus presentation does not have to be synchronized ',...
  'thereby allowing a uniform and unbiased sampling of peri-stimulus time.'];
p5 = [...
  'Event-related designs may be stochastic or deterministic.  Stochastic ',...
  'designs involve one of a number of trial-types occurring with a ',...
  'specified probably at successive intervals in time.  These ',...
  'probabilities can be fixed (stationary designs) or time-dependent ',...
  '(modulated or non-stationary designs).  The most efficient designs ',...
  'obtain when the probabilities of every trial type are equal. ',...
  'A critical issue in stochastic designs is whether to include null events ',...
  'If you wish to estimate the evoke response to a specific event ',...
  'type (as opposed to differential responses) then a null event must be ',...
  'included (even if it is not modeled explicitly).'];
r1 = [...
  '  * Friston KJ, Holmes A, Poline J-B, Grasby PJ, Williams SCR, Frackowiak ',...
  'RSJ & Turner R (1995) Analysis of fMRI time-series revisited. NeuroImage ',...
  '2:45-53'];
r2 = [...
  '  * Worsley KJ and Friston KJ (1995) Analysis of fMRI time-series revisited - ',...
  'again. NeuroImage 2:178-181'];
r3 = [...
  '  * Friston KJ, Frith CD, Frackowiak RSJ, & Turner R (1995) Characterising ',...
  'dynamic brain responses with fMRI: A multivariate approach NeuroImage - ',...
  '2:166-172'];
r4 = [...
  '  * Frith CD, Turner R & Frackowiak RSJ (1995) Characterising evoked ',...
  'hemodynamics with fMRI Friston KJ, NeuroImage 2:157-165'];
r5 = [...
  '  * Josephs O, Turner R and Friston KJ (1997) Event-related fMRI, Hum. Brain ',...
  'Map. 0:00-00'];
conf.help = {'fMRI Statistics','',p1,'',p2,'',p3,'',p4,'',p5,...
             '','Referencess:',r1,'',r2,'',r3,'',r4,'',r5};

return;
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------

function t = check_dir(job)
t = {};
%d = pwd;
%try,
%    cd(job.dir{1});
%catch,
%    t = {['Cannot Change to directory "' job.dir{1} '".']};
%end;

%disp('Checking...');
%disp(fullfile(job.dir{1},'SPM.mat'));

% Should really include a check for a "virtual" SPM.mat
% if exist(fullfile(job.dir{1},'SPM.mat'),'file'),
%     t = {'SPM files exist in the analysis directory.'};
% end;
return;
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
function my_cd(varargin)
job = varargin{1};
if ~isempty(job)
    try
    cd(char(job));
    fprintf('Changing directory to: %s\n',char(job));
    catch
        error('Failed to change directory. Aborting run.')
    end
end
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function run_est(job)
% Set up the design matrix and run a design.

spm_defaults;
global defaults
defaults.modality='FMRI';

% Load SPM.mat file
%------------------------------------------------------------
SPM=[];
load(job.spmmat{:});

original_dir = pwd;
my_cd(job.dir);

%-Ask about overwriting files from previous analyses...
%-------------------------------------------------------------------
if exist(fullfile(job.dir{1},'SPM.mat'),'file')
    str = {	'Current directory contains existing SPM file:',...
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
files = {'mask.???','ResMS.???','RVP.???',...
    'beta_????.???','con_????.???','ResI_????.???',...
    'ess_????.???', 'spm?_????.???'};
for i=1:numel(files),
    if any(files{i} == '*' | files{i} == '?')
        [j,unused] = spm_list_files(pwd,files{i});
        for i=1:size(j,1),
            spm_unlink(deblank(j(i,:)));
        end
    else
        spm_unlink(files{i});
    end
end

% Bayesian 2nd level estimation
if isfield(job.method,'Bayesian2')
    SPM=spm_spm_Bayes(SPM);
    my_cd(original_dir); % Change back
    fprintf('Done\n')
    return
end

% REML estimation
if isfield(job.method,'Classical'),
    SPM.xVi.form = job.method.Classical.cvi;
    SPM = spm_spm(SPM);
    
    % Automatically set up contrasts for factorial designs
    if isfield(SPM,'factor')
        cons=spm_design_contrasts(SPM);
        
        % Create F-contrasts
        for i=1:length(cons),
            con=cons(i).c;
            name=cons(i).name;
            STAT='F';
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
            end;
            spm_contrasts(SPM,length(SPM.xCon));
        end
        
        % Create t-contrasts
        for i=1:length(cons),
            % Create a t-contrast for each row of each F-contrast
            % The resulting contrast image can be used in a 2nd-level analysis
            Fcon=cons(i).c;
            nrows=size(Fcon,1);
            STAT='T';
            for r=1:nrows,
                con=Fcon(r,:);     
                
                % Change name 
                str=cons(i).name;
                sp1=min(find(str==' '));
                if strcmp(str(1:11),'Interaction')
                    name=['Positive ',str,'_',int2str(r)];
                else
                    name=['Positive',str(sp1:end),'_',int2str(r)];
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
                end;
                spm_contrasts(SPM,length(SPM.xCon));
            end
        end
    end
    
    my_cd(original_dir); % Change back
    fprintf('Done\n')
    return
end

% Bayesian 1st level estimation options

% Analyse specific slices or whole volume
if isfield(job.method.Bayesian.space,'Slices')
    SPM.PPM.space_type='Slices';
    SPM.PPM.AN_slices=job.method.Bayesian.space.Slices;
else
    SPM.PPM.space_type='Volume';
end

% Regression coefficient priors
switch job.method.Bayesian.signal
    case 'GMRF',
        SPM.PPM.priors.W='Spatial - GMRF';
    case 'LORETA',
        SPM.PPM.priors.W='Spatial - LORETA';
    case 'Global',
        SPM.PPM.priors.W='Voxel - Shrinkage';
    case 'Uninformative',
        SPM.PPM.priors.W='Voxel - Uninformative';
    otherwise
        disp('Unkown prior for W in spm_config_fmri_stats');
end

% Number of AR coefficients
SPM.PPM.AR_P=job.method.Bayesian.ARP;

% AR coefficient priors
if isfield(job.method.Bayesian.noise,'GMRF')
    SPM.PPM.priors.W='Spatial - GMRF';
elseif isfield(job.method.Bayesian.noise,'LORETA')
    SPM.PPM.priors.W='Spatial - LORETA';
elseif isfield(job.method.Bayesian.noise,'tissue_type')
    SPM.PPM.priors.W='Discrete';
    SPM.PPM.priors.SY=job.method.Bayesian.noise.tissue_type;
end

% Define an empty contrast
NullCon.name=[];
NullCon.c=[];
NullCon.STAT = 'P';
NullCon.X0=[];
NullCon.iX0=[];
NullCon.X1o=[];
NullCon.eidf=1;
NullCon.Vcon=[];
NullCon.Vspm=[];

SPM.xCon=[];
% Set up contrasts for 2nd-level ANOVA
if strcmp(job.method.Bayesian.anova.second,'Yes')
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

% Set up user-specified simple contrasts
ncon=length(job.method.Bayesian.gcon);
K=size(SPM.xX.X,2);
for c = 1:ncon,
    DxCon=NullCon;
    DxCon.name = job.method.Bayesian.gcon(c).name;
    convec=sscanf(job.method.Bayesian.gcon(c).convec,'%f');
    if length(convec)==K
        DxCon.c = convec;
    else
        disp('Error in spm_config_fmri_est: contrast does not match design');
        return
    end
    
    if isempty(SPM.xCon),
        SPM.xCon = DxCon;
    else
        SPM.xCon(end+1) = DxCon;
    end;
end
    
% 1st level Bayesian ANOVA ?
bayes_anova=0;
if strcmp(job.method.Bayesian.anova.first,'Yes')
    bayes_anova=1;
    SPM.PPM.update_F=1; % Compute evidence for each model
    SPM.PPM.compute_det_D=1; 
end
SPM = spm_spm_vb(SPM);
if bayes_anova
    % We don't want to estimate contrasts for each different model
    SPM.xCon=[];
    spm_vb_ppm_anova(SPM);
end

my_cd(original_dir); % Change back
fprintf('Done\n')
return


%-------------------------------------------------------------------------
function vf = vfiles_stats(job)
direc = job.dir{1};
vf    = {fullfile(direc,'SPM.mat')};

% Should really create a few vfiles for beta images etc here as well.

