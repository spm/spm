# Features Module Command List
This module is responsible for for the preparation of a features matrix (normally a covariance) required as a part of an inverse model for source reconstruction.
This also allows us to inspect the eigenspectrum of the matrix and regularise to allow for more stable solutions.

#### Top-level commands
- [BF](#BF)
- [whatconditions](#whatconditions)
- [woi](#woi)
- [modality](#modality)
- [fuse](#fuse)
- [cross_terms](#cross_terms)
- [bootstrap](#cross_terms)
- [visualise](#visualise)

#### Plugins for matrix generation

- [contcov](#contcov)
- [cov](#cov)
- [cov_bysamples](#cov_bysamples)
- [csd](#csd)
- [identity](#identity)
- regmulticov
- [tdcov](#tdcov)
- [vbfa](#vbfa)

#### Plugins for regularisation

- [none](#none)
- [clifftrunc](#clifftrunc)
- [mantrunc](#mantrunc)
- [manual](#manual) (classic Tikhonov regrularisation)
- [minkatrunc](#minkatrunc)
- roi
- tikhonov_rankdef

## Commands
#### BF
Path to BF.mat file

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: cellstr
matlabbatch{1}.spm.tools.beamforming.features.BF = {'/path/to/BF/'};

% DAiSS-Wizard
% Default: REQUIRED
% Input Type: str or cellstr
S.BF = 'path/to/BF';
```
#### whatconditions
Specify which trials you'd like to be included when generating the features matrix. Can be either all trials or a subset of conditions
##### Case: All Trials
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: logical
matlabbatch{1}.spm.tools.beamforming.features.whatconditions.all = 1;

% DAiSS-Wizard
% Default: 'all'
% Input Type: str or cellstr
S.conditions = 'all';
```
##### Case: Specific conditions
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: cellstr
matlabbatch{1}.spm.tools.beamforming.features.whatconditions.conditions = {'condtion_01','condition_02');

% DAiSS-Wizard
% Default: 'all'
% Input Type: cellstr
S.conditions = {'condtion_01','condition_02');
```

#### woi
Which window(s) of interest in each trial do you want to use? Can be an n x 2 matrix to specify multiple windows. Specify in milliseconds.
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: numeric
matlabbatch{1}.spm.tools.beamforming.features.woi = [0 1000]; % for the first second of each trial

% DAiSS-Wizard
% Default: [-Inf Inf]
% Input Type: numeric
S.woi = [0 1000]; % for the first second of each trial
```

#### modality
Which modality/modalities do you want to build a features matrix with? Combinw the options below in a cell array to specify multiple modalities.

Options:
- MEG
- MEGPLANAR
- EEG

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: cellstr
matlabbatch{1}.spm.tools.beamforming.features.modality = {'MEG','EEG'}; % for concurrent MEG/EEG


% DAiSS-Wizard
% Default: {'MEG'}
% Input Type: str or cellstr
S.modality = {'EEG'}; % for just EEG sensors
```

#### fuse
Fuse sensors from different modalities together (assume scaling/pre-whitening has already been performed prior to DAiSS).

Options:
- **no**.       Do not combine. *default for DAiSS-Wizard*
- **meg**.      Use this for Elekta/MEGIN systems where planar gradiometers and magnetometers are both being used.
- **all**.      Comvine all in the selected modalities together!

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: str
matlabbatch{1}.spm.tools.beamforming.features.fuse = 'no'; % for concurrent MEG/EEG


% DAiSS-Wizard
% Default: 'no'
% Input Type: str
S.fuse = 'no'; % for just EEG sensors
```

#### cross_terms
If fusing modalities, set the cross terms between modalities in the matrix to zero.

Options:
- **meegeeg**.      Set terms between MEG and EEG to zero. (Default in DAiSS-Wizard)
- **all**.          Set only different kinds of MEG sensors to zero (might be useful for MEGIN/Elekta)
- **no**.           Don't set cross-terms to zero

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: str
matlabbatch{1}.spm.tools.beamforming.features.cross_terms = 'meegeeg'; 


% DAiSS-Wizard
% Default: 'megeeg'
% Input Type: str
S.cross_terms = 'megeeg'; % for just EEG sensors
```

#### bootstrap
Generate the matrix using a bootstap selection of trials. e.g. instead of trials `[1 2 3 4 5 6 7 8 9]` you would make a matrix out of a random subset of the trials, but keeping the total number of trials the same, such as `[1 1 1 4 5 6 7 7 9]`. 

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: logical
matlabbatch{1}.spm.tools.beamforming.features.bootstrap = false; 

% DAiSS-Wizard
% Default: false
% Input Type: logical
S.bootstrap = false;
```

#### visualise
Visualise the eignespectrum of the generated covariance matrix, can be useful for identifying the number of components to cut/keep when regularising.

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: logical
matlabbatch{1}.spm.tools.beamforming.features.visualise = true; 

% DAiSS-Wizard
% Default: true
% Input Type: logical
S.visualise = true;
```

## Matrix generation plugins
### contcov
A.K.A robust covariance. Performs the operation $C=\frac{1}{n}YY^T$. This function assumes the data has been filtered prior to this step. No extra flags to call as such but need to explicitly call the method. Can work for epoched or continuous data.

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: N/A
matlabbatch{1}.spm.tools.beamforming.features.plugins.contcov = {}; 

% DAiSS-Wizard
% Default: N/A
% Input Type: N/A
S.method = 'contcov';
```

***

### cov
Band-pass filtered covariance matrix. Data is filtered using a discrete cosine transform (DCT) prior to covarianc calculation. **Warning:** this method will only work on epoched data due to the $n_{samples} \times n_{samples}$ matrix required to do the filtering using a DCT.

#### foi
Which requency band(s) of interest do you want to filter in? $n_{bands} \times 2$ to specify multiple windows. 
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: numeric
matlabbatch{1}.spm.tools.features.features.plugin.cov.foi = [13 30]; % for the first second of each trial

% DAiSS-Wizard
% Default: [0 Inf]
% Input Type: numeric
S.method = 'cov'; % for the first second of each trial
S.cov.foi = [13 30];
```

#### taper
Apply a hanning window to minimise the effects of edge-of-epoch effects.
Options:
- **hanning**: Apply the hanning window
- **none**: do nothing
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: str
matlabbatch{1}.spm.tools.features.features.plugin.cov.taper = 'hanning'; % for the first second of each trial

% DAiSS-Wizard
% Default: 'hanning'
% Input Type: numeric
S.method = 'cov'; % for the first second of each trial
S.cov.taper = 'hanning';
```

***

### cov_bysamples
Similar to the Robust Covariance ([contcov](#contcov)) module, but a binary mask is used to speficy which time windows are used. (Useful if you have states allocated from e.g. a HMM where the windows of activation are all different lengths). There are two ways to call this, either by using the **samples** command or by having a channel called 'Class' in your meeg object which can have multiple states encoded within in (assumes mutual exclusivity).

#### samples
A $1 \times n_{samples per trial} \times n_{trials}$ binary array to select which time windows.
  
##### Example calling from matlabbatch/wizard with an explicit variable.
```matlab

% matlabbatch
% Default: N/A
% Input Type: numeric
matlabbatch{1}.spm.tools.features.features.plugin.cov_bysamples.samples = samples; % 1 x t x ntrials

% DAiSS-Wizard
% Default: N/A
% Input Type: numeric
S.method = 'cov_bysamples'; % for the first second of each trial
S.cov_bysamples.samples = samples;
```

##### Example: using a channel called Class
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: numeric
matlabbatch{1}.spm.tools.features.features.plugin.cov_bysamples = {[]}; 

% DAiSS-Wizard
% Default: REQUIRED
% Input Type: numeric
S.method = 'cov_bysamples'; % for the first second of each trial
```

***

### csd
Cross-spectral density matrix. Utilises 'ft_specest_mtmfft' to perform spectral analysis
#### foi 
Specify the frequency band of interest for the CSD
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: numeric (1 x 2 array)
matlabbatch{1}.spm.tools.features.features.plugin.csd.foi = [13 30]; 

% DAiSS-Wizard
% Default: REQUIRED
% Input Type: numeric (1 x 2 array)
S.method = 'csd'; % for the first second of each trial
S.csd.foi = [13 30];
```
#### taper 
Select which kind of taper to convolve with data during spectral analysis

Options:
- hanning
- rectwin
- dpss
- sine

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: str
matlabbatch{1}.spm.tools.features.features.plugin.csd.taper = 'dpss'; 

% DAiSS-Wizard
% Default: 'dpss'
% Input Type: str
S.method = 'csd'; % for the first second of each trial
S.csd.taper = 'dpss';
```
#### keepreal 
Do we only want to keep the real part of the CSD, setting to false keeps the CSD complex.
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: logical
matlabbatch{1}.spm.tools.features.features.plugin.csd.keepreal = false; 

% DAiSS-Wizard
% Default: false
% Input Type: logical
S.method = 'csd'; % for the first second of each trial
S.csd.keepreal = false;
```

#### han
Hanning window trials to minimise edge effects
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: logical
matlabbatch{1}.spm.tools.features.features.plugin.csd.keepreal = false; 

% DAiSS-Wizard
% Default: true
% Input Type: logical
S.method = 'csd'; % for the first second of each trial
S.csd.han = true;
```

***

### identity
Returns an identity matrix. No options to call.
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: logical
matlabbatch{1}.spm.tools.features.features.plugin.identity = {[]};

% DAiSS-Wizard
% Default: true
% Input Type: logical
S.method = 'identity'; % for the first second of each trial
```

***

### regmulticov

***

### tdcov

Band-pass filtered covariance generation, with an additional temporal decomposition step. 
Used for the source inversions which use the Fristonian version of Free Energy (eg. EBB). See [Lopez et al. (2014)](https://doi.org/10.1016/j.neuroimage.2013.09.002) for more information.

#### foi
Specify the frequency band of interest to filter data prior to covariance generation.

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: numeric
matlabbatch{1}.spm.tools.features.features.plugin.tdcov.foi = [0 Inf];

% DAiSS-Wizard
% Default: [0 Inf]
% Input Type: numeric
S.method = 'tdcov'; % for the first second of each trial
S.todcov.foi = [0 inf];
```

#### ntmodes
Number of temporal modes to decompose each trial down to. Set to 0 for an automatic determination.
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: numeric
matlabbatch{1}.spm.tools.features.features.plugin.tdcov.ntmodes = 4;

% DAiSS-Wizard
% Default: 4
% Input Type: numeric
S.method = 'tdcov'; % for the first second of each trial
S.todcov.ntmodes = 4;
```

#### taper
Apply a hanning window to minimise the effects of edge-of-epoch effects.
Options:
- **hanning**: Apply the hanning window
- **none**: do nothing
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: str
matlabbatch{1}.spm.tools.features.features.plugin.tdcov.taper = 'hanning'; 

% DAiSS-Wizard
% Default: 'hanning'
% Input Type: str
S.method = 'tdcov'; % for the first second of each trial
S.tdcov.taper = 'hanning';
```

***

### vbfa
Variational Bayes Factorial Analysis. Used in conjuction with the Champagne source inversion to estimate the noise covariance. See [Wipf et al. (2010)](https://doi.org/10.1016/j.neuroimage.2009.06.083) for more information.

#### nl
Factor dimensionality
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: numeric
matlabbatch{1}.spm.tools.features.features.plugin.vbfa.nl = 5;

% DAiSS-Wizard
% Default: 5
% Input Type: numeric
S.method = 'vbfa'; 
S.vbfa.taper = 5;
```

#### nem
Number of iterations of EM-optimisation to find factors.
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: numeric
matlabbatch{1}.spm.tools.features.features.plugin.vbfa.nl = 50;

% DAiSS-Wizard
% Default: 50
% Input Type: numeric
S.method = 'vbfa'; 
S.vbfa.taper = 50;
```

## Regularisation plugins

More often than not the features matrix has to be inverted. However there are various reasons where the matrix itself is rank deficient (i.e. there is less content in the matrix than its size implies). Inverting a rank deifcient matrix can lead to numerical errors which result in a poor source reconstuction. The methods below attempt to minimise the effect of this.

### none

Ignore everything I just said above, this does nothing! Note DAiSS doesn't explicity have a *no regularisation* function, but the wizard scripts do.
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: numeric
matlabbatch{1}.spm.tools.features.features.regularisation.manual.lambda = 0;

% DAiSS-Wizard
% Default: 'none'
% Input Type: str
S.reg = 'none';
```

***

### clifftrunc

Autmatically determines the number of principal components in a matrix to keep. 
Searches for the first large jump between eigenvalues in log-space, where it appears the eigenvalues have dropped of a cliff. See [Westerner et al. (2022); Fig. 1](https://www.sciencedirect.com/science/article/pii/S1053811921010612?via%3Dihub#fig0001) for examples of where the rank is estimated. Works well for data where a rank-deficiency is expected due to a regression (e.g. SSS, tSSS, SSP, ICA etc.)

#### zscore
The jumps in eigenvalue differences are z-transformed prior to finding the first large jump. Set the critical z-score to stop search at. Set to a value less than 0 to find the greatest jump.
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: numeric
matlabbatch{1}.spm.tools.features.features.regularisation.clifftrunc.zscore = -1;

% DAiSS-Wizard
% Default: -1
% Input Type: numeric
S.reg = 'clifftrunc';
S.clifftrunc.zscore = -1;
```

***

### mantrunc

Manually select the number of principal components to keep. 

#### pcadim
Number of components to keep.
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: numeric
matlabbatch{1}.spm.tools.features.features.regularisation.mantrunc.pcadim = 100;

% DAiSS-Wizard
% Default: 100
% Input Type: numeric
S.reg = 'mantrunc';
S.mantrunc.pcadim = 100;
```

***

### manual
The classic Tikhonov regularisation, where a an identity matrix is added to the matrix to reduce the condition of the matrix. In this case $C_r = C + \mu I$ where $$\mu=\frac{\lambda}{n}\sum_{i=1}^{n}C_{ii}$$

#### lambda
The fractional proportion of the average of the diagonal elements (see equation above).
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: numeric
matlabbatch{1}.spm.tools.features.features.regularisation.manual.lambda = 0.05;

% DAiSS-Wizard
% Default: 0
% Input Type: numeric
S.reg = 'manual';
S.manual.lambda = 0.05; % 5% regularisation
```

***

### minkatrunc
Identifies the rank of the covariance matrix using Bayesian assumptions. See [Minka (2000)](https://proceedings.neurips.cc/paper/2000/file/7503cfacd12053d309b6bed5c89de212-Paper.pdf) for more details. Works well for SQUID-MEG data, but we've noticed that for OPM data which has undergone any methods which reduce the data's rank (such as [HFC](https://doi.org/10.1016/j.neuroimage.2021.118484)) it fails to correctly identify the rank. Use with caution.

#### reduce
Reduces the size of the matrix to match the rank of the data. i.e. A 200 x 200 features matrix with rank 50 will be reduced to size 50 x 50. (Recommended).
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: logical
matlabbatch{1}.spm.tools.features.features.regularisation.minktatrunc.reduce = 1;

% DAiSS-Wizard
% Default: 1
% Input Type: logical
S.reg = 'minkatrunc';
```
