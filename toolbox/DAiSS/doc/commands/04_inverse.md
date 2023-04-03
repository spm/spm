# Inverse Module Command List
This module is responsible for generating the source reconstruction weights for all sources defined.

#### Top-level commands
- [BF](#BF)

#### Plugins for source reconstruction

- [champagne](#champagne)
- deflect
- [dics](#dics)
- [ebb](#ebb)
- [eloreta](#eloreta)
- [lcmv](#lcmv)
- lcmv_multicov
- [minimumnorm](#minimumnorm)
- [nutmeg](#nutmeg) (for dSPM, sLORETA etc.)

## Commands
#### BF
Path to BF.mat file

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: cellstr
matlabbatch{1}.spm.tools.beamforming.inverse.BF = {'/path/to/BF/'};

% DAiSS-Wizard
% Default: REQUIRED
% Input Type: str or cellstr
S.BF = 'path/to/BF';
```

## Inverse plugins
### champagne

A Bayesian source reconstuction method. Used in conjucntion with [VBFA features method](03_features.md#vbfa) See [Wipf et al. (2010)](https://doi.org/10.1016/j.neuroimage.2009.06.083) or [Owen et al. (2012)](https://doi.org/10.1016/j.neuroimage.2011.12.027) for more information. 

#### nem

Number of iterations of EM used.

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: numeric
matlabbatch{1}.spm.tools.beamforming.inverse.plugin.champagne.nem = 100;

% DAiSS-Wizard
% Default: 100
% Input Type: numeric
S.method = 'champagne'
S.champagne.nem = 100;
```

#### vcs

Define the voxel covariance structure. specify a value in the range of 0-2;

Options:
- [0] Scalar
- [1] Diagonal
- [2] General *(Recommended)*

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: numeric
matlabbatch{1}.spm.tools.beamforming.inverse.plugin.champagne.vcs = 2;

% DAiSS-Wizard
% Default: 2
% Input Type: numeric
S.method = 'champagne'
S.champagne.vcs = 2;
```


#### nupd
How is the noise covariance estimate, spefcify a value in the range of 0-2

Options:
- [0] Use provided from [VBFA](03_features.md#vbfa) *(Recommended)*
- [1] Learn scalar
- [2] Learn diagonal

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: numeric
matlabbatch{1}.spm.tools.beamforming.inverse.plugin.champagne.nupd = 0;

% DAiSS-Wizard
% Default: 0
% Input Type: numeric
S.method = 'champagne'
S.champagne.nupd = 0;
```

***

### deflect

***

### dics

Dynamic Imaging of Coherent Source (DICS; [Gross et al. (2001)](https://doi.org/10.1073/pnas.98.2.694)) is a spatial filtering method similar to LCMV beamforming, but uses a cross-spectral density matrix in its formulation rather than a covariance. **Note:** this should be paired with the [csd](03_features.md#csd) features method.

#### fixedori

Select whether you want the dipole orientation optimised to give maximal power. Inputs are either strings **yes** or **no**.

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: string
matlabbatch{1}.spm.tools.beamforming.inverse.plugin.dics.fixedori = 'yes';

% DAiSS-Wizard
% Default: 'yes'
% Input Type: string
S.method = 'dics'
S.dics.fixedori = 'yes';
```

***

### ebb

Empirical Bayesian Beamformer. Uses beamformer-like assumptions of source variance in the generative model prior. See [Belardinelli et al. (2012)](https://doi.org/10.1371/journal.pone.0051985) for more information about the original variant of EBB and [O'Neill et al. (2021)](https://doi.org/10.1038/s41598-021-96933-0) for information about correlated EBB.

**WARNING 1:** the version of EBB here will not give numerically identical results to the version of the EBB routine in [spm_eeg_invert](https://github.com/spm/spm/blob/main/spm_eeg_invert.m) owing to some differences in preprocessing. 

**WARNING 2:** use the [tdcov](03_features.md#tdcov) features method in conjuntion of this for the most useful model evidence results.

EBB has a lot of options. Most should default to a simple behaviour in both matlabbatch and with DAiSS-Wizard.

- [keeplf](#ebb_keeplf)
- [iid](#ebb_iid)
- [corr](#ebb_corr)
- [onlycorr](#ebb_onlycorr)
- [diags](#ebb_diags)
- [mixmethod](#ebb_mixmethod)
- [pairs](#ebb_pairs)
- [reml](#ebb_reml)
- [noise](#ebb_noise)


#### keeplf <a name="ebb_keeplf"></a>

Do you want to keep the reduced lead fields? (Might be useful if a 3D lead field vector is provided, as this will be reduced optimised to 1D with EBB).
```matlab

% matlabbatch
% Default: false
% Input Type: logical
matlabbatch{1}.spm.tools.beamforming.inverse.plugin.ebb.keeplf = 0;

% DAiSS-Wizard
% Default: false
% Input Type: logical
S.method = 'ebb'
S.ebb.keeplf = 0;
```
#### iid <a name="ebb_iid"></a>
Set the source covariance matrix to be an identity matrix. This makes it act like an adaptive classic minimum-norm solution.
```matlab

% matlabbatch
% Default: false
% Input Type: logical
matlabbatch{1}.spm.tools.beamforming.inverse.plugin.ebb.iid = 0;

% DAiSS-Wizard
% Default: false
% Input Type: logical
S.method = 'ebb'
S.ebb.iid = 0;
```
#### corr <a name="ebb_corr"></a>
Adds correlated source assumptions to the prior in addition to the uncorrelated source assmptions of EBB, this is the implementation of cEBB used in [O'Neill et al. (2021)](https://doi.org/10.1038/s41598-021-96933-0).
Note: this automatically assumes the correlated sources are bilateral (i.e. mirrored along the central fissure). For a more flexible implementation of correlated sources see the [pairs](#ebb_pairs) option.
```matlab

% matlabbatch
% Default: false
% Input Type: logical
matlabbatch{1}.spm.tools.beamforming.inverse.plugin.ebb.corr = 0;

% DAiSS-Wizard
% Default: false
% Input Type: logical
S.method = 'ebb'
S.ebb.corr = 0;
```
#### onlycorr <a name="ebb_onlycorr"></a>
Similar to [corr](#ebb_corr) but the uncorrelated source priors are not included.
```matlab

% matlabbatch
% Default: false
% Input Type: logical
matlabbatch{1}.spm.tools.beamforming.inverse.plugin.ebb.onlycorr = 0;

% DAiSS-Wizard
% Default: false
% Input Type: logical
S.method = 'ebb'
S.ebb.onlycorr = 0;
```

#### diags <a name="ebb_diags"></a>
For correlated sources, where does the variance get populated on the source matrix?

Options:
- **off**: using the off-diagonals is similar to how Multiple Sparse Priors models correlated sources, this does force the singals to be cross-projected.
- **on**: (default) on-diagonals do not force signals from both sources to be cross-projected, but rather adjusts lost variance from uncorrelated assumptions.
- **both**: life is short, choose everything.

**Warning:** the use of off-diagonals, seems to always to give a better model evidence score when compared to an uncorrelated source model, EVEN IF there all the sources are truly uncorrelated. It's behaviour is similar to [Cheesoid thinking everything smells of cheese](https://www.youtube.com/watch?v=6Ofesgkhikg&t=190s). (Don't worry, this doesn't happen the Multiple Spare Priors - this is just a quirk of the free energy calculation with large monolithic priors).

```matlab

% matlabbatch
% Default: 'on'
% Input Type: str
matlabbatch{1}.spm.tools.beamforming.inverse.plugin.ebb.diags = 'on';

% DAiSS-Wizard
% Default: 'on'
% Input Type: str
S.method = 'ebb'
S.ebb.diags = 'on';
```

#### mixmethod <a name="ebb_mixmethod"></a>
How do you want the uncorrelated and correlated source priors to be combined together?

Options:
- **sum**: (default) this simply adds the priors together to have one source covariance matrix.
- **reml**: keeps the two priors seperate and scales them relative to each other when passed to ReML. This is largely untested re. whether this is useful or not.

```matlab

% matlabbatch
% Default: 'sum'
% Input Type: str
matlabbatch{1}.spm.tools.beamforming.inverse.plugin.ebb.mixmethod = 'sum';

% DAiSS-Wizard
% Default: ''
% Input Type: str
S.method = 'ebb'
S.ebb.mixmethod = 'sum';
```

#### pairs <a name="ebb_pairs"></a>
Path to a mat file containing a binary adjacency matrix of correlated pairs. If a matrix is not supplied, it will automatically look for a the homologous regions.
**TIP:** If you want a source to not be correlated, pair it with itself (i.e. put a 1 on the diagonal element.

```matlab

% Example: in a simple 10-source model, set sources 5 and 8 to be correlated;
A = eye(10); 
A(5,5) = 0;
A(8,8) = 0;
A(5,8) = 1;
A(8,5) = 1;

save(A,'foo.mat') % variable name not important, as long as only variable in mat file!

% matlabbatch
% Default: ''
% Input Type: str
matlabbatch{1}.spm.tools.beamforming.inverse.plugin.ebb.corr = 1;
matlabbatch{1}.spm.tools.beamforming.inverse.plugin.ebb.pairs = 'foo.mat';

% DAiSS-Wizard
% Default: ''
% Input Type: str
S.method = 'ebb'
S.ebb.corr = 1;
S.ebb.pairs - 'foo.mat';
```

#### reml <a name="ebb_reml"></a>
EBB uses a ReML optimsation to do the final mixing of all source/noise priors and returns the corresponding model evidence score. Information as to how ReML works here can be found in [REF]. How wide do we want to case the hyperparameter net?

Options:
- **Loose**: (default) allows the hyperparameters to freely scale how they want with little penalty. Similar behaviour to how the core of SPM uses ReML.
- **Strict**: Ensures that the noise and source covariance matricies are scaled in a ~1:100 ratio. Might be useful if the strict option returns odd results.

```matlab

% matlabbatch
% Default: 'Loose'
% Input Type: str
matlabbatch{1}.spm.tools.beamforming.inverse.plugin.ebb.reml = 'Loose';

% DAiSS-Wizard
% Default: 'Loose'
% Input Type: str
S.method = 'ebb'
S.ebb.reml = 'loose';
```

#### noise <a name="ebb_noise"></a>
Add a noise matrix to the generative model. 

The noise assumption is that sensor noise is i.i.d but here a DAiSS BF.mat file with a noise covariance matrix can be supplied.

```matlab

% matlabbatch
% Default: ''
% Input Type: str
matlabbatch{1}.spm.tools.beamforming.inverse.plugin.ebb.noise = '/path/to/noise/BF.mat';

% DAiSS-Wizard
% Default: ''
% Input Type: str
S.method = 'ebb'
S.ebb.noise = '/path/to/noise/BF.mat';
```

***

### eloreta
eLORETA is a variant of minimum-norm, which boasts a zero-error dipole localisation. See [Pascual-Maqui et al. (2007)](http://arxiv.org/pdf/0710.3341) for more details.

#### regularisation

Optional regularization parameter (example: 0.05 corresponds to 5% of the average of the eigenvalues of some matrix to be inverted).

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: numeric
matlabbatch{1}.spm.tools.beamforming.inverse.plugin.ebb.regularisation = 0.05;

% DAiSS-Wizard
% Default: 0.05
% Input Type: numeric
S.method = 'eloreta'
S.ebb.regularisation = 0.05;
```

***

### lcmv

Linear Constrained Minimal Variance (LCMV) Beamforming is a spatial filtering method of source reconstruction. It can be quite useful as either a general purpose source reconsruction method, or in cases where you think your sensor-level SNR is not that high, due to being a fairly robust at rejecting non-neural interference. For more information please see either [Brookes et al (2008)](https://doi.org/10.1016/j.neuroimage.2007.09.050) or [Westner et al (2022)](https://doi.org/10.1016/j.neuroimage.2021.118789). 

#### orient

Select whether you want the dipole orientation optimised to give maximal power. Optimisation based on [Sekihara's eigenvalue decomposition method](https://doi.org/10.1109/TBME.2004.827926).

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: logical
matlabbatch{1}.spm.tools.beamforming.inverse.plugin.lcmv.orient = true;

% DAiSS-Wizard
% Default: true
% Input Type: logical
S.method = 'lcmv'
S.lcmv.orient = true;
```

#### keeplf

Do you want to keep the oriented lead fields if the orient option above was set to true. 

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: logical
matlabbatch{1}.spm.tools.beamforming.inverse.plugin.lcmv.orient = true;

% DAiSS-Wizard
% Default: false
% Input Type: logical
S.method = 'lcmv'
S.lcmv.keeplf = false;
```

***

### lcmv_multicov

***

### minimumnorm

Minimum norm estimator based on the [DeFleCT](https://imaging.mrc-cbu.cam.ac.uk/meg/AnalyzingData/DeFleCT_SpatialFiltering_Tools) tools developed by Hauk and Stenroos. See [Hauk and Stenroos (2014)](https://doi.org/10.1002/hbm.22279) for more informaiton.

#### snr

The assumed ratio of variances of signal and noise, used for setting the regularisation parameter.

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: numeric
matlabbatch{1}.spm.tools.beamforming.inverse.plugin.minimumnorm.snr = 5;

% DAiSS-Wizard
% Default: 5
% Input Type: numeric
S.method = 'minimumnorm'
S.minimumnorm.snr = 5;
```

#### trunc

The number of (smallest) singular values of the covariance matrix that are set to zero before making the whitener. For example, if the data has been SSP-projected, it needs to be at least the number of components projected away. **Note:** This flag is useful if you have provided an [unregularised covariance](03_features.md#none) but know that its going to be rank deficient.

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: numeric
matlabbatch{1}.spm.tools.beamforming.inverse.plugin.minimumnorm.trunc = 0;

% DAiSS-Wizard
% Default: 0
% Input Type: numeric
S.method = 'minimumnorm'
S.minimumnorm.trunc = 0;
```

***

### nutmeg

A selection of weigted minumim norm solutions implemented in [NUTMEG](https://www.nitrc.org/plugins/mwiki/index.php/nutmeg:MainPage). 

#### method

Which flavour of weighted minimum norm do you want to use?

Options:
- **dSPM** based on [Dale et al (2000)](https://doi.org/10.1016/S0896-6273(00)81138-1)
- **sLORETA** based on *paper tbc*
- **swLORETA** based on [Palmero-Soler et al (2007)](https://doi.org/10.1088/0031-9155/52/7/002)
