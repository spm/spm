# Output Module Command List
This module is responsible for generating source-level summary images, or exporting source reconstucted timeseries for futher analysis.

#### Top-level commands
- [BF](#BF)

#### Semi-universal inputs (DAiSS-wizard only!)
- [conditions](#semi_conditions)
- [contrast](#semi_contrast)
- [foi](#semi_foi)
- [woi](#semi_woi)

#### Supported Output Plugins
- [image_dics](#image_dics)
- image_mv
- image_power

## Commands
#### BF
Path to BF.mat file

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: cellstr
matlabbatch{1}.spm.tools.beamforming.output.BF = {'/path/to/BF/'};

% DAiSS-Wizard
% Default: REQUIRED
% Input Type: str or cellstr
S.BF = 'path/to/BF';
```

## Semi-universal inputs (DAiSS-Wizard only)
These inputs are not required by every output method, but it sometimes help that many of the common ones can be called simply.

#### conditions<a name="semi_conditions"></a>
Specify which trials you'd like to be included when generating the output. Can be either all trials or a subset of conditions
##### Case: All Trials
```matlab
% DAiSS-Wizard
% Default: 'all'
% Input Type: str or cellstr
S.conditions = 'all';
```
##### Case: Specific conditions
```matlab
% DAiSS-Wizard
% Default: 'all'
% Input Type: cellstr
S.conditions = {'condtion_01','condition_02');
```

#### contrast<a name="semi_contrast"></a>
Specify the contrast between conditions
```matlab
% DAiSS-Wizard
% Default: 1
% Input Type: numeric
S.conditions = [1 -1]; % subtract second condition from first
```

#### foi<a name="semi_foi"></a>
Which requency band(s) of interest do you want to filter in? <img src="https://render.githubusercontent.com/render/math?math=n_{bands} \times 2"> to specify multiple windows. 
```matlab
% DAiSS-Wizard
% Default: [0 Inf]
% Input Type: numeric
S.foi = [13 30];
```

#### woi<a name="semi_woi"></a>
Which window(s) of interest in each trial do you want to use? Can be an n x 2 matrix to specify multiple windows. Specify in milliseconds.
```matlab
% DAiSS-Wizard
% Default: [-Inf Inf]
% Input Type: numeric
S.woi = [0 1000]; % for the first second of each trial
```

## Inverse plugins
### image_dics
Dynamic Imaging of Coherent Source (DICS; [Gross et al. (2001)](https://doi.org/10.1073/pnas.98.2.694)) is used to measure the coherence between reference signal (such as an EMG) and the brain, or simply used to generate power images.

List of image_dics commands
- reference
- powmethod
- whatconditions
- sametrials
- woi
- contrast
- logpower
- foi
- taper
- result
- scale
- modality

#### whatconditions
Specify which trials you'd like to be included when generating the features matrix. Can be either all trials or a subset of conditions
##### Case: All Trials
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: logical
matlabbatch{1}.spm.tools.beamforming.output.plugin.image_dics.whatconditions.all = 1;

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
matlabbatch{1}.spm.tools.beamforming.output.plugin.image_dics.whatconditions.conditions = {'condtion_01','condition_02');

% DAiSS-Wizard
% Default: 'all'
% Input Type: cellstr
S.conditions = {'condtion_01','condition_02');
```

#### sametrials
Take the same trials as used for filter computation. This is useful for bootstraping.
```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: logical
matlabbatch{1}.spm.tools.beamforming.output.plugin.image_dics.sametrials = true;

% matlabbatch
% Default: false
% Input Type: logical
S.method = 'image_dics';
S.image_dics.sametrials = true;
```
