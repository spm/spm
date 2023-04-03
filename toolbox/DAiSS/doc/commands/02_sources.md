# Sources Module Command List

The source module generates the source space used for forward/inverse calculations and calculated the forward solutions. 


#### List of available top-level commands
- [BF](#BF)
- [reduce_rank](#reduce_rank)
- [keep3d](#keep3d)
- [normalise_lf](#normalise_lf)
- [visualise](#visualise)

#### List of available plugins (and their own commands)
- [grid](#grid)
  - [resolution](#resolution)
  - [space](#space)
  - [constrain](#constrain)  
- [mesh](#mesh)
  - [orient](#orient)
  - [fdownsample](#fdownample)
  - [symmetric](#symmetric)
  - [flip](#flip)

## Commands

#### BF
Path to BF.mat file

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: cellstr
matlabbatch{1}.spm.tools.beamforming.sources.BF = {'/path/to/BF/'};

% DAiSS-Wizard
% Default: REQUIRED
% Input Type: str or cellstr
S.BF = 'path/to/BF';
```

#### reduce_rank
Specify how many degrees of freedom we want in our forward solutions (for a given source) for both MEG and EEG simultaneously. Typically for most use cases we specify 2 for MEG (due to the silent radial compenent) and 3 for EEG. However there are some cases where this might be different (e.g. magnetic dipole modelling with MEG).
A two element vector is required, the first element is the MEG degrees of freedom and the second the EEG. (Both elements requied even if you are only using one modality).

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: int or double array (size: 1 x 2)
matlabbatch{1}.spm.tools.beamforming.sources.reduce_rank = [2 3];

% DAiSS-Wizard
% Default: [2 3]
% Input Type: int or double array (size: 1 x 2)
S.reduce_rank = [2 3];
```

#### keep_3d 
If there has been a rank reduction specified by [reduce_rank](#reduce_rank) then setting this to true ensures that there are still 3 reported lead field patterns for each source (rather than collapsed to to the number of degrees of freedom previously specified).

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: logical
matlabbatch{1}.spm.tools.beamforming.sources.keep_3d = true;

% DAiSS-Wizard
% Default: true
% Input Type: logical
S.keep_3d = true;
```

#### normalise_lf
Normalised each lead field such that the 2-norm of each is 1. Can be useful for some rudimentary depth correction (e.g. making beamformers less biased to deep sources, and the MNE-type solutions for superficial sources). 

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: false
matlabbatch{1}.spm.tools.beamforming.sources.normalise_lf = false;

% DAiSS-Wizard
% Default: false
% Input Type: logical
S.normalise_lf = false;
```

#### visualise
Would you like to see what the source space and sensors aligned looks like? (**WARNING**: Using the UK spelling here, bf_wizard_sources might check for this for an Americanization in the future, but using raw matlabbatch coding will probably trip over if spelled visualize!)

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: true
matlabbatch{1}.spm.tools.beamforming.sources.visualise = true;

% DAiSS-Wizard
% Default: true
% Input Type: logical
S.visualise = true;
```

## Plugins (source space generation mathods)

### grid

This pugin facilitates seting up a regular volumetric grid within a spefied region of the brain. Good for spatial filter methods (such as beamformers) but can be used in other cases if you want.

#### resolution
Select the resolution of the grid (in mm)

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: 5
matlabbatch{1}.spm.tools.beamforming.sources.plugin.grid.resolution = 5;

% DAiSS-Wizard
% Default: 5
% Input Type: numeric
S.method = 'grid'; 
S.grid.resolution = 5;
```

#### space
Select the space which the grid is contructed before moving to the invdividual.

Options
- MNI template _(Default, recommended)._ Grid constucted in MNI space and then warped back into the space defined in the [data](01_data.md) module. 
- MNI-Aligned
- Head
- Native

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: str
matlabbatch{1}.spm.tools.beamforming.sources.plugin.grid.space = 'MNI template';

% DAiSS-Wizard
% Default: 'MNI template'
% Input Type: str
S.method = 'grid'; 
S.grid.template = 'MNI template';
```

#### constrain
Select the boundary of which sources which exist outside of it (post-warp) are not included in source analysis.

Options
- iskull _(Default, recommended)._ Uses the inner skull compartment. 
- scalp. Uses the scalp mesh.

```matlab

% matlabbatch
% Default: 'iskull'
% Input Type: str
matlabbatch{1}.spm.tools.beamforming.sources.plugin.grid.constrain = 'iskull';

% DAiSS-Wizard
% Default: 'iskull'
% Input Type: str
S.method = 'grid'; 
S.grid.constrain = 'iskull';
```

### mesh

Generates the source space to be on the vertices of a mesh representing the cortex (or totally bonkers and custom like [a dragon perahps](https://github.com/tierneytim/OPM#c5)). Similar to how MEEG source inversion in the core of SPM operates. Good for distributed soure soltions).

#### orient

Decide how to constrain the orientaiton of the sources at each mesh vertex. 

Options:
- unoriented _(Default)._ No contraint on orientation, allows for a full vector analysis or a data driven orientation optimisation later. 
- original. Forces the orientation of the source to be normal to the mesh surface. 
- downsample. Behaves similarly to 'original' (useful for downsampled source spaces. see fdownsample)


```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: str
matlabbatch{1}.spm.tools.beamforming.sources.plugin.mesh.orient = 'unoriented';

% DAiSS-Wizard
% Default: 'unoriented'
% Input Type: str
S.method = 'mesh'; 
S.mesh.orient = 'unoriented';
```

#### fdownsample
Select the downsampling factor of the mesh to reduce the number of sources modelled (particulrly useful if the cortex mesh is a full resolution freesurfer reconstruction). a value of *n* would mean that every n<sup>th</sup> vertex is included in the source model.

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: numeric
matlabbatch{1}.spm.tools.beamforming.sources.plugin.mesh.fdownsample = 1;

% DAiSS-Wizard
% Default: 1
% Input Type: numeric
S.method = 'mesh'; 
S.mesh.fdownsample = 1;
```

#### symmetric
Do you want the source space to be summetric in the mid-saggital plane? If so, which hemisphere do you want mirrored?

Options:
- no *(default)*. Do not mirror.
- left. Mirror the left hemisphere.
- right. Mirror the right hemisphere.

```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: str
matlabbatch{1}.spm.tools.beamforming.sources.plugin.mesh.symmetric = 'no';

% DAiSS-Wizard
% Default: 'no'
% Input Type: str
S.method = 'mesh'; 
S.mesh.symmetric = 'no';
```

#### flip
Do you want to flip the source space around the mid-saggital plane?


```matlab

% matlabbatch
% Default: REQUIRED
% Input Type: logical or str
matlabbatch{1}.spm.tools.beamforming.sources.plugin.mesh.flip = false;

% DAiSS-Wizard
% Default: false
% Input Type: logical or str
S.method = 'mesh'; 
S.mesh.symmetric = false;
```
