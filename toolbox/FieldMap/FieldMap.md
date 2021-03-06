# FieldMap Toolbox - Description and Usage

This manual describes how to use the FieldMap toolbox for creating 
unwrapped field maps and unwarping geometrically distorted EPI 
images. The toolbox is designed to be interactive so that the user
can see the effect of applying different field maps and unwarping
parameters to EPI images. However, once a set of parameters has been
established for a specific scanning protocol, the routines used by 
the toolbox can also be scripted. This is described in
`FieldMap_ngui.m`.

Input parameters and the mode in which the toolbox works can be 
customised using the defaults file called pm_defaults.m.
For information on these other topics type the name (e.g. 
pm_defaults.m) in the "topic" text field of the spm_help window. 
Other useful help files (for typing in "topic" field of help window)
include:

* `FieldMap_principles.man`: for an explanation of the theory behind the field-map based
   unwarping.
* `pm_unwrap.m`: for specifics on the implementation of phase-unwrapping.
* `pm_ff_unwarp.m`: for specifics on watershed based unwrapping.
* `pm_initial_regions.m`, `pm_create_connectogram.m` and 
`pm_merge_regions.m`: for specifics of region-merging based unwrapping (default).

The toolbox is organised in two parts. The first part deals with the
calculation of an unwrapped field map given MRI data from a field map
acquisition. The second part converts the calculated field map into
a voxel displacement map (VDM) which can then be used to unwarp a 
selected EPI image.

This document gives a description of the Toolbox followed by a 
usage manual.

## DESCRIPTION OF THE TOOLBOX

### Part 1 - Loading and unwrapping field map data

#### Loading field map data

The field map acquisition can be based on an EPI or non-EPI 
sequence. The toolbox can take different kinds of input to 
allow for the different ways in which field map data may be
acquired. An unwrapped field map (in Hz) can be created from 
any of the following data situations:

1) Two pairs of real and imaginary image volumes, one for a 
shorter and one for a longer echo time (ie 4 images).

2) A phase and magnitude image (2 images). The phase image has
been created by the vendor sequence from two echo times
acquisitions (ie 2 images).

3) Two pairs of phase and magnitude images, one for each echo time.
(ie 4 images).
The input data format is specified by the default variable
`pm_def.INPUT_DATA_FORMAT` which works as follows:

a) `pm_def.INPUT_DATA_FORMAT = 'RI'` - the input must be two real 
and imaginary image pairs.

b) `pm_def.INPUT_DATA_FORMAT = 'PM'` - the input must be either one 
or two phase and magnitude image pairs. The units for the phase
images MUST BE RADIANS BETWEEN +pi and -pi. NB - this is commonly 
not radians and WILL therefore require scaling.

In both cases the input data must be in the Analyze format. 
Converting DICOM files to Analyze can be done using the DICOM
toolbox in SPM2. DICOM and other image formats can be converted to 
Analyze using MRIcro:
http://www.cla.sc.edu/psyc/faculty/rorden/mricro.html
Example field map data sets can be found in the examples directory 
of the distribution. 

The short and long echo times associated with the field map
acquisition must be specified in ms using the following parameters:
pm_def.SHORT_ECHO_TIME and pm_def.LONG_ECHO_TIME. Both of these
values are required even if a single phase and magnitude image is 
used as input.

#### Unwrapping

Once the measured field map data has been loaded, the unwrapped
field map can be calculated. This involves the following steps for
which certain parameters must be set:

a) The map of phase changes associated with the measured field map
is calculated from the input data.

b) The resulting phase map is unwrapped using the method specified
by pm_def.UNWRAPPING_METHOD = 'Mark3D' or 'Mark2D' or 'Huttonish'. 
For a description of these different methods see m_unwrap.m or
FieldMap_principles.man. The default option is 'Mark3D'.

c) A mask is created so that unwrapping only occurs in regions where 
there is signal. If necessary, this mask can be expanded so that any 
voxel that hasn't been unwrapped and is less than pm_def.PAD/2
voxels away from an unwrapped one will be replaced by an average of 
the surrounding unwrapped voxels. This can be done by setting the
parameter pm_def.PAD to a value greater than 0. The default
value is 0 but a value > 0 (eg 10) may be necessary if normal
smoothing is chosen instead of weighted smoothing (as explained in
the next step).

d) A weighted gaussian smoothing (weighted by the inverse of the 
noise) is performed on the unwrapped phase-map if the parameter 
pm_def.WS = 1. If pm_def.WS = 0, a normal smoothing is done.
The weighted smoothing is particularly slow on large data sets 
ie high resolution. If field maps are acquired at high resolution
then it is recommended to use pm_def.WS = 0 and do some padding
of the intensity mask eg pm_def.PAD = 10. The size of the Gaussian 
filter used to implement either weighted or normal smoothing of 
the unwrapped maps is usually set to pm_def.FWHM = 10.

#### Saving the field map

The calculated field map can be written out to disk as an Analyze 
image and will be saved with the filename: 
fpm_NAME-OF-FIRST-INPUT-IMAGE.img.
Also, at this stage behind the scenes, the field map is converted
to a VDM using the default parameters and saved with the filename
vdm_NAME-OF-FIRST-INPUT-IMAGE.img. This file is overwritten whenever 
the field map is recalculated or when any conversion parameters 
are changed (as described below). The VDM file can be used as the 
optional input field map for Realign & Unwarp in SPM2.

#### Loading a precalculated field map

Alternatively, it is also possible to load a precalculated unwrapped
field map (fpm_*.img). This should be a single image volume with
units of Hz and be in Analyze format. This precalculated field map 
could have been created and saved using the FieldMap toolbox or by 
other means. Examples of precalculated field maps can be found in
the examples directory.

Once calculated or loaded, the field map is displayed in a figure 
window and the field at different points can be explored.

### Part 2 - Convert the field map to a VDM and unwarp the EPI

The second part of the toolbox deals with converting the field map
(in Hz) to a voxel displacement map, VDM, in units of voxels
followed by EPI unwarping. This requires a few parameters to be set 
and works as follows:

#### Converting field map to VDM

a) The field map is multiplied by the total EPI readout time (in
ms) of the EPI image to be unwarped, resulting in a VDM.
This is specified by pm_def.TOTAL_EPI_READOUT_TIME 
(eg typically 10s of ms).

NB: The total EPI readout time is the time taken to acquire all the 
phase encode steps required to cover k-space (ie one image slice). For 
example, if the EPI sequence has 64 phase encode steps, the total
readout time is the time taken to acquire 64 echoes. 
     ie total readout time = number of echoes * echo spacing.
This time does not include i) the duration of of the excitation, 
ii) the delay between the excitation and the start of the acquisition 
or iii) time for fat saturation etc.

b) The VDM is multiplied by +/-1 to indicate whether the K-space 
traversal for the data acquisition has a +ve or -ve blip direction.
This will ensure that the unwarping is performed in the correct
direction and is specified by:
```
pm_def.K_SPACE_TRAVERSAL_BLIP_DIR = +/- 1.
```

c) The toolbox must know if the field map is based on an EPI or 
non-EPI acquisition. If using an EPI-based field map, the VDM must 
be inverted since the field map was acquired in warped space. 
This is specified by pm_def.EPI_BASED_FIELDMAPS = 1 or 0.

d) There is an option to apply Jacobian Modulation to the unwarped 
EPI image. This modulates the intensity of the unwarped
image so that in regions where voxels were compressed, the intensity
is decreased and where voxels were stretched, the intensities are
increased slightly. The modulation involves multiplying the unwarped
EPI by 1 + the 1-d derivative of the VDM in the phase direction.
An intensity adjustment of this nature may improve the
coregistration results between an unwarped EPI and an undistorted
image.

The resulting VDM is used to unwarp a selected EPI. The warped and 
the unwarped EPI are displayed in the figure window so that the
effects of the unwarping can be inspected.

When any of the above conversion parameters are changed or a new EPI
is selected, a new VDM is created and saved with the filename
vdm_NAME-OF-FIRST-INPUT-IMAGE.img.

Any previous copy of the .img file is overwritten, but the
corresponding .mat file is retained. It is done this way because 
the VDM may have already been coregistered to the EPI (as described
below). Then, for an EPI-based VDM, the match between the VDM and
the EPI will still be valid even if any of the above parameters have
been changed. If the VDM is non-EPI-based and any of the above
parameters are changed, the match between the VDM and the EPI may no
longer be valid. In this case a warning is given to the user that it
may be necessary to perform the coregistration again.

#### Coregistering the VDM to the EPI

When a VDM is first created, it is assumed that the EPI and VDM are
in the same space. This may not always be the case. Therefore the
toolbox allows a magnitude image associated with the acquired field
map to be coregistered to the selected EPI. This makes use of the
SPM2 mutual information coregistration routine.

#### Comparing unwarp results with a structural image

Finally, the toolbox allows a structural image to be loaded and
displayed so that the effects of unwarping can be compared with
an undistorted image. The structural image can also be 
coregistered to the unwarped EPI.

## FIELDMAP TOOLBOX USAGE

The following describes the specific buttons and options available
in the toolbox:

### Create field map in Hz

#### Short TE/Long TE (ms) 

Give shortest/longest echo-time in ms. This sets the parameters: 
pm_def.SHORT_ECHO_TIME and pm_def.LONG_ECHO_TIME. NB  - the default
values will remain as defined in pm_defaults.m. 

#### Load Real/Load Imag 
   
Load Analyze images containing:
a) real part of short echo-time image.
b) imaginary part of short echo-time image.
c) real part of long echo-time image.
d) imaginary part of long echo-time image.

**OR**

#### Load Phase/Magnitude -  
    
Load Analyze images containing:

a) phase of short echo-time image OR single phase image.

b) magnitude of short echo-time image OR single magnitude image.

c) phase of long echo-time image OR LEAVE EMPTY if input consists of
a single phase and magnitude image.

d) magnitude of long echo-time image OR LEAVE EMPTY if input
consists of a single phase and magnitude image.

#### Calculate

This creates an unwrapped field map in Hz which is stored in memory.
The field map is also converted to a voxel displacement map (VDM) 
using the default paramaters and is written out with the filename
vdm_NAME-OF-FIRST-INPUT-IMAGE.img overwriting any previous copy of 
the .img file but not the .mat file. The VDM file can be used as the
optional input field map for Realign * Unwarping in SPM2.
The finished field map is displayed in the figure window.

#### Write 

Write out the finished field map (in Hz) as an analyze image with
the filename fpm_NAME-OF-FIRST-INPUT-IMAGE.img. The field map is 
written out with data type = double so that the unwarping results
will be identical whether the field map is created and used 
immediately or written out and used at a later stage.

#### Load (precalculated field map)  
    
This prompts the user to select an Analyze image containing
the finished field map (in Hz). This is displayed in the figure
window. 

#### Field map value (Hz) 

This returns the value of the field map in Hz at the location 
specified by the mouse pointer in the figure window.

### Create voxel displacement map (VDM) and unwarp EPI


When any of the parameters below are changed, a new VDM is created
and written out as vdm_NAME-OF-FIRST-INPUT-IMAGE.img. The
vdm_NAME-OF-FIRST-INPUT-IMAGE.mat file is not updated unless 
'match VDM to EPI' is selected as described below.

#### EPI-based field map - Yes/No

If the field map is based on EPI data, it must be inverted before
being used to unwarp the EPI. If the field map is based on non-EPI 
(eg Flash) data it can be used directly after being sampled in the 
space and resolution of the EPI image.

#### Polarity of phase-encode blips - +ve/-ve

When images are acquired K-space can be traversed using positive
or negative phase-encode blips. This direction will influence the 
geometric distortions in terms of whether the affected regions of 
the image are stretched or compressed. This effect can be seen by 
toggling this value in the toolbox using one of the example data 
sets.

#### Apply Jacobian modulation - Yes/No

Do Jacobian Modulation to adjust the intensities of voxels that 
have been stretched or compressed.

#### Total EPI readout time (ms)

Give the total time in ms for the readout of the EPI echo train.

#### Load EPI image

This prompts the user to select a sample modulus EPI Analyze image.
This image is automatically unwarped using the VDM calculated with
the current parameters. The warped and the unwarped image are
displayed in the orthviews figure window underneath the field map.

#### Match VDM to EPI 

This option allows the field map to be matched to the EPI image
before it is used to unwarp the EPI. It requires a magnitude image
which is in the same space as the field map. This can be achieved 
in one of 3 different ways:

a) If the field map was just created using a real and imaginary
image, then the magnitude image is automatically created using the
square root of the sum of the squares of the real and imaginary
short echo time images. (NB - the choice of short as opposed to long
echo time images is arbitrary). The magnitude image is automatically
written out with the filename mag_NAME-OF-FIRST-INPUT-IMAGE.img. 

b) If a phase and magnitude image are used to create the field map, 
the magnitude short echo time image is used (NB - the choice of 
short as opposed to long echo time image is arbitrary).

c) If a precalculated field map has been loaded then the user is
prompted to select a magnitude image. 

The magnitude image is then coregistered to the EPI image. This
is done slightly differently depending on whether the field map
is based on EPI or non-EPI data.

#### If using an EPI field map:

a) The magnitude image is coregistered to the EPI.

b) The resulting transformation matrix is written to the .mat file
of the VDM and used to sample it in the space of the EPI before
unwarping the EPI and displaying the results in the figure window.

#### If using a non-EPI field map:

a) The magnitude image is sampled in the space of the EPI to be
unwarped.

b) The VDM is inverted so that it can be used to forward warp the
resampled magnitude image.

c) The VDM is resampled in the space of the resampled magnitude
image.

d) The VDM is used to forward warp the magnitude image resulting
in the magnitude image being similarly warped as the EPI. It is 
automatically written with the filename:
wfmag_NAME-OF-FIRST-INPUT-IMAGE.img. 

e) The wfmag_NAME-OF-FIRST-INPUT-IMAGE.img is coregistered to the 
EPI and its .mat file updated with the resulting transformation
matrix.

f) The resulting transformation matrix is written to the .mat file
of the VDM and used to sample it in the space of the EPI before
unwarping the EPI and displaying the results in the figure window. 

#### Write unwarped

Write unwarped EPI Analyze image with the filename uNAME_OF_EPI.img.

#### Load structural

This prompts the user to select a structural image which is
displayed in the figure window below the other images.

#### MatchStructural

This option coregisters the structural image to the unwarped EPI
and writes the resulting transformation matrix to the .mat file of
the selected structural image.

#### Help

This calls spm_help FieldMap.man to display this help file.

#### Quit

This quits the toolbox and closes all windows associated with it.

```
Jesper Andersson & Chloe Hutton
```
