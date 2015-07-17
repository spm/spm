function fieldmap = tbx_cfg_fieldmap
% MATLABBATCH Configuration file for toolbox 'FieldMap'
%__________________________________________________________________________
% Copyright (C) 2008-2015 Wellcome Trust Centre for Neuroimaging

% $Id: tbx_cfg_fieldmap.m 6501 2015-07-17 14:32:09Z spm $


addpath(fullfile(spm('dir'),'toolbox','FieldMap'));

%==========================================================================
% Default values that are common to all fieldmap jobs
%==========================================================================

%--------------------------------------------------------------------------
% et Echo times [short TE long TE]
%--------------------------------------------------------------------------
et         = cfg_entry;
et.tag     = 'et';
et.name    = 'Echo times [short TE long TE]';
et.help    = {'Enter the short and long echo times (in ms) of the data used to acquire the field map.'};
et.strtype = 'e';
et.num     = [1  2];

%--------------------------------------------------------------------------
% maskbrain Mask brain
%--------------------------------------------------------------------------
maskbrain        = cfg_menu;
maskbrain.tag    = 'maskbrain';
maskbrain.name   = 'Mask brain';
maskbrain.help   = {
    'Select masking or no masking of the brain. If masking is selected,'
    'the magnitude image is used to generate a mask of the brain.'
}';
maskbrain.labels = {
                    'Mask brain'
                    'No brain masking'
}';
maskbrain.values = {1 0};

%--------------------------------------------------------------------------
% blipdir Blip direction
%--------------------------------------------------------------------------
blipdir        = cfg_menu;
blipdir.tag    = 'blipdir';
blipdir.name   = 'Blip direction';
blipdir.help   = {'Enter the blip direction. This is the polarity of the phase-encode blips describing the direction in which k-space is traversed along the y-axis during EPI acquisition with respect to the coordinate system used in SPM. In this coordinate system, the phase encode direction corresponds with the y-direction and is defined as positive from the posterior to the anterior of the head.'};
blipdir.labels = {'-1' '1'};
blipdir.values = {-1 1};

%--------------------------------------------------------------------------
% tert Total EPI readout time
%--------------------------------------------------------------------------
tert         = cfg_entry;
tert.tag     = 'tert';
tert.name    = 'Total EPI readout time';
tert.help    = {
                'Enter the total EPI readout time (in ms). This is the time taken to '
                'acquire all of the phase encode steps required to cover k-space (ie one image slice). '
                'For example, if the EPI sequence has 64 phase encode steps, the total readout time is '
                'the time taken to acquire 64 echoes, e.g. '
                'total readout time = number of echoes x echo spacing. '
                'This time does not include i) the duration of the excitation, ii) the delay between, '
                'the excitation and the start of the acquisition or iii) time for fat saturation etc.'
}';
tert.strtype = 'e';
tert.num     = [1  1];

%--------------------------------------------------------------------------
% epifm EPI-based field map?
%--------------------------------------------------------------------------
epifm        = cfg_menu;
epifm.tag    = 'epifm';
epifm.name   = 'EPI-based field map?';
epifm.help   = {'Select non-EPI or EPI based field map. The field map data may be acquired using a non-EPI sequence (typically a gradient echo sequence) or an EPI sequence. The processing will be slightly different for the two cases. If using an EPI-based field map, the resulting Voxel Displacement Map will be inverted since the field map was acquired in distorted space.'};
epifm.labels = {'non-EPI' 'EPI'};
epifm.values = {0 1};

%--------------------------------------------------------------------------
% ajm Jacobian modulation?
%--------------------------------------------------------------------------
ajm        = cfg_menu;
ajm.tag    = 'ajm';
ajm.name   = 'Jacobian modulation?';
ajm.help   = {'Select whether or not to use Jacobian modulation. This will adjust the intensities of voxels that have been stretched or compressed but in general is not recommended for EPI distortion correction'};
ajm.labels = {
              'Do not use'
              'Use'
}';
ajm.values = {0 1};
ajm.def    = @(val)pm_get_defaults('DO_JACOBIAN_MODULATION', val{:});

%--------------------------------------------------------------------------
% method Unwrapping method
%--------------------------------------------------------------------------
method        = cfg_menu;
method.tag    = 'method';
method.name   = 'Unwrapping method';
method.help   = {'Select method for phase unwrapping'};
method.labels = {
                 'Mark3D'
                 'Mark2D'
                 'Huttonish'
}';
method.values = {
                 'Mark3D'
                 'Mark2D'
                 'Huttonish'
}';
method.def    = @(val)pm_get_defaults('UNWRAPPING_METHOD', val{:});

%--------------------------------------------------------------------------
% fwhm FWHM
%--------------------------------------------------------------------------
fwhm         = cfg_entry;
fwhm.tag     = 'fwhm';
fwhm.name    = 'FWHM';
fwhm.help    = {'FWHM of Gaussian filter used to implement weighted smoothing of unwrapped maps.'};
fwhm.strtype = 'e';
fwhm.num     = [1  1];
fwhm.def     = @(val)pm_get_defaults('FWHM', val{:});

%--------------------------------------------------------------------------
% pad pad
%--------------------------------------------------------------------------
pad         = cfg_entry;
pad.tag     = 'pad';
pad.name    = 'pad';
pad.help    = {'Size of padding kernel if required.'};
pad.strtype = 'e';
pad.num     = [1  1];
pad.def     = @(val)pm_get_defaults('PAD', val{:});

%--------------------------------------------------------------------------
% ws Weighted smoothing
%--------------------------------------------------------------------------
ws         = cfg_menu;
ws.tag     = 'ws';
ws.name    = 'Weighted smoothing';
ws.help    = {'Select normal or weighted smoothing.'};
ws.labels = {
             'Weighted Smoothing'
             'No weighted smoothing'
}';
ws.values{1} = 1;
ws.values{2} = 0;
ws.def     = @(val)pm_get_defaults('WS', val{:});

%--------------------------------------------------------------------------
% uflags uflags
%--------------------------------------------------------------------------
uflags         = cfg_branch;
uflags.tag     = 'uflags';
uflags.name    = 'uflags';
uflags.val     = {method fwhm pad ws };
uflags.help    = {'Different options for phase unwrapping and field map processing'};

%--------------------------------------------------------------------------
% template Template image for brain masking
%--------------------------------------------------------------------------
template         = cfg_files;
template.tag     = 'template';
template.name    = 'Template image for brain masking';
template.help    = {'Select template file for segmentation to create brain mask'};
template.filter = 'nii';
template.ufilter = '.*';
template.num     = [1 1];
template.def     = @(val)pm_get_defaults('MFLAGS.TEMPLATE', val{:});

%--------------------------------------------------------------------------
% fwhm FWHM
%--------------------------------------------------------------------------
fwhm         = cfg_entry;
fwhm.tag     = 'fwhm';
fwhm.name    = 'FWHM';
fwhm.help    = {'FWHM of Gaussian filter for smoothing brain mask.'};
fwhm.strtype = 'e';
fwhm.num     = [1  1];
fwhm.def     = @(val)pm_get_defaults('MFLAGS.FWHM', val{:});

%--------------------------------------------------------------------------
% nerode Number of erosions
%--------------------------------------------------------------------------
nerode         = cfg_entry;
nerode.tag     = 'nerode';
nerode.name    = 'Number of erosions';
nerode.help    = {'Number of erosions used to create brain mask.'};
nerode.strtype = 'e';
nerode.num     = [1  1];
nerode.def     = @(val)pm_get_defaults('MFLAGS.NERODE', val{:});

%--------------------------------------------------------------------------
% ndilate Number of dilations
%--------------------------------------------------------------------------
ndilate         = cfg_entry;
ndilate.tag     = 'ndilate';
ndilate.name    = 'Number of dilations';
ndilate.help    = {'Number of dilations used to create brain mask.'};
ndilate.strtype = 'e';
ndilate.num     = [1  1];
ndilate.def     = @(val)pm_get_defaults('MFLAGS.NDILATE', val{:});

%--------------------------------------------------------------------------
% thresh Threshold
%--------------------------------------------------------------------------
thresh         = cfg_entry;
thresh.tag     = 'thresh';
thresh.name    = 'Threshold';
thresh.help    = {'Threshold used to create brain mask from segmented data.'};
thresh.strtype = 'e';
thresh.num     = [1  1];
thresh.def     = @(val)pm_get_defaults('MFLAGS.THRESH', val{:});

%--------------------------------------------------------------------------
% reg Regularization
%--------------------------------------------------------------------------
reg         = cfg_entry;
reg.tag     = 'reg';
reg.name    = 'Regularization';
reg.help    = {'Regularization value used in the segmentation. A larger value helps the segmentation to converge.'};
reg.strtype = 'e';
reg.num     = [1  1];
reg.def     = @(val)pm_get_defaults('MFLAGS.REG', val{:});

%--------------------------------------------------------------------------
% mflags mflags
%--------------------------------------------------------------------------
mflags         = cfg_branch;
mflags.tag     = 'mflags';
mflags.name    = 'mflags';
mflags.val     = {template fwhm nerode ndilate thresh reg };
mflags.help    = {'Different options used for the segmentation and creation of the brain mask.'};

%--------------------------------------------------------------------------
% defaultsval Defaults values
%--------------------------------------------------------------------------
defaultsval         = cfg_branch;
defaultsval.tag     = 'defaultsval';
defaultsval.name    = 'Defaults values';
defaultsval.val     = {et maskbrain blipdir tert epifm ajm uflags mflags };
defaultsval.help    = {'Defaults values'};

%--------------------------------------------------------------------------
% defaultsfile Defaults File
%--------------------------------------------------------------------------
defaultsfile         = cfg_files;
defaultsfile.tag     = 'defaultsfile';
defaultsfile.name    = 'Defaults File';
defaultsfile.help    = {'Select the ''pm_defaults*.m'' file containing the parameters for the field map data. Please make sure that the parameters defined in the defaults file are correct for your field map and EPI sequence. To create your own customised defaults file, either edit the distributed version and/or save it with the name ''pm_defaults_yourname.m''.'};
defaultsfile.filter  = 'm';
defaultsfile.dir     = fileparts(mfilename('fullpath'));
defaultsfile.ufilter = '^pm_defaults.*\.m$';
defaultsfile.num     = [1 1];
defaultsfile.def     = @(val)pm_get_defaults('defaultsfilename', val{:});

%--------------------------------------------------------------------------
% defaults FieldMap defaults
%--------------------------------------------------------------------------
defaults         = cfg_choice;
defaults.tag     = 'defaults';
defaults.name    = 'FieldMap defaults';
defaults.help    = {'FieldMap default values can be entered as a file or set of values.'};
defaults.values  = {defaultsval defaultsfile};

%--------------------------------------------------------------------------
% epi Select EPI to Unwarp
%--------------------------------------------------------------------------
epi         = cfg_files;
epi.tag     = 'epi';
epi.name    = 'Select EPI to Unwarp';
epi.help    = {'Select a single image to distortion correct. The corrected image will be saved with the prefix u. Note that this option is mainly for quality control of correction so that the original and distortion corrected images can be displayed for comparison. To unwarp multiple images please use either Realign & Unwarp or Apply VDM.'};
epi.filter = 'image';
epi.ufilter = '.*';
epi.num     = [1 1];

%--------------------------------------------------------------------------
% session Session
%--------------------------------------------------------------------------
session         = cfg_branch;
session.tag     = 'session';
session.name    = 'Session';
session.val     = {epi};
session.help    = {'Data for this session.'};

%--------------------------------------------------------------------------
% generic EPI Sessions
%--------------------------------------------------------------------------
generic1         = cfg_repeat;
generic1.tag     = 'generic1';
generic1.name    = 'EPI Sessions';
generic1.help    = {'If a single set of field map data will be used for multiple EPI runs/sessions, select the first EPI in each run/session. A VDM file will created for each run/session, matched to the first EPI in each run/session and saved with a unique name extension.'};
generic1.values  = {session};
generic1.val     = {session};
generic1.num     = [1 Inf];

%--------------------------------------------------------------------------
% matchvdm Match VDM to EPI?
%--------------------------------------------------------------------------
matchvdm        = cfg_menu;
matchvdm.tag    = 'matchvdm';
matchvdm.name   = 'Match VDM to EPI?';
matchvdm.help   = {'Match VDM file to EPI image. This will coregister the field map data to the selected EPI for each run/session.'};
matchvdm.labels = {
                   'match VDM'
                   'none'
}';
matchvdm.values = {1 0};

%--------------------------------------------------------------------------
% sessname Name extension for session specific vdm files
%--------------------------------------------------------------------------
sessname         = cfg_entry;
sessname.tag     = 'sessname';
sessname.name    = 'VDM filename extension';
sessname.help    = {'Filename extension for run/session specific VDM files. The extension will be followed by an incremented integer for run/session number.'};
sessname.strtype = 's';
sessname.num     = [1  Inf];
sessname.def     = @(val)pm_get_defaults('sessname', val{:});

%--------------------------------------------------------------------------
% writeunwarped Write unwarped EPI?
%--------------------------------------------------------------------------
writeunwarped        = cfg_menu;
writeunwarped.tag    = 'writeunwarped';
writeunwarped.name   = 'Write unwarped EPI?';
writeunwarped.help   = {'Write out distortion corrected EPI image. The image is saved with the prefix u. Note that this option is mainly for quality control of correction so that the original and distortion corrected images can be displayed for comparison. To unwarp multiple images please use either Realign & Unwarp or Apply VDM.'};
writeunwarped.labels = {
                        'write unwarped EPI'
                        'none'
}';
writeunwarped.values = {1 0};

%--------------------------------------------------------------------------
% anat Select anatomical image for comparison
%--------------------------------------------------------------------------
anat         = cfg_files;
anat.tag     = 'anat';
anat.name    = 'Anatomical image for comparison';
anat.help    = {'Select an anatomical image for comparison with the distortion corrected EPI or leave empty. Note that this option is mainly for quality control of correction.'};
anat.filter  = 'image';
anat.ufilter = '.*';
anat.num     = [0 1];
anat.val     = {''};

%--------------------------------------------------------------------------
% matchanat Match anatomical image to EPI?
%--------------------------------------------------------------------------
matchanat        = cfg_menu;
matchanat.tag    = 'matchanat';
matchanat.name   = 'Match anatomical image to EPI?';
matchanat.help   = {'Match the anatomical image to the distortion corrected EPI. Note that this option is mainly for quality control of correction allowing for visual inspection and comparison of the distortion corrected EPI.'};
matchanat.labels = {
                    'none'
                    'match anat'
}';
matchanat.values = {0 1};

%==========================================================================
% Calculate vdm* file
%==========================================================================

%--------------------------------------------------------------------------
% phase Phase Image
%--------------------------------------------------------------------------
phase         = cfg_files;
phase.tag     = 'phase';
phase.name    = 'Phase Image';
phase.help    = {'Select a single phase image. This should be the result from the subtraction of two phase images (where the subtraction is usually done automatically by the scanner software). The phase image will be scaled between +/- PI.'};
phase.filter  = 'image';
phase.ufilter = '.*';
phase.num     = [1 1];

%--------------------------------------------------------------------------
% magnitude Magnitude Image
%--------------------------------------------------------------------------
magnitude         = cfg_files;
magnitude.tag     = 'magnitude';
magnitude.name    = 'Magnitude Image';
magnitude.help    = {'Select a single magnitude image. This is used for masking the phase information and coregistration with the EPI data. If two magnitude images are available, select the one acquired at the shorter echo time because it will have greater signal'};
magnitude.filter  = 'image';
magnitude.ufilter = '.*';
magnitude.num     = [1 1];

%--------------------------------------------------------------------------
% presubphasemag
%--------------------------------------------------------------------------
presubphasemag      = cfg_branch;
presubphasemag.tag  = 'presubphasemag';
presubphasemag.name = 'Presubtracted Phase and Magnitude Data';
presubphasemag.val  = {phase magnitude};
presubphasemag.help = {'Calculate a voxel displacement map (VDM) from presubtracted phase and magnitude field map data. This option expects a single magnitude image and a single phase image resulting from the subtraction of two phase images (where the subtraction is usually done automatically by the scanner software). The phase image will be scaled between +/- PI.'};


%--------------------------------------------------------------------------
% shortreal Short Echo Real Image
%--------------------------------------------------------------------------
shortreal         = cfg_files;
shortreal.tag     = 'shortreal';
shortreal.name    = 'Short Echo Real Image';
shortreal.help    = {'Select short echo real image.'};
shortreal.filter = 'image';
shortreal.ufilter = '.*';
shortreal.num     = [1 1];

%--------------------------------------------------------------------------
% shortimag Short Echo Imaginary Image
%--------------------------------------------------------------------------
shortimag         = cfg_files;
shortimag.tag     = 'shortimag';
shortimag.name    = 'Short Echo Imaginary Image';
shortimag.help    = {'Select short echo imaginary image.'};
shortimag.filter = 'image';
shortimag.ufilter = '.*';
shortimag.num     = [1 1];

%--------------------------------------------------------------------------
% longreal Long Echo Real Image
%--------------------------------------------------------------------------
longreal         = cfg_files;
longreal.tag     = 'longreal';
longreal.name    = 'Long Echo Real Image';
longreal.help    = {'Select long echo real image.'};
longreal.filter = 'image';
longreal.ufilter = '.*';
longreal.num     = [1 1];

%--------------------------------------------------------------------------
% longimag Long Echo Imaginary Image
%--------------------------------------------------------------------------
longimag         = cfg_files;
longimag.tag     = 'longimag';
longimag.name    = 'Long Echo Imaginary Image';
longimag.help    = {'Select long echo imaginary image.'};
longimag.filter = 'image';
longimag.ufilter = '.*';
longimag.num     = [1 1];

%--------------------------------------------------------------------------
% realimag
%--------------------------------------------------------------------------
realimag      = cfg_branch;
realimag.tag  = 'realimag';
realimag.name = 'Real and Imaginary Data';
realimag.val  = {shortreal shortimag longreal longimag};
realimag.help = {'Calculate a voxel displacement map (VDM) from real and imaginary field map data. This option expects two real and imaginary pairs of data of two different echo times. The phase images will be scaled between +/- PI.'};

%--------------------------------------------------------------------------
% shortphase Short Echo Phase Image
%--------------------------------------------------------------------------
shortphase         = cfg_files;
shortphase.tag     = 'shortphase';
shortphase.name    = 'Short Echo Phase Image';
shortphase.help    = {'Select short echo phase image.'};
shortphase.filter  = 'image';
shortphase.ufilter = '.*';
shortphase.num     = [1 1];

%--------------------------------------------------------------------------
% shortmag Short Echo Magnitude Image
%--------------------------------------------------------------------------
shortmag         = cfg_files;
shortmag.tag     = 'shortmag';
shortmag.name    = 'Short Echo Magnitude Image';
shortmag.help    = {'Select short echo magnitude image.'};
shortmag.filter  = 'image';
shortmag.ufilter = '.*';
shortmag.num     = [1 1];

%--------------------------------------------------------------------------
% longphase Long Echo Phase Image
%--------------------------------------------------------------------------
longphase         = cfg_files;
longphase.tag     = 'longphase';
longphase.name    = 'Long Echo Phase Image';
longphase.help    = {'Select long echo phase image.'};
longphase.filter  = 'image';
longphase.ufilter = '.*';
longphase.num     = [1 1];

%--------------------------------------------------------------------------
% longmag Long Echo Magnitude Image
%--------------------------------------------------------------------------
longmag         = cfg_files;
longmag.tag     = 'longmag';
longmag.name    = 'Long Echo Magnitude Image';
longmag.help    = {'Select long echo magnitude image.'};
longmag.filter  = 'image';
longmag.ufilter = '.*';
longmag.num     = [1 1];

%--------------------------------------------------------------------------
% phasemag
%--------------------------------------------------------------------------
phasemag      = cfg_branch;
phasemag.tag  = 'phasemag';
phasemag.name = 'Phase and Magnitude Data';
phasemag.val  = {shortphase shortmag longphase longmag};
phasemag.help = {'Calculate a voxel displacement map (VDM) from double phase and magnitude field map data. This option expects two phase and magnitude pairs of data of two different echo times.'};

%--------------------------------------------------------------------------
% precalcfieldmap Precalculated field map
%--------------------------------------------------------------------------
precalcfieldmap1         = cfg_files;
precalcfieldmap1.tag     = 'precalcfieldmap';
precalcfieldmap1.name    = 'Precalculated field map';
precalcfieldmap1.help    = {'Select a precalculated field map. This should be a processed field map (ie phase unwrapped, masked if necessary and scaled to Hz), for example as generated by the FieldMap toolbox and are stored with fpm_* prefix.'};
precalcfieldmap1.filter  = 'image';
precalcfieldmap1.ufilter = '.*';
precalcfieldmap1.num     = [1 1];

%--------------------------------------------------------------------------
% magfieldmap Magnitude image in same space as field map
%--------------------------------------------------------------------------
magfieldmap         = cfg_files;
magfieldmap.tag     = 'magfieldmap';
magfieldmap.name    = 'Magnitude image in same space as field map';
magfieldmap.help    = {'Select magnitude image which is in the same space as the field map to do matching to EPI.'};
magfieldmap.filter  = 'image';
magfieldmap.ufilter = '.*';
magfieldmap.num     = [0 1];

%--------------------------------------------------------------------------
% precalcfieldmap
%--------------------------------------------------------------------------
precalcfieldmap      = cfg_branch;
precalcfieldmap.tag  = 'precalcfieldmap';
precalcfieldmap.name = 'Precalculated FieldMap (in Hz)';
precalcfieldmap.val  = {precalcfieldmap1 magfieldmap};
precalcfieldmap.help = {'Calculate a voxel displacement map (VDM) from a precalculated field map. This option expects a processed field map (ie phase unwrapped, masked if necessary and scaled to Hz). Precalculated field maps can be generated by the FieldMap toolbox and stored as fpm_* files.'};

%--------------------------------------------------------------------------
% data
%--------------------------------------------------------------------------
data         = cfg_choice;
data.tag     = 'data';
data.name    = 'Field map';
data.val     = {};
data.help    = {''};
data.values  = {presubphasemag realimag phasemag precalcfieldmap};

%--------------------------------------------------------------------------
% subj Subject
%--------------------------------------------------------------------------
subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {data defaults generic1 matchvdm sessname writeunwarped anat matchanat};
subj.help    = {'Data for this subject or field map session.'};

%--------------------------------------------------------------------------
% generic Data
%--------------------------------------------------------------------------
generic        = cfg_repeat;
generic.tag    = 'generic';
generic.name   = 'Data';
generic.help   = {'Subjects or sessions for which individual field map data has been acquired.'};
generic.values = {subj};
generic.val    = {subj};
generic.num    = [1 Inf];

%--------------------------------------------------------------------------
% calculatevdm calculate vdm* file 
%--------------------------------------------------------------------------
calculatevdm      = cfg_exbranch;
calculatevdm.tag  = 'calculatevdm';
calculatevdm.name = 'Calculate VDM';
calculatevdm.val  = {generic};
calculatevdm.help = {
    'Generate unwrapped field maps which are converted to voxel displacement maps (VDM) that can be used to unwarp geometrically distorted EPI images.'
    'The resulting VDM files are saved with the prefix vdm and can be applied to images using Apply VDM or in combination with Realign & Unwarp to calculate and correct for the combined effects of static and movement-related susceptibility induced distortions.'
    };
calculatevdm.prog = @FieldMap_calculatevdm;
calculatevdm.vout = @vout_calculatevdm;

%==========================================================================
% Apply vdm* file
%==========================================================================

%--------------------------------------------------------------------------
% scans Images
%--------------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Images';
scans.help    = {'Select scans for this session. These are assumed to be realigned to the first in the time series (e.g. using Realign: Estimate) but do not need to be resliced'};
scans.filter  = 'image';
scans.ufilter = '.*';
scans.num     = [1 Inf];

%--------------------------------------------------------------------------
% vdmfile selected vdm files
%--------------------------------------------------------------------------
vdmfile         = cfg_files;
vdmfile.tag     = 'vdmfile';
vdmfile.name    = 'Fieldmap (vdm* file)';
vdmfile.help    = {'Select VDM (voxel displacement map) for this session (e.g. created via FieldMap toolbox). This is assumed to be in alignment with the images selected for resampling (note this can be achieved via the FieldMap toolbox).'};
vdmfile.filter  = 'image';
vdmfile.ufilter = '.*';
vdmfile.num     = [1 1];

%--------------------------------------------------------------------------
% vdmapply session Session
%--------------------------------------------------------------------------
data         = cfg_branch;
data.tag     = 'data';
data.name    = 'Session';
data.val     = {scans vdmfile};
data.help    = {'Data for this session.'};

%--------------------------------------------------------------------------
% generic Data
%--------------------------------------------------------------------------
generic2         = cfg_repeat;
generic2.tag     = 'generic2';
generic2.name    = 'Data';
generic2.help    = {'Subjects or sessions for which VDM file is being applied to images.'};
generic2.values  = {data};
generic2.val     = {data};
generic2.num     = [1 Inf];

%--------------------------------------------------------------------------
% Reslice options for applyvdm
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% wrap Wrapping
%--------------------------------------------------------------------------
pedir         = cfg_menu;
pedir.tag     = 'pedir';
pedir.name    = 'Distortion direction';
pedir.help    = {'In which direction are the distortions? Any single dimension can be corrected therefore input data may have been acquired with phase encode directions in Y (most typical), X or Z'};
pedir.labels  = {
               'Posterior-Anterior (Y)'
               'Right-Left (X)'
               'Foot-Head (Z)'}';
pedir.values  = {2 1 3};
pedir.def     = @(val)pm_get_defaults('pedir', val{:});

%--------------------------------------------------------------------------
%  Which images to reslice?
%--------------------------------------------------------------------------
applyvdmwhich        = cfg_menu;
applyvdmwhich.tag    = 'which';
applyvdmwhich.name   = 'Resliced images';
applyvdmwhich.help   = {
                   'All Images (1..n) '
                   '  This applies the VDM and reslices all the images. '
                   'All Images + Mean Image '
                   '   This applies the VDM reslices all the images and creates a mean of the resliced images.'
}';
applyvdmwhich.labels = {
                  ' All Images (1..n)'
                  ' All Images + Mean Image'
}';
applyvdmwhich.values = {[2 0] [2 1]};
applyvdmwhich.def    = @(val)spm_get_defaults('realign.write.which', val{:});

%--------------------------------------------------------------------------
% rinterp Interpolation
%--------------------------------------------------------------------------
rinterp         = cfg_menu;
rinterp.tag     = 'rinterp';
rinterp.name    = 'Interpolation';
rinterp.help    = {'The method by which the images are sampled when being written in a different space. Nearest Neighbour is fastest, but not recommended for image realignment. Trilinear Interpolation is probably OK for PET, but not so suitable for fMRI because higher degree interpolation generally gives better results/* \cite{thevenaz00a,unser93a,unser93b}*/. Although higher degree methods provide better interpolation, but they are slower because they use more neighbouring voxels.'};
rinterp.labels  = {
                  'Nearest neighbour'
                  'Trilinear'
                  '2nd Degree B-spline '
                  '3rd Degree B-Spline'
                  '4th Degree B-Spline'
                  '5th Degree B-Spline '
                  '6th Degree B-Spline'
                  '7th Degree B-Spline'
}';
rinterp.values  = {0 1 2 3 4 5 6 7};
rinterp.def     = @(val)spm_get_defaults('realign.write.interp', val{:});

%--------------------------------------------------------------------------
% wrap Wrapping
%--------------------------------------------------------------------------
wrap        = cfg_menu;
wrap.tag    = 'wrap';
wrap.name   = 'Wrapping';
wrap.help   = {
               'This indicates which directions in the volumes the values should wrap around in.  For example, in MRI scans, the images wrap around in the phase encode direction, so (e.g.) the subject''s nose may poke into the back of the subject''s head. These are typically:'
               '    No wrapping - for PET or images that have already been spatially transformed. Also the recommended option if you are not really sure.'
               '    Wrap in  Y  - for (un-resliced) MRI where phase encoding is in the Y direction (voxel space) etc.'
}';
wrap.labels = {
               'No wrap'
               'Wrap X'
               'Wrap Y'
               'Wrap X & Y'
               'Wrap Z '
               'Wrap X & Z'
               'Wrap Y & Z'
               'Wrap X, Y & Z'
}';
wrap.values = {[0 0 0] [1 0 0] [0 1 0] [1 1 0] [0 0 1] [1 0 1] [0 1 1]...
               [1 1 1]};
wrap.def    = @(val)spm_get_defaults('realign.write.wrap', val{:});

%--------------------------------------------------------------------------
% mask Masking
%--------------------------------------------------------------------------
mask         = cfg_menu;
mask.tag     = 'mask';
mask.name    = 'Masking';
mask.help    = {'Because of subject motion, different images are likely to have different patterns of zeros from where it was not possible to sample data. With masking enabled, the program searches through the whole time series looking for voxels which need to be sampled from outside the original images. Where this occurs, that voxel is set to zero for the whole set of images (unless the image format can represent NaN, in which case NaNs are used where possible).'};
mask.labels  = {
               'Mask images'
               'Dont mask images'
}';
mask.values  = {1 0};
mask.def     = @(val)spm_get_defaults('realign.write.mask', val{:});

%--------------------------------------------------------------------------
% prefix Filename Prefix
%--------------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the distortion corrected image file(s). Default prefix is ''u''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.def     = @(val)spm_get_defaults('unwarp.write.prefix', val{:});

%--------------------------------------------------------------------------
% applyvdmroptions Reslicing Options
%--------------------------------------------------------------------------
applyvdmroptions      = cfg_branch;
applyvdmroptions.tag  = 'roptions';
applyvdmroptions.name = 'Reslice Options';
applyvdmroptions.val  = {pedir applyvdmwhich rinterp wrap mask prefix};
applyvdmroptions.help = {'Apply VDM reslice options'};

%--------------------------------------------------------------------------
% applyvdm Apply vdm* file to EPI files
%--------------------------------------------------------------------------
applyvdm      = cfg_exbranch;
applyvdm.tag  = 'applyvdm';
applyvdm.name = 'Apply VDM ';
applyvdm.val  = {generic2 applyvdmroptions};
applyvdm.help = {
    'Apply VDM (voxel displacement map) to resample voxel values in selected image(s). This allows a VDM to be applied to any images which are assumed to be already realigned (e.g. including EPI fMRI time series and DTI data).'
    'The VDM can be created from a field map acquisition using the FieldMap toolbox and comprises voxel shift values which describe geometric distortions occuring as a result of magnetic susceptbility artefacts. Distortions along any single dimension can be corrected therefore input data may have been acquired with phase encode directions in X, Y (most typical) and Z.'
    'The selected images are assumed to be realigned to the first in the time series (e.g. using Realign: Estimate) but do not need to be resliced. The VDM is assumed to be in alignment with the images selected for resampling (note this can be achieved via the FieldMap toolbox). The resampled images are written to the input subdirectory with the same (prefixed) filename.'
    'e.g. The typical processing steps for fMRI time series would be 1) Realign: Estimate, 2) FieldMap to create VDM, 3) Apply VDM.'
    'Note that this routine is a general alternative to using the VDM in combination with Realign & Unwarp which estimates and corrects for the combined effects of static and movement-related susceptibility induced distortions. Apply VDM can be used when dynamic distortions are not (well) modelled by Realign & Unwarp (e.g. for fMRI data acquired with R->L phase-encoding direction, high field fMRI data or DTI data).'
    }';
applyvdm.prog = @FieldMap_applyvdm;
applyvdm.vout = @vout_applyvdm;

%--------------------------------------------------------------------------
% fieldmap FieldMap
%--------------------------------------------------------------------------
fieldmap        = cfg_choice;
fieldmap.tag    = 'fieldmap';
fieldmap.name   = 'FieldMap';
fieldmap.help   = {'The FieldMap toolbox generates unwrapped field maps which are converted to voxel displacement maps (VDM) that can be used to unwarp geometrically distorted EPI images. For references and an explanation of the theory behind the field map based unwarping, see FieldMap_principles.man. The resulting VDM files are saved with the prefix vdm and can be applied to images using Apply VDM or in combination with Realign & Unwarp to calculate and correct for the combined effects of static and movement-related susceptibility induced distortions.'};
fieldmap.values = {calculatevdm applyvdm};


%==========================================================================
function out = FieldMap_calculatevdm(job)
for i=1:numel(job.subj)
   out(i) = FieldMap_Run(job.subj(i));
end

%==========================================================================
function dep = vout_calculatevdm(job)
depnum = 1;
for k=1:numel(job.subj)
    for l=1:numel(job.subj(k).session)
        dep(depnum)            = cfg_dep;
        dep(depnum).sname      = sprintf('Voxel displacement map (Subj %d, Session %d)',k,l);
        dep(depnum).src_output = substruct('()',{k},'.','vdmfile','{}',{l});
        dep(depnum).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
        depnum = depnum + 1;
    end
end

%==========================================================================
function dep = vout_applyvdm(job)
for k=1:numel(job.data)   
    if job.roptions.which(1) > 0
        cdep(1)            = cfg_dep;
        cdep(1).sname      = sprintf('VDM corrected images (Sess %d)', k);
        cdep(1).src_output = substruct('.','sess', '()',{k}, '.','rfiles');
        cdep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if k == 1
        dep = cdep;
    else
        dep = [dep cdep];
    end
end
if ~strcmp(job.roptions.which,'<UNDEFINED>') && job.roptions.which(2)
    if exist('dep','var')
        dep(end+1) = cfg_dep;
    else
        dep = cfg_dep;
    end
    dep(end).sname      = 'Mean Image';
    dep(end).src_output = substruct('.','rmean');
    dep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end
