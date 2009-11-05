function fieldmap = tbx_cfg_fieldmap
% MATLABBATCH Configuration file for toolbox 'FieldMap'
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: tbx_cfg_fieldmap.m 3536 2009-11-05 15:00:32Z chloe $

addpath(fullfile(spm('dir'),'toolbox','FieldMap'));

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Default values that are common to all fieldmap jobs
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
% et Echo times [short TE long TE]
% ---------------------------------------------------------------------
et         = cfg_entry;
et.tag     = 'et';
et.name    = 'Echo times [short TE long TE]';
et.help    = {'Enter the short and long echo times (in ms) of the data used to acquire the field map.'};
et.strtype = 'e';
et.num     = [1  2];
% ---------------------------------------------------------------------
% maskbrain Mask brain
% ---------------------------------------------------------------------
maskbrain         = cfg_menu;
maskbrain.tag     = 'maskbrain';
maskbrain.name    = 'Mask brain';
maskbrain.help    = {
                     'Select masking or no masking of the brain. If masking is selected,'
                     'the magnitude image is used to generate a mask of the brain.'
}';
maskbrain.labels = {
                    'Mask brain'
                    'No brain masking'
}';
maskbrain.values{1} = 1;
maskbrain.values{2} = 0;
% ---------------------------------------------------------------------
% blipdir Blip direction
% ---------------------------------------------------------------------
blipdir         = cfg_menu;
blipdir.tag     = 'blipdir';
blipdir.name    = 'Blip direction';
blipdir.help    = {'Enter the blip direction. This is the polarity of the phase-encode blips describing the direction in which k-space is traversed along the y-axis during EPI acquisition with respect to the coordinate system used in SPM. In this coordinate system, the phase encode direction corresponds with the y-direction and is deﬁned as positive from the posterior to the anterior of the head.'};
blipdir.labels = {
                  '-1'
                  '1'
}';
blipdir.values{1} = -1;
blipdir.values{2} = 1;
% ---------------------------------------------------------------------
% tert Total EPI readout time
% ---------------------------------------------------------------------
tert         = cfg_entry;
tert.tag     = 'tert';
tert.name    = 'Total EPI readout time';
tert.help    = {
                'Enter the total EPI readout time (in ms). This is the time taken to '
                'acquire all of the phase encode steps required to cover k-space (ie one image slice). '
                'For example, if the EPI sequence has 64 phase encode steps, the total readout time is '
                'the time taken to acquire 64 echoes, e.g. '
                'total readout time = number of echoes × echo spacing. '
                'This time does not include i) the duration of the excitation, ii) the delay between, '
                'the excitation and the start of the acquisition or iii) time for fat saturation etc.'
}';
tert.strtype = 'e';
tert.num     = [1  1];
% ---------------------------------------------------------------------
% epifm EPI-based field map?
% ---------------------------------------------------------------------
epifm         = cfg_menu;
epifm.tag     = 'epifm';
epifm.name    = 'EPI-based field map?';
epifm.help    = {'Select non-EPI or EPI based field map. The field map data may be acquired using a non-EPI sequence (typically a gradient echo sequence) or an EPI sequence. The processing will be slightly different for the two cases. If using an EPI-based ﬁeld map, the resulting Voxel Displacement Map will be inverted since the ﬁeld map was acquired in distorted space.'};
epifm.labels = {
                'non-EPI'
                'EPI'
}';
epifm.values{1} = 0;
epifm.values{2} = 1;
% ---------------------------------------------------------------------
% ajm Jacobian modulation?
% ---------------------------------------------------------------------
ajm         = cfg_menu;
ajm.tag     = 'ajm';
ajm.name    = 'Jacobian modulation?';
ajm.help    = {'Select whether or not to use Jacobian modulation. This will adjust the intensities of voxels that have been stretched or compressed but in general is not recommended for EPI distortion correction'};
ajm.labels = {
              'Do not use'
              'Use'
}';
ajm.values{1} = 0;
ajm.values{2} = 1;
ajm.def       = @(val)pm_get_defaults('DO_JACOBIAN_MODULATION', val{:});
% ---------------------------------------------------------------------
% method Unwrapping method
% ---------------------------------------------------------------------
method         = cfg_menu;
method.tag     = 'method';
method.name    = 'Unwrapping method';
method.help    = {'Select method for phase unwrapping'};
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
method.def     = @(val)pm_get_defaults('UNWRAPPING_METHOD', val{:});
% ---------------------------------------------------------------------
% fwhm FWHM
% ---------------------------------------------------------------------
fwhm         = cfg_entry;
fwhm.tag     = 'fwhm';
fwhm.name    = 'FWHM';
fwhm.help    = {'FWHM of Gaussian filter used to implement weighted smoothing of unwrapped maps.'};
fwhm.strtype = 'e';
fwhm.num     = [1  1];
fwhm.def     = @(val)pm_get_defaults('FWHM', val{:});
% ---------------------------------------------------------------------
% pad pad
% ---------------------------------------------------------------------
pad         = cfg_entry;
pad.tag     = 'pad';
pad.name    = 'pad';
pad.help    = {'Size of padding kernel if required.'};
pad.strtype = 'e';
pad.num     = [1  1];
pad.def     = @(val)pm_get_defaults('PAD', val{:});
% ---------------------------------------------------------------------
% ws Weighted smoothing
% ---------------------------------------------------------------------
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
% ---------------------------------------------------------------------
% uflags uflags
% ---------------------------------------------------------------------
uflags         = cfg_branch;
uflags.tag     = 'uflags';
uflags.name    = 'uflags';
uflags.val     = {method fwhm pad ws };
uflags.help    = {'Different options for phase unwrapping and field map processing'};
% ---------------------------------------------------------------------
% template Template image for brain masking
% ---------------------------------------------------------------------
template         = cfg_files;
template.tag     = 'template';
template.name    = 'Template image for brain masking';
template.help    = {'Select template file for segmentation to create brain mask'};
template.filter = 'nii';
template.ufilter = '.*';
template.num     = [1 1];
template.def     = @(val)pm_get_defaults('MFLAGS.TEMPLATE', val{:});
% ---------------------------------------------------------------------
% fwhm FWHM
% ---------------------------------------------------------------------
fwhm         = cfg_entry;
fwhm.tag     = 'fwhm';
fwhm.name    = 'FWHM';
fwhm.help    = {'FWHM of Gaussian filter for smoothing brain mask.'};
fwhm.strtype = 'e';
fwhm.num     = [1  1];
fwhm.def     = @(val)pm_get_defaults('MFLAGS.FWHM', val{:});
% ---------------------------------------------------------------------
% nerode Number of erosions
% ---------------------------------------------------------------------
nerode         = cfg_entry;
nerode.tag     = 'nerode';
nerode.name    = 'Number of erosions';
nerode.help    = {'Number of erosions used to create brain mask.'};
nerode.strtype = 'e';
nerode.num     = [1  1];
nerode.def     = @(val)pm_get_defaults('MFLAGS.NERODE', val{:});
% ---------------------------------------------------------------------
% ndilate Number of dilations
% ---------------------------------------------------------------------
ndilate         = cfg_entry;
ndilate.tag     = 'ndilate';
ndilate.name    = 'Number of dilations';
ndilate.help    = {'Number of dilations used to create brain mask.'};
ndilate.strtype = 'e';
ndilate.num     = [1  1];
ndilate.def     = @(val)pm_get_defaults('MFLAGS.NDILATE', val{:});
% ---------------------------------------------------------------------
% thresh Threshold
% ---------------------------------------------------------------------
thresh         = cfg_entry;
thresh.tag     = 'thresh';
thresh.name    = 'Threshold';
thresh.help    = {'Threshold used to create brain mask from segmented data.'};
thresh.strtype = 'e';
thresh.num     = [1  1];
thresh.def     = @(val)pm_get_defaults('MFLAGS.THRESH', val{:});
% ---------------------------------------------------------------------
% reg Regularization
% ---------------------------------------------------------------------
reg         = cfg_entry;
reg.tag     = 'reg';
reg.name    = 'Regularization';
reg.help    = {'Regularization value used in the segmentation. A larger value helps the segmentation to converge.'};
reg.strtype = 'e';
reg.num     = [1  1];
reg.def     = @(val)pm_get_defaults('MFLAGS.REG', val{:});
% ---------------------------------------------------------------------
% mflags mflags
% ---------------------------------------------------------------------
mflags         = cfg_branch;
mflags.tag     = 'mflags';
mflags.name    = 'mflags';
mflags.val     = {template fwhm nerode ndilate thresh reg };
mflags.help    = {'Different options used for the segmentation and creation of the brain mask.'};
% ---------------------------------------------------------------------
% defaultsval Defaults values
% ---------------------------------------------------------------------
defaultsval         = cfg_branch;
defaultsval.tag     = 'defaultsval';
defaultsval.name    = 'Defaults values';
defaultsval.val     = {et maskbrain blipdir tert epifm ajm uflags mflags };
defaultsval.help    = {'Defaults values'};
% ---------------------------------------------------------------------
% defaultsfile Defaults File
% ---------------------------------------------------------------------
defaultsfile         = cfg_files;
defaultsfile.tag     = 'defaultsfile';
defaultsfile.name    = 'Defaults File';
defaultsfile.help    = {'Select the ''pm_defaults*.m'' file containing the parameters for the field map data. Please make sure that the parameters defined in the defaults file are correct for your field map and EPI sequence. To create your own defaults file, either edit the distributed version and/or save it with the name ''pm_defaults_yourname.m''.'};
defaultsfile.filter  = 'm';
[deffilepath, tmp]   = fileparts(mfilename('fullpath'));
defaultsfile.dir     = deffilepath;
defaultsfile.ufilter = '^pm_defaults.*\.m$';
defaultsfile.num     = [1 1];
defaultsfile.def     = @(val)pm_get_defaults('defaultsfilename', val{:});
% ---------------------------------------------------------------------
% defaults FieldMap defaults
% ---------------------------------------------------------------------
defaults         = cfg_choice;
defaults.tag     = 'defaults';
defaults.name    = 'FieldMap defaults';
defaults.help    = {'FieldMap default values can be entered as a file or set of values.'};
defaults.values  = {defaultsval defaultsfile };
% ---------------------------------------------------------------------
% epi Select EPI to Unwarp
% ---------------------------------------------------------------------
epi         = cfg_files;
epi.tag     = 'epi';
epi.name    = 'Select EPI to Unwarp';
epi.help    = {'Select an image to distortion correct. The corrected image will be saved with the prefix u. The original and the distortion corrected images can be displayed for comparison.'};
epi.filter = 'image';
epi.ufilter = '.*';
epi.num     = [1 1];
% ---------------------------------------------------------------------
% session Session
% ---------------------------------------------------------------------
session         = cfg_branch;
session.tag     = 'session';
session.name    = 'Session';
session.val     = {epi };
session.help    = {'Data for this session.'};
% ---------------------------------------------------------------------
% generic EPI Sessions
% ---------------------------------------------------------------------
generic1         = cfg_repeat;
generic1.tag     = 'generic';
generic1.name    = 'EPI Sessions';
generic1.help    = {'If a single VDM file will be used for multiple sessions, select the first EPI in each session. A copy of the VDM file will be matched to the first EPI in each session and save with a seprate name.'};
generic1.values  = {session };
generic1.num     = [1 Inf];
% ---------------------------------------------------------------------
% matchvdm Match VDM to EPI?
% ---------------------------------------------------------------------
matchvdm         = cfg_menu;
matchvdm.tag     = 'matchvdm';
matchvdm.name    = 'Match VDM to EPI?';
matchvdm.help    = {'Match VDM file to EPI image. This option will coregister the field map data to the selected EPI before doing distortion correction.'};
matchvdm.labels = {
                   'match vdm'
                   'none'
}';
matchvdm.values{1} = 1;
matchvdm.values{2} = 0;
% ---------------------------------------------------------------------
% sessname Name extension for session specific vdm files
% ---------------------------------------------------------------------
sessname         = cfg_entry;
sessname.tag     = 'sessname';
sessname.name    = 'Name extension for session specific vdm files';
sessname.help    = {'This will be the name extension followed by an incremented integer for session specific vdm files.'};
sessname.strtype = 's';
sessname.num     = [1  Inf];
sessname.def     = @(val)pm_get_defaults('sessname', val{:});
% ---------------------------------------------------------------------
% writeunwarped Write unwarped EPI?
% ---------------------------------------------------------------------
writeunwarped         = cfg_menu;
writeunwarped.tag     = 'writeunwarped';
writeunwarped.name    = 'Write unwarped EPI?';
writeunwarped.help    = {'Write out distortion corrected EPI image. The image is saved with the prefix u.'};
writeunwarped.labels = {
                        'write unwarped EPI'
                        'none'
}';
writeunwarped.values{1} = 1;
writeunwarped.values{2} = 0;
% ---------------------------------------------------------------------
% anat Select anatomical image for comparison
% ---------------------------------------------------------------------
anat         = cfg_files;
anat.tag     = 'anat';
anat.name    = 'Select anatomical image for comparison';
anat.help    = {'Select an anatomical image for comparison with the distortion corrected EPI.'};
anat.filter = 'image';
anat.ufilter = '.*';
anat.num     = [0 1];
% ---------------------------------------------------------------------
% matchanat Match anatomical image to EPI?
% ---------------------------------------------------------------------
matchanat         = cfg_menu;
matchanat.tag     = 'matchanat';
matchanat.name    = 'Match anatomical image to EPI?';
matchanat.help    = {'Match the anatomical image to the distortion corrected EPI. This can be done to allow for visual inspection and comparison of the distortion corrected EPI.'};
matchanat.labels = {
                    'match anat'
                    'none'
}';
matchanat.values{1} = 1;
matchanat.values{2} = 0;

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Different kinds of input jobs
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Presubtracted phase and magnitude
% ---------------------------------------------------------------------
% phase Phase Image
% ---------------------------------------------------------------------
phase         = cfg_files;
phase.tag     = 'phase';
phase.name    = 'Phase Image';
phase.help    = {'Select a single phase image. This should be the result from the subtraction of two phase images (where the subtraction is usually done automatically by the scanner software). The phase image will be scaled between +/- PI.'};
phase.filter = 'image';
phase.ufilter = '.*';
phase.num     = [1 1];
% ---------------------------------------------------------------------
% magnitude Magnitude Image
% ---------------------------------------------------------------------
magnitude         = cfg_files;
magnitude.tag     = 'magnitude';
magnitude.name    = 'Magnitude Image';
magnitude.help    = {'Select a single magnitude image'};
magnitude.filter = 'image';
magnitude.ufilter = '.*';
magnitude.num     = [1 1];

% ---------------------------------------------------------------------
% subj Subject
% ---------------------------------------------------------------------
subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {phase magnitude defaults generic1 matchvdm sessname writeunwarped anat matchanat };
subj.help    = {'Data for this subject or field map session.'};
% ---------------------------------------------------------------------
% generic Data
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Data';
generic.help    = {'Subjects or sessions for which individual field map data has been acquired.'};
generic.values  = {subj };
generic.num     = [1 Inf];
% ---------------------------------------------------------------------
% presubphasemag Presubtracted Phase and Magnitude Data
% ---------------------------------------------------------------------
presubphasemag         = cfg_exbranch;
presubphasemag.tag     = 'presubphasemag';
presubphasemag.name    = 'Presubtracted Phase and Magnitude Data';
presubphasemag.val     = {generic };
presubphasemag.help    = {'Calculate a voxel displacement map (VDM) from presubtracted phase and magnitude field map data. This option expects a single magnitude image and a single phase image resulting from the subtraction of two phase images (where the subtraction is usually done automatically by the scanner software). The phase image will be scaled between +/- PI.'};
presubphasemag.prog = @fieldmap_presubphasemag;
presubphasemag.vout = @vout;

% ---------------------------------------------------------------------
% Real and imaginary
% ---------------------------------------------------------------------
% shortreal Short Echo Real Image
% ---------------------------------------------------------------------
shortreal         = cfg_files;
shortreal.tag     = 'shortreal';
shortreal.name    = 'Short Echo Real Image';
shortreal.help    = {'Select short echo real image'};
shortreal.filter = 'image';
shortreal.ufilter = '.*';
shortreal.num     = [1 1];
% ---------------------------------------------------------------------
% shortimag Short Echo Imaginary Image
% ---------------------------------------------------------------------
shortimag         = cfg_files;
shortimag.tag     = 'shortimag';
shortimag.name    = 'Short Echo Imaginary Image';
shortimag.help    = {'Select short echo imaginary image'};
shortimag.filter = 'image';
shortimag.ufilter = '.*';
shortimag.num     = [1 1];
% ---------------------------------------------------------------------
% longreal Long Echo Real Image
% ---------------------------------------------------------------------
longreal         = cfg_files;
longreal.tag     = 'longreal';
longreal.name    = 'Long Echo Real Image';
longreal.help    = {'Select long echo real image'};
longreal.filter = 'image';
longreal.ufilter = '.*';
longreal.num     = [1 1];
% ---------------------------------------------------------------------
% longimag Long Echo Imaginary Image
% ---------------------------------------------------------------------
longimag         = cfg_files;
longimag.tag     = 'longimag';
longimag.name    = 'Long Echo Imaginary Image';
longimag.help    = {'Select long echo imaginary image'};
longimag.filter = 'image';
longimag.ufilter = '.*';
longimag.num     = [1 1];
% ---------------------------------------------------------------------
% subj Subject
% ---------------------------------------------------------------------
subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {shortreal shortimag longreal longimag defaults generic1 matchvdm sessname writeunwarped anat matchanat };
subj.help    = {'Data for this subject or field map session.'};
% ---------------------------------------------------------------------
% generic Data
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Data';
generic.help    = {'Subjects or sessions for which individual field map data has been acquired.'};
generic.values  = {subj };
generic.num     = [1 Inf];
% ---------------------------------------------------------------------
% realimag Real and Imaginary Data
% ---------------------------------------------------------------------
realimag         = cfg_exbranch;
realimag.tag     = 'realimag';
realimag.name    = 'Real and Imaginary Data';
realimag.val     = {generic };
realimag.help    = {'Calculate a voxel displacement map (VDM) from real and imaginary field map data. This option expects two real and imaginary pairs of data of two different echo times.'};
realimag.prog = @fieldmap_realimag;
realimag.vout = @vout;

% ---------------------------------------------------------------------
% Phase and magnitude (4 files)
% ---------------------------------------------------------------------
% shortphase Short Echo Phase Image
% ---------------------------------------------------------------------
shortphase         = cfg_files;
shortphase.tag     = 'shortphase';
shortphase.name    = 'Short Echo Phase Image';
shortphase.help    = {'Select short echo phase image'};
shortphase.filter = 'image';
shortphase.ufilter = '.*';
shortphase.num     = [1 1];
% ---------------------------------------------------------------------
% shortmag Short Echo Magnitude Image
% ---------------------------------------------------------------------
shortmag         = cfg_files;
shortmag.tag     = 'shortmag';
shortmag.name    = 'Short Echo Magnitude Image';
shortmag.help    = {'Select short echo magnitude image'};
shortmag.filter = 'image';
shortmag.ufilter = '.*';
shortmag.num     = [1 1];
% ---------------------------------------------------------------------
% longphase Long Echo Phase Image
% ---------------------------------------------------------------------
longphase         = cfg_files;
longphase.tag     = 'longphase';
longphase.name    = 'Long Echo Phase Image';
longphase.help    = {'Select long echo phase image'};
longphase.filter = 'image';
longphase.ufilter = '.*';
longphase.num     = [1 1];
% ---------------------------------------------------------------------
% longmag Long Echo Magnitude Image
% ---------------------------------------------------------------------
longmag         = cfg_files;
longmag.tag     = 'longmag';
longmag.name    = 'Long Echo Magnitude Image';
longmag.help    = {'Select long echo magnitude image'};
longmag.filter = 'image';
longmag.ufilter = '.*';
longmag.num     = [1 1];
% ---------------------------------------------------------------------
% subj Subject
% ---------------------------------------------------------------------
subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {shortphase shortmag longphase longmag defaults generic1 matchvdm sessname writeunwarped anat matchanat };
subj.help    = {'Data for this subject or field map session.'};
% ---------------------------------------------------------------------
% generic Data
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Data';
generic.help    = {'Subjects or sessions for which individual field map data has been acquired.'};
generic.values  = {subj };
generic.num     = [1 Inf];
% ---------------------------------------------------------------------
% phasemag Phase and Magnitude Data
% ---------------------------------------------------------------------
phasemag         = cfg_exbranch;
phasemag.tag     = 'phasemag';
phasemag.name    = 'Phase and Magnitude Data';
phasemag.val     = {generic };
phasemag.help    = {'Calculate a voxel displacement map (VDM) from double phase and magnitude field map data. This option expects two phase and magnitude pairs of data of two different echo times.'};
phasemag.prog = @fieldmap_phasemag;
phasemag.vout = @vout;

% ---------------------------------------------------------------------
% Precalculated fieldmap
% ---------------------------------------------------------------------
% precalcfieldmap Precalculated field map
% ---------------------------------------------------------------------
precalcfieldmap1         = cfg_files;
precalcfieldmap1.tag     = 'precalcfieldmap';
precalcfieldmap1.name    = 'Precalculated field map';
precalcfieldmap1.help    = {'Select a precalculated field map. This should be a processed field map (ie phase unwrapped, masked if necessary and scaled to Hz) , for example as generated by the FieldMap toolbox and stored as an fpm_* file.'};
precalcfieldmap1.filter = 'image';
precalcfieldmap1.ufilter = '.*';
precalcfieldmap1.num     = [1 1];
% ---------------------------------------------------------------------
% magfieldmap Select magnitude image in same space as fieldmap
% ---------------------------------------------------------------------
magfieldmap         = cfg_files;
magfieldmap.tag     = 'magfieldmap';
magfieldmap.name    = 'Select magnitude image in same space as fieldmap';
magfieldmap.help    = {'Select magnitude image which is in the same space as the field map to do matching to EPI.'};
magfieldmap.filter = 'image';
magfieldmap.ufilter = '.*';
magfieldmap.num     = [0 1];

% ---------------------------------------------------------------------
% subj Subject
% ---------------------------------------------------------------------
subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {precalcfieldmap1 magfieldmap defaults generic1 matchvdm sessname writeunwarped anat matchanat };
subj.help    = {'Data for this subject or field map session.'};
% ---------------------------------------------------------------------
% generic Data
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Data';
generic.help    = {'Subjects or sessions for which individual field map data has been acquired.'};
generic.values  = {subj };
generic.num     = [1 Inf];
% ---------------------------------------------------------------------
% precalcfieldmap Precalculated FieldMap (in Hz)
% ---------------------------------------------------------------------
precalcfieldmap         = cfg_exbranch;
precalcfieldmap.tag     = 'precalcfieldmap';
precalcfieldmap.name    = 'Precalculated FieldMap (in Hz)';
precalcfieldmap.val     = {generic };
precalcfieldmap.help    = {'Calculate a voxel displacement map (VDM) from a precalculated field map. This option expects a processed field map (ie phase unwrapped, masked if necessary and scaled to Hz). Precalculated field maps can be generated by the FieldMap toolbox and stored as fpm_* files.'};
precalcfieldmap.prog = @fieldmap_precalcfieldmap;
precalcfieldmap.vout = @vout;

% ---------------------------------------------------------------------
% fieldmap FieldMap
% ---------------------------------------------------------------------
fieldmap         = cfg_choice;
fieldmap.tag     = 'fieldmap';
fieldmap.name    = 'FieldMap';
fieldmap.help    = {
                    'The FieldMap toolbox generates unwrapped field maps which are converted to voxel displacement maps (VDM) that can be used to unwarp geometrically distorted EPI images. For references and an explantion of the theory behind the field map based unwarping, see FieldMap_principles.man. The resulting VDM files are saved with the prefix vdm and can be used in combination with Realign & Unwarp to calculate and correct for the combined effects of static and movement-related susceptibility induced distortions.'
                    ''
}';
fieldmap.values  = {presubphasemag realimag phasemag precalcfieldmap };

%------------------------------------------------------------------------
function out=fieldmap_presubphasemag(job)
for i=1:numel(job.subj),
   out(i)=FieldMap_Run(job.subj(i));
end
%------------------------------------------------------------------------
function out=fieldmap_realimag(job)
for i=1:numel(job.subj),
   out(i)=FieldMap_Run(job.subj(i));
end
%------------------------------------------------------------------------
function out=fieldmap_phasemag(job)
for i=1:numel(job.subj),
   out(i)=FieldMap_Run(job.subj(i));
end
%------------------------------------------------------------------------
function out=fieldmap_precalcfieldmap(job)
for i=1:numel(job.subj),
   out(i)=FieldMap_Run(job.subj(i));
end

%------------------------------------------------------------------------
function dep = vout(job)

depnum=1;
for k=1:numel(job.subj)
 for l=1:numel(job.subj(k).session)
   dep(depnum)            = cfg_dep;
   dep(depnum).sname      = sprintf('Voxel displacement map (Subj %d, Session %d)',k,l);
   dep(depnum).src_output = substruct('()',{k},'.','vdmfile','{}',{l});
   dep(depnum).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
   depnum=depnum+1;
 end
end;



        
        


