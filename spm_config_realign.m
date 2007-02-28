function opts = spm_config_realign
% Configuration file for realign jobs
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_config_realign.m 751 2007-02-28 10:56:59Z volkmar $


%_______________________________________________________________________

quality.type    = 'entry';
quality.name    = 'Quality';
quality.tag     = 'quality';
quality.strtype = 'r';
quality.num     = [1 1];
quality.def     = 'realign.estimate.quality';
quality.extras  = [0 1];
quality.help    = {[...
'Quality versus speed trade-off.  Highest quality (1) gives most ',...
'precise results, whereas lower qualities gives faster realignment. ',...
'The idea is that some voxels contribute little to the estimation of ',...
'the realignment parameters. This parameter is involved in selecting ',...
'the number of voxels that are used.']};

%------------------------------------------------------------------------

weight.type   = 'files';
weight.name   = 'Weighting';
weight.tag    = 'weight';
weight.filter = 'image';
weight.num    = [0 1];
weight.val    = {{}};
weight.help   = {[...
'The option of providing a weighting image to weight each voxel ',...
'of the reference image differently when estimating the realignment ',...
'parameters.  The weights are proportional to the inverses of the ',...
'standard deviations. ',...
'For example, when there is a lot of extra-brain motion - e.g., during ',...
'speech, or when there are serious artifacts in a particular region of ',...
'the images.']};

%------------------------------------------------------------------------

interp.type   = 'menu';
interp.name   = 'Interpolation';
interp.tag    = 'interp';
interp.labels = {'Trilinear (1st Degree)','2nd Degree B-Spline',...
'3rd Degree B-Spline ','4th Degree B-Spline','5th Degree B-Spline',...
'6th Degree B-Spline','7th Degree B-Spline'};
interp.values = {1,2,3,4,5,6,7};
interp.def   = 'realign.estimate.interp';
interp.help  = {...
['The method by which the images are sampled when estimating the optimum transformation. ',...
'Higher degree interpolation methods provide the better interpolation, but they are slower ',...
'because they use more neighbouring voxels /* \cite{thevenaz00a,unser93a,unser93b}*/. ']};
 
%------------------------------------------------------------------------
 
wrap.type   = 'menu';
wrap.name   = 'Wrapping';
wrap.tag    = 'wrap';
wrap.labels = {'No wrap','Wrap X','Wrap Y','Wrap X & Y','Wrap Z',...
                'Wrap X & Z','Wrap Y & Z','Wrap X, Y & Z'};
wrap.values = {[0 0 0],[1 0 0],[0 1 0],[1 1 0],[0 0 1],[1 0 1],[0 1 1],[1 1 1]};
wrap.def    = 'realign.estimate.wrap';
wrap.help   = {...
['This indicates which directions in the volumes the values should wrap around in.  ',...
'For example, in MRI scans, the images wrap around in the phase encode direction, ',...
'so (e.g.) the subject''s nose may poke into the back of the subject''s head. ',...
'These are typically:'],...
['    No wrapping - for PET or images that have already ',...
 '                  been spatially transformed. Also the recommended option if ',...
 '                  you are not really sure.'],...
['    Wrap in  Y  - for (un-resliced) MRI where phase encoding ',...
 '                  is in the Y direction (voxel space).']};
 
%------------------------------------------------------------------------

fwhm.type    = 'entry';
fwhm.name    = 'Smoothing (FWHM)';
fwhm.tag     = 'fwhm';
fwhm.num     = [1 1];
fwhm.def     = 'realign.estimate.fwhm';
fwhm.strtype = 'e';
p1           = [...
'The FWHM of the Gaussian smoothing kernel (mm) applied to the ',...
'images before estimating the realignment parameters.'];
p2           = '    * PET images typically use a 7 mm kernel.';
p3           = '    * MRI images typically use a 5 mm kernel.';
fwhm.help    = {p1,'',p2,'',p3};

%------------------------------------------------------------------------

sep.type = 'entry';
sep.name = 'Separation';
sep.tag  = 'sep';
sep.num  = [1 1];
sep.strtype = 'e';
sep.def  = 'realign.estimate.sep';
%sep.val = {4};
sep.help = {[...
'The separation (in mm) between the points sampled in the ',...
'reference image.  Smaller sampling distances gives more accurate ',...
'results, but will be slower.']};

%------------------------------------------------------------------------

rtm.type   = 'menu';
rtm.name   = 'Num Passes';
rtm.tag    = 'rtm';
rtm.labels = {'Register to first','Register to mean'};
rtm.values = {0,1};
rtm.def    = 'realign.estimate.rtm';
p1         = [...
'Register to first: Images are registered to the first image in the series. ',...
'Register to mean:   A two pass procedure is used in order to register the ',...
'images to the mean of the images after the first realignment.'];
p2         = [...
'PET images are typically registered to the mean. This is because PET data are ',...
'more noisy than fMRI and there are fewer of them, so time is less of an issue.'];
p3         = [...
'MRI images are typically registered to the first image.  The more accurate way ',...
'would be to use a two pass procedure, but this probably wouldn''t improve the results ',...
'so much and would take twice as long to run.'];
rtm.help    = {p1,'',p2,'',p3};

%------------------------------------------------------------------------

% global defaults
% if ~isempty(defaults) && isfield(defaults,'modality') ...
%                       && strcmp(lower(defaults.modality),'pet'),
%     fwhm.val = {7};
%     rtm.val  = {1};
% else
%     fwhm.val = {5};
%     rtm.val  = {0};
% end;

eoptions.type = 'branch';
eoptions.name = 'Estimation Options';
eoptions.tag  = 'eoptions';
eoptions.val  = {quality,sep,fwhm,rtm,interp,wrap,weight};
eoptions.help = {[...
'Various registration options. ',...
'If in doubt, simply keep the default values.']};

%------------------------------------------------------------------------

which.type = 'menu';
which.name = 'Resliced images';
which.tag  = 'which';
which.labels = {' All Images (1..n)','Images 2..n',...
                ' All Images + Mean Image',' Mean Image Only'};
which.values = {[2 0],[1 0],[2 1],[0 1]};
which.val    = {[2 1]};
which.help = {...
['All Images (1..n) : ',...
'  This reslices all the images - including the first image selected ',...
'  - which will remain in its original position.'],...
'',...
['Images 2..n : ',...
'   Reslices images 2..n only. Useful for if you wish to reslice ',...
'   (for example) a PET image to fit a structural MRI, without ',...
'   creating a second identical MRI volume.'],...
'',...
['All Images + Mean Image : ',...
'   In addition to reslicing the images, it also creates a mean of the ',...
'   resliced image.'],...
'',...
['Mean Image Only : ',...
'   Creates the mean resliced image only.']};

%------------------------------------------------------------------------
 
interp.type = 'menu';
interp.name = 'Interpolation';
interp.tag  = 'interp';
interp.labels = {'Nearest neighbour','Trilinear','2nd Degree B-Spline',...
'3rd Degree B-Spline','4th Degree B-Spline','5th Degree B-Spline',...
'6th Degree B-Spline','7th Degree B-Spline','Fourier Interpolation'};
interp.values = {0,1,2,3,4,5,6,7,Inf};
interp.def  = 'realign.write.interp';
interp.help = {...
['The method by which the images are sampled when being written in a ',...
'different space.',...
'Nearest Neighbour is fastest, but not recommended for image realignment. ',...
'Bilinear Interpolation is probably OK for PET, but not so suitable for fMRI because ',...
'higher degree interpolation generally gives better results/* \cite{thevenaz00a,unser93a,unser93b}*/. ',...
'Although higher degree methods provide better interpolation, but they are slower ',...
'because they use more neighbouring voxels. ',...
'Fourier Interpolation/* \cite{eddy96,cox99}*/ is another option, but note that it ',...
'is only implemented for purely rigid body transformations.  Voxel sizes must all be ',...
'identical and isotropic.']};
 
%------------------------------------------------------------------------
 
wrap.type = 'menu';
wrap.name = 'Wrapping';
wrap.tag  = 'wrap';
wrap.labels = {'No wrap','Wrap X','Wrap Y','Wrap X & Y','Wrap Z',...
'Wrap X & Z','Wrap Y & Z','Wrap X, Y & Z'};
wrap.values = {[0 0 0],[1 0 0],[0 1 0],[1 1 0],[0 0 1],[1 0 1],[0 1 1],[1 1 1]};
wrap.def    = 'realign.write.wrap';
wrap.help = {...
['This indicates which directions in the volumes the values should wrap around in.  ',...
'For example, in MRI scans, the images wrap around in the phase encode direction, ',...
'so (e.g.) the subject''s nose may poke into the back of the subject''s head. ',...
'These are typically:'],...
['    No wrapping - for PET or images that have already ',...
'                  been spatially transformed.'],...
['    Wrap in  Y  - for (un-resliced) MRI where phase encoding ',...
'                  is in the Y direction (voxel space).']};
 
%------------------------------------------------------------------------
 
mask.type = 'menu';
mask.name = 'Masking';
mask.tag  = 'mask';
mask.labels = {'Mask images','Dont mask images'};
mask.values = {1,0};
mask.def    = 'realign.write.mask';
mask.help = {[...
'Because of subject motion, different images are likely to have different ',...
'patterns of zeros from where it was not possible to sample data. ',...
'With masking enabled, the program searches through the whole time series ',...
'looking for voxels which need to be sampled from outside the original ',...
'images. Where this occurs, that voxel is set to zero for the whole set ',...
'of images (unless the image format can represent NaN, in which case ',...
'NaNs are used where possible).']};
 
%------------------------------------------------------------------------
 
roptions.type = 'branch';
roptions.name = 'Reslice Options';
roptions.tag  = 'roptions';
roptions.val  = {which,interp,wrap,mask};
roptions.help = {'Various reslicing options. If in doubt, simply keep the default values.'};

%------------------------------------------------------------------------

scans.type = 'files';
scans.name = 'Session';
scans.tag  = 'data';
scans.num  = [1 Inf];
scans.filter = 'image';
scans.help   = {[...
'Select scans for this session. ',...
'In the coregistration step, the sessions are first realigned to ',...
'each other, by aligning the first scan from each session to the ',...
'first scan of the first session.  Then the images within each session ',...
'are aligned to the first image of the session. ',...
'The parameter estimation is performed this way because it is assumed ',...
'(rightly or not) that there may be systematic differences ',...
'in the images between sessions.']};

%------------------------------------------------------------------------

data.type   = 'repeat';
data.name   = 'Data';
% data.tag  = 'data';
data.values = {scans};
data.num    = [1 Inf];
data.help   = {[...
'Add new sessions for this subject. ',...
'In the coregistration step, the sessions are first realigned to ',...
'each other, by aligning the first scan from each session to the ',...
'first scan of the first session.  Then the images within each session ',...
'are aligned to the first image of the session. ',...
'The parameter estimation is performed this way because it is assumed ',...
'(rightly or not) that there may be systematic differences ',...
'in the images between sessions.']};

%------------------------------------------------------------------------

est.type = 'branch';
est.name = 'Realign: Estimate';
est.tag  = 'estimate';
est.val  = {data,eoptions};
est.prog = @estimate;
est.vfiles = @vfiles_estimate;

p1 = [...
'This routine realigns a time-series of images acquired from the same ',...
'subject using a least squares approach and a 6 parameter (rigid body) ',...
'spatial transformation/* \cite{friston95a}*/.  The first image in the list specified by the ',...
'user is used as a reference to which all subsequent scans are realigned. ',...
'The reference scan does not have to the the first chronologically and ',...
'it may be wise to chose a "representative scan" in this role.'];

p2 = [...
'The aim is primarily to remove movement artefact in fMRI and PET ',...
'time-series (or more generally longitudinal studies). ',...
'The headers are modified for each of the input images, such that. ',...
'they reflect the relative orientations of the data. ',...
'The details of the transformation are displayed in the results window ',...
'as plots of translation and rotation. ',...
'A set of realignment parameters are saved for each session, named ',...
'rp_*.txt. These can be modelled as confounds within the general linear model/* \cite{friston95a}*/.'];

est.help = {p1,'',p2};

%------------------------------------------------------------------------

scans.type = 'files';
scans.name = 'Images';
scans.tag  = 'data';
scans.num  = [1 Inf];
scans.filter = 'image';
scans.help   = {'Select scans to reslice to match the first.'};

%------------------------------------------------------------------------

write.type = 'branch';
write.name = 'Realign: Reslice';
write.tag  = 'write';
write.val  = {scans,roptions};
write.help = {[...
'This function reslices a series of registered images such that they ',...
'match the first image selected voxel-for-voxel. The resliced images ',...
'are named the same as the originals, except that they are prefixed ',...
'by ''r''.']};
write.prog   = @reslice;
write.vfiles = @vfiles_reslice;
%------------------------------------------------------------------------

estwrit.type = 'branch';
estwrit.name = 'Realign: Estimate & Reslice';
estwrit.tag  = 'estwrite';
estwrit.val  = {data,eoptions,roptions};
p1 = [...
'This routine realigns a time-series of images acquired from the same ',...
'subject using a least squares approach and a 6 parameter (rigid body)',...
'spatial transformation/* \cite{friston95a}*/.  The first image in the list specified by the ',...
'user is used as a reference to which all subsequent scans are realigned. ',...
'The reference scan does not have to the the first chronologically and ',...
'it may be wise to chose a "representative scan" in this role.'];

p2 = [...
'The aim is primarily to remove movement artefact in fMRI and PET ',...
'time-series (or more generally longitudinal studies) /* \cite{ashburner97bir}*/. ',...
'The headers are modified for each of the input images, such that. ',...
'they reflect the relative orientations of the data. ',...
'The details of the transformation are displayed in the results window ',...
'as plots of translation and rotation. ',...
'A set of realignment parameters are saved for each session, named ',...
'rp_*.txt. After realignment, the images are resliced ',...
'such that they match the first image selected voxel-for-voxel. ',...
'The resliced images are named the same as the originals, except that ',...
'they are prefixed by ''r''.'];

estwrit.help = {p1,'',p2};
estwrit.prog   = @estwrite_fun;
estwrit.vfiles = @vfiles_estwrit;

%------------------------------------------------------------------------

opts.type = 'repeat';
opts.name = 'Realign';
opts.tag  = 'realign';
opts.values = {est,write,estwrit};
opts.num  = [1 Inf];
opts.modality = {'PET','FMRI','VBM'};
opts.help = {...
'Within-subject registration of image time series.'};

%------------------------------------------------------------------------

%------------------------------------------------------------------------

return;

%------------------------------------------------------------------------

%------------------------------------------------------------------------
function estimate(varargin)
job           = varargin{1};
P             = {};
for i=1:length(job.data),
	P{i}  = strvcat(job.data{i});
end;
flags.quality = job.eoptions.quality;
flags.fwhm    = job.eoptions.fwhm;
flags.sep     = job.eoptions.sep;
flags.rtm     = job.eoptions.rtm;
flags.PW      = strvcat(job.eoptions.weight);
flags.interp  = job.eoptions.interp;
flags.wrap    = job.eoptions.wrap;
spm_realign(P,flags);
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function reslice(varargin)
job          = varargin{1};
P            = strvcat(job.data);
flags.mask   = job.roptions.mask;
flags.mean   = job.roptions.which(2);
flags.interp = job.roptions.interp;
flags.which  = job.roptions.which(1);
flags.wrap   = job.roptions.wrap;
spm_reslice(P,flags);
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function estwrite_fun(varargin)
job           = varargin{1};
P             = {};
for i=1:length(job.data),
    P{i} = strvcat(job.data{i});
end;
flags.quality = job.eoptions.quality;
flags.fwhm    = job.eoptions.fwhm;
flags.sep     = job.eoptions.sep;
flags.rtm     = job.eoptions.rtm;
flags.PW      = strvcat(job.eoptions.weight);
flags.interp  = job.eoptions.interp;
flags.wrap    = job.eoptions.wrap;
spm_realign(P,flags);

P            = strvcat(P);
flags.mask   = job.roptions.mask;
flags.mean   = job.roptions.which(2);
flags.interp = job.roptions.interp;
flags.which  = job.roptions.which(1);
flags.wrap   = job.roptions.wrap;
spm_reslice(P,flags);
return;

%------------------------------------------------------------------------
 
%------------------------------------------------------------------------
function vf = vfiles_reslice(job)
P = job.data;
if numel(P)>0 && iscell(P{1}),
    P = cat(1,P{:});
end;

switch job.roptions.which(1),
case 0,
    vf = {};
case 1,
    vf = cell(numel(P)-1,1);
    for i=1:length(vf),
        [pth,nam,ext,num] = spm_fileparts(P{i+1});
        vf{i} = fullfile(pth,['r', nam, ext, num]);
    end;
otherwise,
    vf = cell(numel(P),1);
    for i=1:length(vf),
        [pth,nam,ext,num] = spm_fileparts(P{i});
        vf{i} = fullfile(pth,['r', nam, ext, num]);
    end;
end;
if job.roptions.which(2),
    [pth,nam,ext,num] = spm_fileparts(P{1});
    vf = {vf{:}, fullfile(pth,['mean', nam, ext, num])};
end;

%------------------------------------------------------------------------
 
%------------------------------------------------------------------------
function vf = vfiles_estimate(job)
P = job.data;
vf = {};
if numel(P) > 0
    if ~iscell(P{1})
	P = {P};
    end;
    for k = 1:numel(P)
	[pth,nam,ext,num] = spm_fileparts(P{k}{1});
	vf{k} = fullfile(pth, sprintf('rp_%s.txt', nam));
    end;
end;

%------------------------------------------------------------------------
 
%------------------------------------------------------------------------
function vf = vfiles_estwrit(job)
vf1 = vfiles_estimate(job);
vf2 = vfiles_reslice(job);
vf = {vf1{:}, vf2{:}};

