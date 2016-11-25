function realignunwarp = spm_cfg_realignunwarp
% SPM Configuration file for Realign & Unwarp
%__________________________________________________________________________
% Copyright (C) 2005-2016 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_realignunwarp.m 6952 2016-11-25 16:03:13Z guillaume $


%--------------------------------------------------------------------------
% scans Images
%--------------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Images';
scans.help    = {
    'Select scans for this session.'
    'In the coregistration step, the sessions are first realigned to each other, by aligning the first scan from each session to the first scan of the first session.  Then the images within each session are aligned to the first image of the session. The parameter estimation is performed this way because it is assumed (rightly or not) that there may be systematic differences in the images between sessions.'
    }';
scans.filter  = 'image';
scans.ufilter = '.*';
scans.num     = [1 Inf];
scans.preview = @(f) spm_check_registration(char(f));

%--------------------------------------------------------------------------
% pmscan Phase map (vdm* file)
%--------------------------------------------------------------------------
pmscan         = cfg_files;
pmscan.tag     = 'pmscan';
pmscan.name    = 'Phase map (vdm* file)';
pmscan.help    = {
    'Select pre-calculated phase map, or leave empty for no phase correction.'
    'The vdm* file is assumed to be already in alignment with the first scan of the first session.'
    }';
pmscan.filter  = 'image';
pmscan.ufilter = '^vdm5_.*';
pmscan.num     = [0 1];
pmscan.val     = {''};
pmscan.preview = @(f) spm_image('Display',char(f));

%--------------------------------------------------------------------------
% data Session
%--------------------------------------------------------------------------
data         = cfg_branch;
data.tag     = 'data';
data.name    = 'Session';
data.val     = {scans pmscan};
data.help    = {
    'Only add similar session data to a realign+unwarp branch, i.e., choose Data or Data+phase map for all sessions, but don''t use them interchangeably.'
    ''
    'In the coregistration step, the sessions are first realigned to each other, by aligning the first scan from each session to the first scan of the first session.  Then the images within each session are aligned to the first image of the session. The parameter estimation is performed this way because it is assumed (rightly or not) that there may be systematic differences in the images between sessions.'
    }';

%--------------------------------------------------------------------------
% generic Data
%--------------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Data';
generic.help    = {'Data sessions to unwarp.'};
generic.values  = {data};
generic.num     = [1 Inf];

%--------------------------------------------------------------------------
% quality Quality
%--------------------------------------------------------------------------
quality         = cfg_entry;
quality.tag     = 'quality';
quality.name    = 'Quality';
quality.help    = {
    'Quality versus speed trade-off.'
    'Highest quality (1) gives most precise results, whereas lower qualities gives faster realignment. The idea is that some voxels contribute little to the estimation of the realignment parameters. This parameter is involved in selecting the number of voxels that are used.'
    }';
quality.strtype = 'r';
quality.num     = [1 1];
quality.extras  = [0 1];
quality.def     = @(val)spm_get_defaults('realign.estimate.quality', val{:});

%--------------------------------------------------------------------------
% sep Separation
%--------------------------------------------------------------------------
sep         = cfg_entry;
sep.tag     = 'sep';
sep.name    = 'Separation';
sep.help    = {
    'The separation (in mm) between the points sampled in the reference image.'
    'Smaller sampling distances gives more accurate results, but will be slower.'};
sep.strtype = 'r';
sep.num     = [1 1];
sep.def     = @(val)spm_get_defaults('realign.estimate.sep', val{:});

%--------------------------------------------------------------------------
% fwhm Smoothing (FWHM)
%--------------------------------------------------------------------------
fwhm         = cfg_entry;
fwhm.tag     = 'fwhm';
fwhm.name    = 'Smoothing (FWHM)';
fwhm.help    = {
    'The FWHM of the Gaussian smoothing kernel (mm) applied to the images before estimating the realignment parameters.'
    ''
    '    * PET images typically use a 7 mm kernel.'
    ''
    '    * MRI images typically use a 5 mm kernel.'
    }';
fwhm.strtype = 'r';
fwhm.num     = [1 1];
fwhm.def     = @(val)spm_get_defaults('realign.estimate.fwhm', val{:});

%--------------------------------------------------------------------------
% rtm Num Passes
%--------------------------------------------------------------------------
rtm         = cfg_menu;
rtm.tag     = 'rtm';
rtm.name    = 'Num Passes';
rtm.help    = {
    'Register to first: Images are registered to the first image in the series. Register to mean: A two pass procedure is used in order to register the images to the mean of the images after the first realignment.'
    ''
    '    * PET images are typically registered to the mean.'
    ''
    '    * MRI images are typically registered to the first image.'
    }';
rtm.labels  = {
               'Register to first'
               'Register to mean'
}';
rtm.values  = {0 1};
rtm.def     = @(val)spm_get_defaults('unwarp.estimate.rtm', val{:});

%--------------------------------------------------------------------------
% einterp Interpolation
%--------------------------------------------------------------------------
einterp         = cfg_menu;
einterp.tag     = 'einterp';
einterp.name    = 'Interpolation';
einterp.help    = {
    'The method by which the images are sampled when estimating the optimum transformation.'
    'Higher degree interpolation methods provide the better interpolation, but they are slower because they use more neighbouring voxels /* \cite{thevenaz00a,unser93a,unser93b}*/.'
    }';
einterp.labels  = {
                   'Nearest neighbour'
                   'Trilinear'
                   '2nd Degree B-spline'
                   '3rd Degree B-Spline'
                   '4th Degree B-Spline'
                   '5th Degree B-Spline'
                   '6th Degree B-Spline'
                   '7th Degree B-Spline'
}';
einterp.values  = {0 1 2 3 4 5 6 7};
einterp.def     = @(val)spm_get_defaults('realign.estimate.interp', val{:});

%--------------------------------------------------------------------------
% ewrap Wrapping
%--------------------------------------------------------------------------
ewrap         = cfg_menu;
ewrap.tag     = 'ewrap';
ewrap.name    = 'Wrapping';
ewrap.help    = {
    'This indicates which directions in the volumes the values should wrap around in.'
    'These are typically:'
    '* No wrapping - for images that have already been spatially transformed.'
    '* Wrap in Y   - for (un-resliced) MRI where phase encoding is in the Y direction (voxel space).'
    }';
ewrap.labels = {
                'No wrap'
                'Wrap X'
                'Wrap Y'
                'Wrap X & Y'
                'Wrap Z '
                'Wrap X & Z'
                'Wrap Y & Z'
                'Wrap X, Y & Z'
}';
ewrap.values = {[0 0 0] [1 0 0] [0 1 0] [1 1 0] [0 0 1] [1 0 1] [0 1 1]...
                [1 1 1]};
ewrap.def    = @(val)spm_get_defaults('realign.estimate.wrap', val{:});

%--------------------------------------------------------------------------
% weight Weighting
%--------------------------------------------------------------------------
weight         = cfg_files;
weight.tag     = 'weight';
weight.name    = 'Weighting';
weight.val     = {''};
weight.help    = {
    'Optional weighting image to weight each voxel of the reference image differently when estimating the realignment parameters.'
    'The weights are proportional to the inverses of the standard deviations. For example, when there is a lot of extra-brain motion - e.g., during speech, or when there are serious artifacts in a particular region of the images.'
    }';
weight.filter  = 'image';
weight.ufilter = '.*';
weight.num     = [0 1];
weight.preview = @(f) spm_image('Display',char(f));

%--------------------------------------------------------------------------
% eoptions Estimation Options
%--------------------------------------------------------------------------
eoptions         = cfg_branch;
eoptions.tag     = 'eoptions';
eoptions.name    = 'Estimation Options';
eoptions.val     = {quality sep fwhm rtm einterp ewrap weight };
eoptions.help    = {'Various registration options that could be modified to improve the results. Whenever possible, the authors of SPM try to choose reasonable settings, but sometimes they can be improved.'};

%--------------------------------------------------------------------------
% basfcn Basis Functions
%--------------------------------------------------------------------------
basfcn         = cfg_menu;
basfcn.tag     = 'basfcn';
basfcn.name    = 'Basis Functions';
basfcn.help    = {
    'Number of basis functions to use for each dimension.'
    'If the third dimension is left out, the order for that dimension is calculated to yield a roughly equal spatial cut-off in all directions. Default: [12 12 *]'
    }';
basfcn.labels  = {
                  '8x8x*'
                  '10x10x*'
                  '12x12x*'
                  '14x14x*'
}';
basfcn.values  = {[8 8] [10 10] [12 12] [14 14]};
basfcn.def     = @(val)spm_get_defaults('unwarp.estimate.basfcn', val{:});

%--------------------------------------------------------------------------
% regorder Regularisation
%--------------------------------------------------------------------------
regorder         = cfg_menu;
regorder.tag     = 'regorder';
regorder.name    = 'Regularisation';
regorder.help    = {
    'Unwarp looks for the solution that maximises the likelihood (minimises the variance) while simultaneously maximising the smoothness of the estimated field (c.f. Lagrange multipliers). This parameter determines how to balance the compromise between these (i.e. the value of the multiplier). Test it on your own data (if you can be bothered) or go with the defaults. '
    ''
    'Regularisation of derivative fields is based on the regorder''th (spatial) derivative of the field. The choices are 0, 1, 2, or 3.  Default: 1'
    }';
regorder.labels = {
                   '0'
                   '1'
                   '2'
                   '3'
}';
regorder.values = {0 1 2 3};
regorder.def    = @(val)spm_get_defaults('unwarp.estimate.regorder', val{:});

%--------------------------------------------------------------------------
% lambda Reg. Factor
%--------------------------------------------------------------------------
lambda         = cfg_menu;
lambda.tag     = 'lambda';
lambda.name    = 'Reg. Factor';
lambda.help    = {'Regularisation factor. Default: Medium.'};
lambda.labels  = {
                 'A little'
                 'Medium'
                 'A lot'
}';
lambda.values  = {10000 100000 1000000};
lambda.def     = @(val)spm_get_defaults('unwarp.estimate.regwgt', val{:});

%--------------------------------------------------------------------------
% jm Jacobian deformations
%--------------------------------------------------------------------------
jm         = cfg_menu;
jm.tag     = 'jm';
jm.name    = 'Jacobian deformations';
jm.help    = {
    'Option to include Jacobian intensity modulation when estimating the fields.'
    '"Jacobian intensity modulation" refers to the dilution/concentration of intensity that ensue as a consequence of the distortions. Think of a semi-transparent coloured rubber sheet that you hold against a white background. If you stretch a part of the sheet (induce distortions) you will see the colour fading in that particular area. In theory it is a brilliant idea to include also these effects when estimating the field (see e.g. Andersson et al, NeuroImage 20:870-888). In practice for this specific problem it is NOT a good idea. Default: No'
    }';
jm.labels  = {
              'Yes'
              'No'
}';
jm.values  = {1 0};
jm.def     = @(val)spm_get_defaults('unwarp.estimate.jm', val{:});

%--------------------------------------------------------------------------
% fot First-order effects
%--------------------------------------------------------------------------
fot         = cfg_entry;
fot.tag     = 'fot';
fot.name    = 'First-order effects';
fot.help    = {
    'Vector of first order effects to model.'
    ''
    'Theoretically (ignoring effects of shimming) one would expect the field to depend only on subject out-of-plane rotations. Hence the default choice ("Pitch and Roll", i.e., [4 5]). Go with that unless you have very good reasons to do otherwise'
    ''
    'Movements to be modelled are referred to by number. 1= x translation; 2= y translation; 3= z translation 4 = x rotation,  5 = y rotation and 6 = z rotation.'
    ''
    'To model pitch & roll enter: [4 5]'
    ''
    'To model all movements enter: [1:6]'
    ''
    'Otherwise enter a customised set of movements to model'
    }';
fot.strtype = 'n';
fot.num     = [1 Inf];
fot.def     = @(val)spm_get_defaults('unwarp.estimate.foe', val{:});

%--------------------------------------------------------------------------
% sot Second-order effects
%--------------------------------------------------------------------------
sot         = cfg_entry;
sot.tag     = 'sot';
sot.name    = 'Second-order effects';
sot.help    = {
    'List of second order terms to model second derivatives of.'
    'This is entered as  a vector of movement parameters similar to first order effects, or leave blank for NONE.'
    ''
    'Movements to be modelled are referred to by number:'
    ''
    '1= x translation; 2= y translation; 3= z translation 4 = x rotation,  5 = y rotation and 6 = z rotation.'
    ''
    'To model the interaction of pitch & roll enter: [4 5]'
    ''
    'To model all movements enter: [1:6]'
    ''
    'The vector will be expanded into an n x 2 matrix of effects. For example [4 5] will be expanded to:'
    ''
    '[ 4 4'
    ''
    '  4 5'
    ''
    '  5 5 ]'
    }';
sot.strtype = 'n';
sot.num     = [Inf Inf];
sot.def     = @(val)spm_get_defaults('unwarp.estimate.soe', val{:});

%--------------------------------------------------------------------------
% uwfwhm Smoothing for unwarp (FWHM)
%--------------------------------------------------------------------------
uwfwhm         = cfg_entry;
uwfwhm.tag     = 'uwfwhm';
uwfwhm.name    = 'Smoothing for unwarp (FWHM)';
uwfwhm.help    = {'FWHM (mm) of smoothing filter applied to images prior to estimation of deformation fields.'};
uwfwhm.strtype = 'r';
uwfwhm.num     = [1 1];
uwfwhm.def     = @(val)spm_get_defaults('unwarp.estimate.fwhm', val{:});

%--------------------------------------------------------------------------
% rem Re-estimate movement params
%--------------------------------------------------------------------------
rem         = cfg_menu;
rem.tag     = 'rem';
rem.name    = 'Re-estimate movement params';
rem.help    = {'Re-estimation means that movement-parameters should be re-estimated at each unwarping iteration.'};
rem.labels  = {
              'Yes'
              'No'
}';
rem.values  = {1 0};
rem.def     = @(val)spm_get_defaults('unwarp.estimate.rem', val{:});

%--------------------------------------------------------------------------
% noi Number of Iterations
%--------------------------------------------------------------------------
noi         = cfg_entry;
noi.tag     = 'noi';
noi.name    = 'Number of Iterations';
noi.help    = {'Maximum number of iterations.'};
noi.strtype = 'n';
noi.num     = [1 1];
noi.def     = @(val)spm_get_defaults('unwarp.estimate.noi', val{:});

%--------------------------------------------------------------------------
% expround Taylor expansion point
%--------------------------------------------------------------------------
expround         = cfg_menu;
expround.tag     = 'expround';
expround.name    = 'Taylor expansion point';
expround.help    = {
    'Point in position space to perform Taylor-expansion around.'
    'Choices are (''First'', ''Last'' or ''Average''). ''Average'' should (in principle) give the best variance reduction. If a field-map acquired before the time-series is supplied then expansion around the ''First'' MIGHT give a slightly better average geometric fidelity.'
    };
expround.labels = {
                   'Average'
                   'First'
                   'Last'
}';
expround.values = {
                   'Average'
                   'First'
                   'Last'
}';
expround.def     = @(val)spm_get_defaults('unwarp.estimate.expround', val{:});

%--------------------------------------------------------------------------
% uweoptions Unwarp Estimation Options
%--------------------------------------------------------------------------
uweoptions      = cfg_branch;
uweoptions.tag  = 'uweoptions';
uweoptions.name = 'Unwarp Estimation Options';
uweoptions.val  = {basfcn regorder lambda jm fot sot uwfwhm rem noi expround};
uweoptions.help = {'Various registration & unwarping estimation options.'};

%--------------------------------------------------------------------------
% uwwhich Reslices images (unwarp)?
%--------------------------------------------------------------------------
uwwhich         = cfg_menu;
uwwhich.tag     = 'uwwhich';
uwwhich.name    = 'Resliced images (unwarp)?';
uwwhich.help    = {
    'Specify the images to reslice.'
    'All Images (1..n) '
    '  This reslices and unwarps all the images. '
    ' '
    'All Images + Mean Image '
    '  In addition to reslicing the images, it also creates a mean of the resliced images.'
    }';
uwwhich.labels = {
                  ' All Images (1..n)'
                  ' All Images + Mean Image'
}';
uwwhich.values = {[2 0] [2 1]};
uwwhich.def    = @(val)spm_get_defaults('realign.write.which', val{:});

%--------------------------------------------------------------------------
% rinterp Interpolation
%--------------------------------------------------------------------------
rinterp         = cfg_menu;
rinterp.tag     = 'rinterp';
rinterp.name    = 'Interpolation';
rinterp.help    = {
    'The method by which the images are sampled when being written in a different space.'
    'Nearest Neighbour is fastest, but not recommended for image realignment. Trilinear Interpolation is probably OK for PET, but not so suitable for fMRI because higher degree interpolation generally gives better results/* \cite{thevenaz00a,unser93a,unser93b}*/. Although higher degree methods provide better interpolation, but they are slower because they use more neighbouring voxels.'
    }';
rinterp.labels = {
                  'Nearest neighbour'
                  'Trilinear'
                  '2nd Degree B-spline '
                  '3rd Degree B-Spline'
                  '4th Degree B-Spline'
                  '5th Degree B-Spline '
                  '6th Degree B-Spline'
                  '7th Degree B-Spline'
}';
rinterp.values = {0 1 2 3 4 5 6 7};
rinterp.def    = @(val)spm_get_defaults('realign.write.interp', val{:});

%--------------------------------------------------------------------------
% wrap Wrapping
%--------------------------------------------------------------------------
wrap         = cfg_menu;
wrap.tag     = 'wrap';
wrap.name    = 'Wrapping';
wrap.help    = {
    'This indicates which directions in the volumes the values should wrap around in.'
    'These are typically:'
    '* No wrapping - for images that have already been spatially transformed.'
    '* Wrap in Y  - for (un-resliced) MRI where phase encoding is in the Y direction (voxel space).'
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
prefix.help    = {'Specify the string to be prepended to the filenames of the smoothed image file(s). Default prefix is ''u''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.def     = @(val)spm_get_defaults('unwarp.write.prefix', val{:});

%--------------------------------------------------------------------------
% uwroptions Unwarp Reslicing Options
%--------------------------------------------------------------------------
uwroptions         = cfg_branch;
uwroptions.tag     = 'uwroptions';
uwroptions.name    = 'Unwarp Reslicing Options';
uwroptions.val     = {uwwhich rinterp wrap mask prefix};
uwroptions.help    = {'Various registration & unwarping estimation options.'};

%--------------------------------------------------------------------------
% realignunwarp Realign & Unwarp
%--------------------------------------------------------------------------
realignunwarp         = cfg_exbranch;
realignunwarp.tag     = 'realignunwarp';
realignunwarp.name    = 'Realign & Unwarp';
realignunwarp.val     = {generic eoptions uweoptions uwroptions};
realignunwarp.help    = {
    'Within-subject registration and unwarping of time series.'
    ''
    'The realignment part of this routine realigns a time-series of images acquired from the same subject using a least squares approach and a 6 parameter (rigid body) spatial transformation.  The first image in the list specified by the user is used as a reference to which all subsequent scans are realigned. The reference scan does not have to the the first chronologically and it may be wise to chose a "representative scan" in this role.'
    ''
    'The aim is primarily to remove movement artefact in fMRI and PET time-series (or more generally longitudinal studies). This affects the header of each of the input images which contains details about the voxel-to-world mapping. The details of the transformation are displayed in the results window as plots of translation and rotation. A set of realignment parameters are saved for each session, named rp_*.txt.'
    ''
    'In the coregistration step, the sessions are first realigned to each other, by aligning the first scan from each session to the first scan of the first session.  Then the images within each session are aligned to the first image of the session. The parameter estimation is performed this way because it is assumed (rightly or not) that there may be systematic differences in the images between sessions.'
    ''
    'See spm_uw_estimate.m for a detailed description of the implementation.'
    ''
    'Even after realignment there is considerable variance in fMRI time series that covary with, and is most probably caused by, subject movements/* \cite{ja_geometric}*/. It is also the case that this variance is typically large compared to experimentally induced variance. Anyone interested can include the estimated movement parameters as covariates in the design matrix, and take a look at an F-contrast encompassing those columns. It is quite dramatic. The result is loss of sensitivity, and if movements are correlated to task specificity. I.e. we may mistake movement induced variance for true activations. The problem is well known, and several solutions have been suggested. A quite pragmatic (and conservative) solution is to include the estimated movement parameters (and possibly squared) as covariates in the design matrix. Since we typically have loads of degrees of freedom in fMRI we can usually afford this. The problems occur when movements are correlated with the task, since the strategy above will discard "good" and "bad" variance alike (i.e. remove also "true" activations).'
    ''
    'The "covariate" strategy described above was predicated on a model where variance was assumed to be caused by "spin history" effects, but will work pretty much equally good/bad regardless of what the true underlying cause is. Others have assumed that the residual variance is caused mainly by errors introduced by the interpolation kernel in the resampling step of the realignment. One has tried to solve this through higher order resampling (huge Sinc kernels, or k-space resampling). Unwarp is based on a different hypothesis regarding the residual variance. EPI images are not particularly faithful reproductions of the object, and in particular there are severe geometric distortions in regions where there is an air-tissue interface (e.g. orbitofrontal cortex and the anterior medial temporal lobes). In these areas in particular the observed image is a severely warped version of reality, much like a funny mirror at a fair ground. When one moves in front of such a mirror ones image will distort in different ways and ones head may change from very elongated to seriously flattened. If we were to take digital snapshots of the reflection at these different positions it is rather obvious that realignment will not suffice to bring them into a common space.'
    ''
    'The situation is similar with EPI images, and an image collected for a given subject position will not be identical to that collected at another. We call this effect susceptibility-by-movement interaction. Unwarp is predicated on the assumption that the susceptibility-by-movement interaction is responsible for a sizable part of residual movement related variance.'
    ''
    'Assume that we know how the deformations change when the subject changes position (i.e. we know the derivatives of the deformations with respect to subject position). That means that for a given time series and a given set of subject movements we should be able to predict the "shape changes" in the object and the ensuing variance in the time series. It also means that, in principle, we should be able to formulate the inverse problem, i.e. given the observed variance (after realignment) and known (estimated) movements we should be able to estimate how deformations change with subject movement. We have made an attempt at formulating such an inverse model, and at solving for the "derivative fields". A deformation field can be thought of as little vectors at each position in space showing how that particular location has been deflected. A "derivative field" is then the rate of change of those vectors with respect to subject movement. Given these "derivative fields" we should be able to remove the variance caused by the susceptibility-by-movement interaction. Since the underlying model is so restricted we would also expect experimentally induced variance to be preserved. Our experiments have also shown this to be true.'
    ''
    'In theory it should be possible to estimate also the "static" deformation field, yielding an unwarped (to some true geometry) version of the time series. In practise that doesn''t really seem to work. Hence, the method deals only with residual movement related variance induced by the susceptibility-by-movement interaction. This means that the time-series will be undistorted to some "average distortion" state rather than to the true geometry. If one wants additionally to address the issue of anatomical fidelity one should combine Unwarp with a measured fieldmap.'
    ''
    'The description above can be thought of in terms of a Taylor expansion of the field as a function of subject movement. Unwarp alone will estimate the first (and optionally second, see below) order terms of this expansion. It cannot estimate the zeroth order term (the distortions common to all scans in the time series) since that doesn''t introduce (almost) any variance in the time series. The measured fieldmap takes the role of the zeroth order term. Refer to the FieldMap toolbox and the documents FieldMap.man and FieldMap_principles.man for a description of how to obtain fieldmaps in the format expected by Unwarp.'
    ''
    'If we think of the field as a function of subject movement it should in principle be a function of six variables since rigid body movement has six degrees of freedom. However, the physics of the problem tells us that the field should not depend on translations nor on rotation in a plane perpendicular to the magnetic flux. Hence it should in principle be sufficient to model the field as a function of out-of-plane rotations (i.e. pitch and roll). One can object to this in terms of the effects of shimming (object no longer immersed in a homogenous field) that introduces a dependence on all movement parameters. In addition SPM/Unwarp cannot really tell if the transversal slices it is being passed are really perpendicular to the flux or not. In practice it turns out thought that it is never (at least we haven''t seen any case) necessary to include more than Pitch and Roll. This is probably because the individual movement parameters are typically highly correlated anyway, which in turn is probably because most heads that we scan are attached to a neck around which rotations occur. On the subject of Taylor expansion we should mention that there is the option to use a second-order expansion (through the defaults) interface. This implies estimating also the rate-of-change w.r.t. to some movement parameter of the rate-of-change of the field w.r.t. some movement parameter (colloquially known as a second derivative). It can be quite interesting to watch (and it is amazing that it is possible) but rarely helpful/necessary.'
    ''
    'In the defaults there is also an option to include Jacobian intensity modulation when estimating the fields. "Jacobian intensity modulation" refers to the dilution/concentration of intensity that ensue as a consequence of the distortions. Think of a semi-transparent coloured rubber sheet that you hold against a white background. If you stretch a part of the sheet (induce distortions) you will see the colour fading in that particular area. In theory it is a brilliant idea to include also these effects when estimating the field (see e.g. Andersson et al, NeuroImage 20:870-888). In practice for this specific problem it is NOT a good idea.'
    ''
    'It should be noted that this is a method intended to correct data afflicted by a particular problem. If there is little movement in your data to begin with this method will do you little good. If on the other hand there is appreciable movement in your data (>1deg) it will remove some of that unwanted variance. If, in addition, movements are task related it will do so without removing all your "true" activations. The method attempts to minimise total (across the image volume) variance in the data set. It should be realised that while (for small movements) a rather limited portion of the total variance is removed, the susceptibility-by-movement interaction effects are quite localised to "problem" areas. Hence, for a subset of voxels in e.g. frontal-medial and orbitofrontal cortices and parts of the temporal lobes the reduction can be quite dramatic (>90). The advantages of using Unwarp will also depend strongly on the specifics of the scanner and sequence by which your data has been acquired. When using the latest generation scanners distortions are typically quite small, and distortion-by-movement interactions consequently even smaller. A small check list in terms of distortions is '
    'a) Fast gradients->short read-out time->small distortions '
    'b) Low field (i.e. <3T)->small field changes->small distortions '
    'c) Low res (64x64)->short read-out time->small distortions '
    'd) SENSE/SMASH->short read-out time->small distortions '
    'If you can tick off all points above chances are you have minimal distortions to begin with and Unwarp might not be of use to you.'
    }';
realignunwarp.prog     = @spm_run_realignunwarp;
realignunwarp.vout     = @vout_realignunwarp;
realignunwarp.modality = {'PET','FMRI'};


%==========================================================================
function dep = vout_realignunwarp(job)
for k=1:numel(job.data)
    cdep(1)            = cfg_dep;
    cdep(1).sname      = sprintf('Realignment Param File (Sess %d)', k);
    cdep(1).src_output = substruct('.','sess', '()',{k}, '.','rpfile');
    cdep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
    cdep(2)            = cfg_dep;
    cdep(2).sname      = sprintf('Unwarp Params File (Sess %d)', k);
    cdep(2).src_output = substruct('.','sess', '()',{k}, '.','dsfile');
    cdep(2).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});
    if job.uwroptions.uwwhich(1) == 2
        cdep(3)            = cfg_dep;
        cdep(3).sname      = sprintf('Unwarped Images (Sess %d)', k);
        cdep(3).src_output = substruct('.','sess', '()',{k}, '.','uwrfiles');
        cdep(3).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if k == 1
        dep = cdep;
    else
        dep = [dep cdep];
    end
end

if job.uwroptions.uwwhich(2)
    dep(end+1)          = cfg_dep;
    dep(end).sname      = 'Unwarped Mean Image';
    dep(end).src_output = substruct('.','meanuwr');
    dep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end
