function [SPM] = spm_fmri_spm_ui(SPM)
% Setting up the general linear model for fMRI time-series
% FORMAT [SPM] = spm_fmri_spm_ui(SPM)
%
% creates SPM with the following fields
%
%       xY: [1x1 struct] - data structure
%    nscan: [double]     - vector of scans per session
%      xBF: [1x1 struct] - Basis function structure   (see spm_fMRI_design)
%     Sess: [1x1 struct] - Session structure          (see spm_fMRI_design)
%       xX: [1x1 struct] - Design matrix structure    (see spm_fMRI_design)
%      xGX: [1x1 struct] - Global variate structure
%      xVi: [1x1 struct] - Non-sphericity structure
%       xM: [1x1 struct] - Masking structure
%    xsDes: [1x1 struct] - Design description structure
%
%
%     SPM.xY
%             P: [n x ? char]       - filenames
%            VY: [n x 1 struct]     - filehandles
%            RT: Repeat time
%
%    SPM.xGX
%
%       iGXcalc: {'none'|'Scaling'} - Global normalization option
%       sGXcalc: 'mean voxel value' - Calculation method
%        sGMsca: 'session specific' - Grand mean scaling
%            rg: [n x 1 double]     - Global estimate
%            GM: 100                - Grand mean
%           gSF: [n x 1 double]     - Global scaling factor
%
%    SPM.xVi
%            Vi: {[n x n sparse]..} - covariance components
%          form: {'none'|'AR(1)'}   - form of non-sphericity
%
%     SPM.xM
%             T: [n x 1 double]     - Masking index
%            TH: [n x 1 double]     - Threshold
%             I: 0
%            VM:                    - Mask filehandles
%            xs: [1x1 struct]       - cellstr description
%
%__________________________________________________________________________
%
% spm_fmri_spm_ui configures the design matrix, data specification and
% filtering that specify the ensuing statistical analysis. These arguments
% are passed to spm_spm that then performs the actual parameter estimation.
%
% The design matrix defines the experimental design and the nature of
% hypothesis testing to be implemented.  The design matrix has one row for
% each scan and one column for each effect or explanatory variable (e.g.
% regressor or stimulus function).  The parameters are estimated in a least
% squares sense using the general linear model.  Specific profiles within
% these parameters are tested using a linear compound or contrast with the
% T or F statistic.  The resulting statistical map constitutes an SPM.  The
% SPM{T}/{F} is then characterized in terms of focal or regional
% differences by assuming that (under the null hypothesis) the components
% of the SPM (i.e. residual fields) behave as smooth stationary Gaussian
% fields.
%
% spm_fmri_spm_ui allows you to (i) specify a statistical model in terms of
% a design matrix, (ii) associate some data with a pre-specified design [or
% (iii) specify both the data and design] and then proceed to estimate the
% parameters of the model.
% Inferences can be made about the ensuing parameter estimates (at a first
% or fixed-effect level) in the results section, or they can be re-entered
% into a second (random-effect) level analysis by treating the session or
% subject-specific [contrasts of] parameter estimates as new summary data.
% Inferences at any level are obtained by specifying appropriate T or F
% contrasts in the results section to produce SPMs and tables of p values
% and statistics.
%
% spm_fmri_spm calls spm_fMRI_design which allows you to configure a design
% matrix in terms of events or epochs.
%
% spm_fMRI_design allows you to build design matrices with separable
% session-specific partitions.  Each partition may be the same (in which
% case it is only necessary to specify it once) or different.  Responses
% can be either event- or epoch related, The only distinction is the
% duration of the underlying input or stimulus function. Mathematically
% they are both modelled by convolving a series of delta (stick) or box
% functions (u), indicating the onset of an event or epoch with a set of
% basis functions.  These basis functions model the hemodynamic
% convolution, applied by the brain, to the inputs.  This convolution can
% be first-order or a generalized convolution modelled to second order (if
% you specify the Volterra option). [The same inputs are used by the
% hemodynamic model or or dynamic causal models which model the convolution
% explicitly in terms of hidden state variables (see spm_hdm_ui and
% spm_dcm_ui).]
% Basis functions can be used to plot estimated responses to single events
% once the parameters (i.e. basis function coefficients) have been
% estimated.  The importance of basis functions is that they provide a
% graceful transition between simple fixed response models (like the
% box-car) and finite impulse response (FIR) models, where there is one
% basis function for each scan following an event or epoch onset.  The nice
% thing about basis functions, compared to FIR models, is that data
% sampling and stimulus presentation does not have to be synchronized
% thereby allowing a uniform and unbiased sampling of peri-stimulus time.
%
% Event-related designs may be stochastic or deterministic.  Stochastic
% designs involve one of a number of trial-types occurring with a specified
% probably at successive intervals in time.  These probabilities can be
% fixed (stationary designs) or time-dependent (modulated or non-stationary
% designs).  The most efficient designs obtain when the probabilities of
% every trial type are equal.
% A critical issue in stochastic designs is whether to include null events
% If you wish to estimate the evoke response to a specific event type (as
% opposed to differential responses) then a null event must be included
% (even if it is not modelled explicitly).
%
% The choice of basis functions depends upon the nature of the inference
% sought.  One important consideration is whether you want to make
% inferences about compounds of parameters (i.e.  contrasts).  This is the
% case if (i) you wish to use a SPM{T} to look separately at activations
% and deactivations or (ii) you with to proceed to a second (random-effect)
% level of analysis.  If this is the case then (for event-related studies)
% use a canonical hemodynamic response function (HRF) and derivatives with
% respect to latency (and dispersion).  Unlike other bases, contrasts of
% these effects have a physical interpretation and represent a parsimonious
% way of characterising event-related responses.  Bases such as a Fourier
% set require the SPM{F} for inference.
%
% See spm_fMRI_design for more details about how designs are specified.
%
% Serial correlations in fast fMRI time-series are dealt with as described
% in spm_spm.  At this stage you need to specify the filtering that will be
% applied to the data (and design matrix) to give a generalized least
% squares (GLS) estimate of the parameters required.  This filtering is
% important to ensure that the GLS estimate is efficient and that the error
% variance is estimated in an unbiased way.
%
% The serial correlations will be estimated with a ReML (restricted maximum
% likelihood) algorithm using an autoregressive AR(1) model during
% parameter estimation.  This estimate assumes the same correlation
% structure for each voxel, within each session.  The ReML estimates are
% then used to correct for non-sphericity during inference by adjusting the
% statistics and degrees of freedom appropriately.  The discrepancy between
% estimated and actual intrinsic (i.e. prior to filtering) correlations are
% greatest at low frequencies.  Therefore specification of the high-pass
% filter is particularly important.
%
% High-pass filtering is implemented at the level of the filtering matrix K
% (as opposed to entering as confounds in the design matrix).  The default
% cut-off period is 128 seconds. Use 'explore design' to ensure this
% cut-off is not removing too much experimental variance.
% Note that high-pass filtering uses a residual forming matrix (i.e. it is
% not a convolution) and is simply to a way to remove confounds without
% estimating their parameters explicitly.  The constant term is also
% incorporated into this filter matrix.
%
%--------------------------------------------------------------------------
% Refs:
%
% Friston KJ, Holmes A, Poline J-B, Grasby PJ, Williams SCR, Frackowiak
% RSJ & Turner R (1995) Analysis of fMRI time-series revisited. NeuroImage
% 2:45-53
%
% Worsley KJ and Friston KJ (1995) Analysis of fMRI time-series revisited -
% again. NeuroImage 2:178-181
%
% Friston KJ, Frith CD, Frackowiak RSJ, & Turner R (1995) Characterising
% dynamic brain responses with fMRI: A multivariate approach NeuroImage -
% 2:166-172
%
% Frith CD, Turner R & Frackowiak RSJ (1995) Characterising evoked
% hemodynamics with fMRI Friston KJ, NeuroImage 2:157-165
%
% Josephs O, Turner R and Friston KJ (1997) Event-related fMRI, Hum. Brain
% Map. 5:243-248
%
%__________________________________________________________________________
% Copyright (C) 1994-2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fmri_spm_ui.m 7738 2019-12-02 12:45:37Z guillaume $


SVNid = '$Rev: 7738 $';

%==========================================================================
% - D E S I G N   M A T R I X
%==========================================================================
SPM   = spm_fMRI_design(SPM);

%-Session and scan number
%--------------------------------------------------------------------------
nscan = SPM.nscan;
nsess = length(nscan);


%==========================================================================
% - D E S I G N   P A R A M E T E R S
%==========================================================================
SPM.SPMid = spm('FnBanner',mfilename,SVNid);

%-High-pass filtering
%==========================================================================

%-Low frequency confounds
%--------------------------------------------------------------------------
try
    HParam     = [SPM.xX.K(:).HParam];
    if length(HParam) == 1 
        HParam = repmat(HParam,1,nsess);
    end
catch
    error('High-pass filter not specified.');
end

%-Create and set filter structure
%--------------------------------------------------------------------------
for  i = 1:nsess
    K(i) = struct('HParam', HParam(i),...
                  'row',    SPM.Sess(i).row,...
                  'RT',     SPM.xY.RT);
end
SPM.xX.K = spm_filter(K);


%-Intrinsic autocorrelations (Vi) for non-sphericity ReML estimation
%==========================================================================
try
    cVi  = SPM.xVi.form;
catch
    error('Serial correlations not specified.');
end

%-Create Vi structure
%--------------------------------------------------------------------------

if ~ischar(cVi)                 % AR coefficient specified
    %----------------------------------------------------------------------
    SPM.xVi.Vi = spm_Ce(nscan,cVi(1));
    cVi        = ['AR( ' sprintf('%0.1f ',cVi) ')'];

else
    switch lower(cVi)
        
        case {'i.i.d', 'none'}  %  xVi.V is i.i.d
            %--------------------------------------------------------------
            SPM.xVi.V  = speye(sum(nscan));
            cVi        = 'i.i.d';

        case 'fast'
            %--------------------------------------------------------------
            SPM.xVi.Vi = spm_Ce('fast',nscan,SPM.xY.RT);
            cVi        = upper(cVi);
            
        otherwise               % otherwise assume AR(0.2) in xVi.Vi
            %--------------------------------------------------------------
            SPM.xVi.Vi = spm_Ce('ar',nscan,0.2);
            cVi        = 'AR(0.2)';
    end
end

SPM.xVi.form = cVi;


%-Design description - for saving and display
%==========================================================================
for i     = 1:nsess, ntr(i) = length(SPM.Sess(i).U); end
Fstr      = sprintf('[min] Cutoff: %d {s}',min([SPM.xX.K(:).HParam]));
SPM.xsDes = struct(...
    'Basis_functions',      SPM.xBF.name,...
    'Number_of_sessions',   sprintf('%d',nsess),...
    'Trials_per_session',   sprintf('%-3d',ntr),...
    'Interscan_interval',   sprintf('%0.2f {s}',SPM.xY.RT),...
    'High_pass_Filter',     Fstr);


%-Return if design-only specification
%==========================================================================
try, SPM.xY.P; catch, return; end


%==========================================================================
% - C O N F I G U R E   D E S I G N
%==========================================================================

spm('Pointer','Watch');

%-Get image files
%==========================================================================

%-Map files
%--------------------------------------------------------------------------
fprintf('%-40s: ','Mapping files')                                      %-#
VY    = spm_data_hdr_read(SPM.xY.P);
fprintf('%30s\n','...done')                                             %-#

%-Check internal consistency of images
%--------------------------------------------------------------------------
spm_check_orientations(VY);

%-Place mapped files in xY
%--------------------------------------------------------------------------
SPM.xY.VY = VY;


%-Compute Global variate
%==========================================================================
GM    = 100;
q     = length(VY);
g     = zeros(q,1);
fprintf('%-40s: ','Calculating globals')                                %-#
spm_progress_bar('Init',q,'Calculating globals');
if spm_mesh_detect(VY)
    for i = 1:q
        dat = spm_data_read(VY(i));
        g(i) = mean(dat(~isnan(dat)));
        spm_progress_bar('Set',i)
    end
else
    for i = 1:q
        g(i) = spm_global(VY(i));
        spm_progress_bar('Set',i)
    end
end
spm_progress_bar('Clear');
fprintf('%30s\n','...done')                                             %-#

%-Scale if specified (otherwise session specific grand mean scaling)
%--------------------------------------------------------------------------
gSF   = GM./g;
if strcmpi(SPM.xGX.iGXcalc,'none')
    for i = 1:nsess
        gSF(SPM.Sess(i).row) = GM./mean(g(SPM.Sess(i).row));
    end
end

%-Apply gSF to memory-mapped scalefactors to implement scaling
%--------------------------------------------------------------------------
for i = 1:q
    SPM.xY.VY(i).pinfo(1:2,:) = SPM.xY.VY(i).pinfo(1:2,:) * gSF(i);
    if spm_mesh_detect(VY)
        SPM.xY.VY(i).private.private.data{1}.data.scl_slope = ...
            SPM.xY.VY(i).private.private.data{1}.data.scl_slope * gSF(i);
        SPM.xY.VY(i).private.private.data{1}.data.scl_inter = ...
            SPM.xY.VY(i).private.private.data{1}.data.scl_inter * gSF(i);
    else
        SPM.xY.VY(i).private.dat.scl_slope = ...
            SPM.xY.VY(i).private.dat.scl_slope * gSF(i);
        SPM.xY.VY(i).private.dat.scl_inter = ...
            SPM.xY.VY(i).private.dat.scl_inter * gSF(i);
    end
end

%-Place global variates in xGX
%--------------------------------------------------------------------------
SPM.xGX.sGXcalc = 'mean voxel value';
SPM.xGX.sGMsca  = 'session specific';
SPM.xGX.rg      = g;
SPM.xGX.GM      = GM;
SPM.xGX.gSF     = gSF;


%-Masking
%==========================================================================

%-Masking threshold, as proportion of globals
%--------------------------------------------------------------------------
try
    gMT = SPM.xM.gMT;
catch
    gMT = spm_get_defaults('mask.thresh');
end
TH = g.*gSF*gMT;

%-Place masking structure in xM
%--------------------------------------------------------------------------
SPM.xM = struct(...
    'T',   ones(q,1),...
    'TH',  TH,...
    'gMT', gMT,...
    'I',   0,...
    'VM',  {[]},...
    'xs',  struct('Masking','analysis threshold'));


%-Design description - for saving and display
%==========================================================================
xs = struct(...
    'Global_calculation',   SPM.xGX.sGXcalc,...
    'Grand_mean_scaling',   SPM.xGX.sGMsca,...
    'Global_normalisation', SPM.xGX.iGXcalc);
for fn=(fieldnames(xs))', SPM.xsDes.(fn{1}) = xs.(fn{1}); end


%==========================================================================
% - S A V E   A N D   D I S P L A Y
%==========================================================================

%-Save SPM.mat
%--------------------------------------------------------------------------
%if ~nargout
    fprintf('%-40s: ','Saving SPM configuration')                       %-#
    fmt = spm_get_defaults('mat.format');
    s = whos('SPM');
    if s.bytes > 2147483647, fmt = '-v7.3'; end
    save('SPM.mat', 'SPM', fmt);
    fprintf('%30s\n','...SPM.mat saved')                                %-#
%end


%-Display design report
%--------------------------------------------------------------------------
if ~spm('CmdLine') && ~isempty(spm_figure('FindWin','Graphics'))
    fprintf('%-40s: ','Design reporting')                               %-#
    try,   fname = reshape(cellstr(SPM.xY.P),size(SPM.xY.VY));
    catch, fname = {}; end
    spm_DesRep('DesMtx',SPM.xX,fname,SPM.xsDes)
    fprintf('%30s\n','...done')                                         %-#
end

spm('Pointer','Arrow');
