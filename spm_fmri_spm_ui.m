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
% (see also spm_spm_ui)
%
%__________________________________________________________________________
%
% spm_fmri_spm_ui configures the design matrix, data specification and
% filtering that specify the ensuing statistical analysis. These
% arguments are passed to spm_spm that then performs the actual parameter
% estimation.
%
% The design matrix defines the experimental design and the nature of
% hypothesis testing to be implemented.  The design matrix has one row
% for each scan and one column for each effect or explanatory variable.
% (e.g. regressor or stimulus function).  The parameters are estimated in
% a least squares sense using the general linear model.  Specific profiles
% within these parameters are tested using a linear compound or contrast
% with the T or F statistic.  The resulting statistical map constitutes
% an SPM.  The SPM{T}/{F} is then characterized in terms of focal or regional
% differences by assuming that (under the null hypothesis) the components of
% the SPM (i.e. residual fields) behave as smooth stationary Gaussian fields.
%
% spm_fmri_spm_ui allows you to (i) specify a statistical model in terms
% of a design matrix, (ii) associate some data with a pre-specified design
% [or (iii) specify both the data and design] and then proceed to estimate
% the parameters of the model.
% Inferences can be made about the ensuing parameter estimates (at a first
% or fixed-effect level) in the results section, or they can be re-entered
% into a second (random-effect) level analysis by treating the session or
% subject-specific [contrasts of] parameter estimates as new summary data.
% Inferences at any level obtain by specifying appropriate T or F contrasts
% in the results section to produce SPMs and tables of p values and statistics.
%
% spm_fmri_spm calls spm_fMRI_design which allows you to configure a
% design matrix in terms of events or epochs.
%
% spm_fMRI_design allows you to build design matrices with separable
% session-specific partitions.  Each partition may be the same (in which
% case it is only necessary to specify it once) or different.  Responses
% can be either event- or epoch related, The only distinction is the duration
% of the underlying input or stimulus function. Mathematically they are both
% modelled by convolving a series of delta (stick) or box functions (u),
% indicating the onset of an event or epoch with a set of basis
% functions.  These basis functions model the hemodynamic convolution,
% applied by the brain, to the inputs.  This convolution can be first-order
% or a generalized convolution modelled to second order (if you specify the
% Volterra option). [The same inputs are used by the hemodynamic model or
% or dynamic causal models which model the convolution explicitly in terms of
% hidden state variables (see spm_hdm_ui and spm_dcm_ui).]
% Basis functions can be used to plot estimated responses to single events
% once the parameters (i.e. basis function coefficients) have
% been estimated.  The importance of basis functions is that they provide
% a graceful transition between simple fixed response models (like the
% box-car) and finite impulse response (FIR) models, where there is one
% basis function for each scan following an event or epoch onset.  The
% nice thing about basis functions, compared to FIR models, is that data
% sampling and stimulus presentation does not have to be synchronized
% thereby allowing a uniform and unbiased sampling of peri-stimulus time.
%
% Event-related designs may be stochastic or deterministic.  Stochastic
% designs involve one of a number of trial-types occurring with a
% specified probably at successive intervals in time.  These
% probabilities can be fixed (stationary designs) or time-dependent
% (modulated or non-stationary designs).  The most efficient designs
% obtain when the probabilities of every trial type are equal.
% A critical issue in stochastic designs is whether to include null events
% If you wish to estimate the evoke response to a specific event
% type (as opposed to differential responses) then a null event must be
% included (even if it is not modelled explicitly).
%
% The choice of basis functions depends upon the nature of the inference
% sought.  One important consideration is whether you want to make
% inferences about compounds of parameters (i.e.  contrasts).  This is
% the case if (i) you wish to use a SPM{T} to look separately at
% activations and deactivations or (ii) you with to proceed to a second
% (random-effect) level of analysis.  If this is the case then (for
% event-related studies) use a canonical hemodynamic response function
% (HRF) and derivatives with respect to latency (and dispersion).  Unlike
% other bases, contrasts of these effects have a physical interpretation
% and represent a parsimonious way of characterising event-related
% responses.  Bases such as a Fourier set require the SPM{F} for
% inference.
%
% See spm_fMRI_design for more details about how designs are specified.
%
% Serial correlations in fast fMRI time-series are dealt with as
% described in spm_spm.  At this stage you need to specify the filtering
% that will be applied to the data (and design matrix) to give a
% generalized least squares (GLS) estimate of the parameters required.
% This filtering is important to ensure that the GLS estimate is
% efficient and that the error variance is estimated in an unbiased way.
%
% The serial correlations will be estimated with a ReML (restricted
% maximum likelihood) algorithm using an autoregressive AR(1) model
% during parameter estimation.  This estimate assumes the same
% correlation structure for each voxel, within each session.  The ReML
% estimates are then used to correct for non-sphericity during inference
% by adjusting the statistics and degrees of freedom appropriately.  The
% discrepancy between estimated and actual intrinsic (i.e. prior to
% filtering) correlations are greatest at low frequencies.  Therefore
% specification of the high-pass filter is particularly important.
%
% High-pass filtering is implemented at the level of the
% filtering matrix K (as opposed to entering as confounds in the design
% matrix).  The default cut-off period is 128 seconds.  Use 'explore design'
% to ensure this cut-off is not removing too much experimental variance.
% Note that high-pass filtering uses a residual forming matrix (i.e.
% it is not a convolution) and is simply to a way to remove confounds
% without estimating their parameters explicitly.  The constant term
% is also incorporated into this filter matrix.
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
% Map. 0:00-00
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Jean-Baptiste Poline & Christian Buchel
% $Id: spm_fmri_spm_ui.m 4178 2011-01-27 15:12:53Z guillaume $

SVNid  = '$Rev: 4178 $';

%-GUI setup
%--------------------------------------------------------------------------
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','fMRI stats model setup',0);

% get design matrix and/or data
%==========================================================================
if ~nargin

    str = 'specify design or data';
    if spm_input(str,1,'b',{'design','data'},[1 0]);

        % specify a design
        %------------------------------------------------------------------
        if sf_abort, spm_clf(Finter), return, end
        SPM   = spm_fMRI_design;
        spm_fMRI_design_show(SPM);
        return

    else

        % get design
        %------------------------------------------------------------------
        load(spm_select(1,'^SPM\.mat$','Select SPM.mat'));

    end

else

    % get design matrix
    %----------------------------------------------------------------------
    SPM       = spm_fMRI_design(SPM);

end

% get Repeat time
%--------------------------------------------------------------------------
try
    SPM.xY.RT;
catch
    SPM.xY.RT = spm_input('Interscan interval {secs}','+1');
end

% session and scan number
%--------------------------------------------------------------------------
nscan         = SPM.nscan;
nsess         = length(nscan);

% check data are specified
%--------------------------------------------------------------------------
try
    SPM.xY.P;
catch

    % get filenames
    %----------------------------------------------------------------------
    P         = [];
    for i = 1:nsess
        str   = sprintf('select scans for session %0.0f',i);
        q     = spm_select(nscan(i),'image',str);
        P     = strvcat(P,q);
    end

    % place in data field
    %----------------------------------------------------------------------
    SPM.xY.P  = P;

end


% Assemble remaining design parameters
%==========================================================================
SPM.SPMid     = spm('FnBanner',mfilename,SVNid);


% Global normalization
%--------------------------------------------------------------------------
try
    SPM.xGX.iGXcalc;
catch
    spm_input('Global intensity normalisation...',1,'d',mfilename)
    str             = 'remove Global effects';
    SPM.xGX.iGXcalc = spm_input(str,'+1','scale|none',{'Scaling' 'None'});
end
SPM.xGX.sGXcalc = 'mean voxel value';
SPM.xGX.sGMsca  = 'session specific';


% High-pass filtering and serial correlations
%==========================================================================

% low frequency confounds
%--------------------------------------------------------------------------
try
    myLastWarn = 0;
    HParam     = [SPM.xX.K(:).HParam];
    if length(HParam) == 1 
        HParam = HParam*ones(1,nsess);
    elseif length(HParam) ~= nsess
        myLastWarn = 1;
        error('Continue with manual HPF specification in the catch block');
    end
catch
    % specify low frequency confounds
    %----------------------------------------------------------------------
    spm_input('Temporal autocorrelation options','+1','d',mfilename)
    switch spm_input('High-pass filter?','+1','b','none|specify');

        case 'specify'  % default in seconds
            %--------------------------------------------------------------
            HParam = spm_get_defaults('stats.fmri.hpf')*ones(1,nsess);
            str    = 'cutoff period (secs)';
            HParam = spm_input(str,'+1','e',HParam,[1 nsess]);

        case 'none'     % Inf seconds (i.e. constant term only)
            %--------------------------------------------------------------
            HParam = Inf(1,nsess);
    end
    if myLastWarn
        warning('SPM:InvalidHighPassFilterSpec',...
            ['Different number of High-pass filter values and sessions.\n',...
            'HPF filter configured manually. Design setup will proceed.']);
        clear myLastWarn
    end
end

% create and set filter struct
%--------------------------------------------------------------------------
for  i = 1:nsess
    K(i) = struct('HParam', HParam(i),...
                  'row',    SPM.Sess(i).row,...
                  'RT',     SPM.xY.RT);
end
SPM.xX.K = spm_filter(K);

% intrinsic autocorrelations (Vi)
%--------------------------------------------------------------------------
try
    cVi   = SPM.xVi.form;
catch
    % Construct Vi structure for non-sphericity ReML estimation
    %----------------------------------------------------------------------
    str   = 'Correct for serial correlations?';
    cVi   = {'none','AR(1)'};
    cVi   = spm_input(str,'+1','b',cVi);
end

% create Vi struct
%--------------------------------------------------------------------------

if ~ischar(cVi) % AR coefficient specified
    %----------------------------------------------------------------------
    SPM.xVi.Vi = spm_Ce(nscan,cVi(1));
    cVi        = ['AR( ' sprintf('%0.1f ',cVi) ')'];

else
    switch lower(cVi)

        case 'none'     %  xVi.V is i.i.d
            %--------------------------------------------------------------
            SPM.xVi.V  = speye(sum(nscan));
            cVi        = 'i.i.d';


        otherwise       % otherwise assume AR(0.2) in xVi.Vi
            %--------------------------------------------------------------
            SPM.xVi.Vi = spm_Ce(nscan,0.2);
            cVi        = 'AR(0.2)';

    end
end

SPM.xVi.form = cVi;


%==========================================================================
% - C O N F I G U R E   D E S I G N
%==========================================================================
spm_clf(Finter);
spm('FigName','Configuring, please wait...',Finter,CmdLine);
spm('Pointer','Watch');


% get file identifiers
%==========================================================================

%-Map files
%--------------------------------------------------------------------------
fprintf('%-40s: ','Mapping files')                                      %-#
VY    = spm_vol(SPM.xY.P);
fprintf('%30s\n','...done')                                             %-#

%-check internal consistency of images
%--------------------------------------------------------------------------
spm_check_orientations(VY);

%-place in xY
%--------------------------------------------------------------------------
SPM.xY.VY = VY;


%-Compute Global variate
%==========================================================================
GM    = 100;
q     = length(VY);
g     = zeros(q,1);
fprintf('%-40s: %30s','Calculating globals',' ')                        %-#
for i = 1:q
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),sprintf('%4d/%-4d',i,q))%-#
    g(i) = spm_global(VY(i));
end
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')                %-#

% scale if specified (otherwise session specific grand mean scaling)
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
    SPM.xY.VY(i).pinfo(1:2,:) = SPM.xY.VY(i).pinfo(1:2,:)*gSF(i);
end

%-place global variates in global structure
%--------------------------------------------------------------------------
SPM.xGX.rg  = g;
SPM.xGX.GM  = GM;
SPM.xGX.gSF = gSF;


%-Masking structure automatically set to 80% of mean
%==========================================================================
try
    TH = g.*gSF*spm_get_defaults('mask.thresh');
catch
    TH = g.*gSF*0.8;
end
SPM.xM = struct('T',    ones(q,1),...
                'TH',   TH,...
                'I',    0,...
                'VM',   {[]},...
                'xs',   struct('Masking','analysis threshold'));


%-Design description - for saving and display
%==========================================================================
for i     = 1:nsess, ntr(i) = length(SPM.Sess(i).U); end
Fstr      = sprintf('[min] Cutoff: %d {s}',min([SPM.xX.K(:).HParam]));
SPM.xsDes = struct(...
    'Basis_functions',      SPM.xBF.name,...
    'Number_of_sessions',   sprintf('%d',nsess),...
    'Trials_per_session',   sprintf('%-3d',ntr),...
    'Interscan_interval',   sprintf('%0.2f {s}',SPM.xY.RT),...
    'High_pass_Filter',     sprintf('Cutoff: %d {s}',SPM.xX.K(1).HParam),...
    'Global_calculation',   SPM.xGX.sGXcalc,...
    'Grand_mean_scaling',   SPM.xGX.sGMsca,...
    'Global_normalisation', SPM.xGX.iGXcalc);


%-Save SPM.mat
%==========================================================================
fprintf('%-40s: ','Saving SPM configuration')                           %-#
if spm_matlab_version_chk('7') >= 0,
    save('SPM.mat', 'SPM', '-V6');
else
    save('SPM.mat', 'SPM');
end
fprintf('%30s\n','...SPM.mat saved')                                    %-#


%-Display Design report
%==========================================================================
if ~CmdLine
    fprintf('%-40s: ','Design reporting')                               %-#
    fname = cat(1,{SPM.xY.VY.fname}');
    spm_DesRep('DesMtx',SPM.xX,fname,SPM.xsDes)
    fprintf('%30s\n','...done')                                         %-#
end


%-End: Cleanup GUI
%==========================================================================
spm_clf(Finter)
spm('FigName','Stats: configured',Finter,CmdLine);
spm('Pointer','Arrow')


%==========================================================================
%- S U B - F U N C T I O N S
%==========================================================================

function abort = sf_abort
%==========================================================================
if exist(fullfile(pwd,'SPM.mat'),'file')
    str = { 'Current directory contains existing SPM file:',...
            'Continuing will overwrite existing file!'};

    abort = spm_input(str,1,'bd','stop|continue',[1,0],1,mfilename);
    if abort, fprintf('%-40s: %30s\n\n',...
            'Abort...   (existing SPM files)',spm('time')), end
else
    abort = 0;
end
