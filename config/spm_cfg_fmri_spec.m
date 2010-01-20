function fmri_spec = spm_cfg_fmri_spec
% SPM Configuration file for fMRI model specification
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_fmri_spec.m 3691 2010-01-20 17:08:30Z guillaume $

rev = '$Rev: 3691 $';
% ---------------------------------------------------------------------
% dir Directory
% ---------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.help    = {'Select a directory where the SPM.mat file containing the specified design matrix will be written.'};
dir.filter = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];
% ---------------------------------------------------------------------
% units Units for design
% ---------------------------------------------------------------------
units         = cfg_menu;
units.tag     = 'units';
units.name    = 'Units for design';
units.help    = {'The onsets of events or blocks can be specified in either scans or seconds.'};
units.labels = {
                'Scans'
                'Seconds'
}';
units.values = {
                'scans'
                'secs'
}';
% ---------------------------------------------------------------------
% RT Interscan interval
% ---------------------------------------------------------------------
RT         = cfg_entry;
RT.tag     = 'RT';
RT.name    = 'Interscan interval';
RT.help    = {'Interscan interval, TR, (specified in seconds).  This is the time between acquiring a plane of one volume and the same plane in the next volume.  It is assumed to be constant throughout.'};
RT.strtype = 'e';
RT.num     = [1 1];
% ---------------------------------------------------------------------
% fmri_t Microtime resolution
% ---------------------------------------------------------------------
fmri_t         = cfg_entry;
fmri_t.tag     = 'fmri_t';
fmri_t.name    = 'Microtime resolution';
fmri_t.help    = {
                  'The microtime resolution, t, is the number of time-bins per scan used when building regressors. '
                  ''
                  'Do not change this parameter unless you have a long TR and wish to shift regressors so that they are aligned to a particular slice. '
}';
fmri_t.strtype = 'e';
fmri_t.num     = [1 1];
fmri_t.def     = @(val)spm_get_defaults('stats.fmri.fmri_t', val{:});
% ---------------------------------------------------------------------
% fmri_t0 Microtime onset
% ---------------------------------------------------------------------
fmri_t0         = cfg_entry;
fmri_t0.tag     = 'fmri_t0';
fmri_t0.name    = 'Microtime onset';
fmri_t0.help    = {
                   'The microtime onset, t0, is the first time-bin at which the regressors are resampled to coincide with data acquisition.  If t0 = 1 then the regressors will be appropriate for the first slice.  If you want to temporally realign the regressors so that they match responses in the middle slice then make t0 = t/2 (assuming there is a negligible gap between volume acquisitions). '
                   ''
                   'Do not change the default setting unless you have a long TR. '
}';
fmri_t0.strtype = 'e';
fmri_t0.num     = [1 1];
fmri_t0.def     = @(val)spm_get_defaults('stats.fmri.fmri_t0', val{:});
% ---------------------------------------------------------------------
% timing Timing parameters
% ---------------------------------------------------------------------
timing         = cfg_branch;
timing.tag     = 'timing';
timing.name    = 'Timing parameters';
timing.val     = {units RT fmri_t fmri_t0 };
timing.help    = {
                  'Specify various timing parameters needed to construct the design matrix. This includes the units of the design specification and the interscan interval.'
                  ''
                  'Also, with longs TRs you may want to shift the regressors so that they are aligned to a particular slice.  This is effected by changing the microtime resolution and onset. '
}';
% ---------------------------------------------------------------------
% scans Scans
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {'Select the fMRI scans for this session.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [1 Inf];
% ---------------------------------------------------------------------
% name Name
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Name';
name.help    = {'Condition Name'};
name.strtype = 's';
name.num     = [1 Inf];
% ---------------------------------------------------------------------
% onset Onsets
% ---------------------------------------------------------------------
onset         = cfg_entry;
onset.tag     = 'onset';
onset.name    = 'Onsets';
onset.help    = {'Specify a vector of onset times for this condition type. '};
onset.strtype = 'e';
onset.num     = [Inf 1];
% ---------------------------------------------------------------------
% duration Durations
% ---------------------------------------------------------------------
duration         = cfg_entry;
duration.tag     = 'duration';
duration.name    = 'Durations';
duration.help    = {'Specify the event durations. Epoch and event-related responses are modeled in exactly the same way but by specifying their different durations.  Events are specified with a duration of 0.  If you enter a single number for the durations it will be assumed that all trials conform to this duration. If you have multiple different durations, then the number must match the number of onset times.'};
duration.strtype = 'e';
duration.num     = [Inf 1];
% ---------------------------------------------------------------------
% tmod Time Modulation
% ---------------------------------------------------------------------
tmod         = cfg_menu;
tmod.tag     = 'tmod';
tmod.name    = 'Time Modulation';
tmod.help    = {
                'This option allows for the characterisation of linear or nonlinear time effects. For example, 1st order modulation would model the stick functions and a linear change of the stick function heights over time. Higher order modulation will introduce further columns that contain the stick functions scaled by time squared, time cubed etc.'
                ''
                'Interactions or response modulations can enter at two levels.  Firstly the stick function itself can be modulated by some parametric variate (this can be time or some trial-specific variate like reaction time) modeling the interaction between the trial and the variate or, secondly interactions among the trials themselves can be modeled using a Volterra series formulation that accommodates interactions over time (and therefore within and between trial types).'
}';
tmod.labels = {
               'No Time Modulation'
               '1st order Time Modulation'
               '2nd order Time Modulation'
               '3rd order Time Modulation'
               '4th order Time Modulation'
               '5th order Time Modulation'
               '6th order Time Modulation'
}';
tmod.values = {0 1 2 3 4 5 6};
tmod.val    = {0};
% ---------------------------------------------------------------------
% name Name
% ---------------------------------------------------------------------
name1         = cfg_entry;
name1.tag     = 'name';
name1.name    = 'Name';
name1.help    = {'Enter a name for this parameter.'};
name1.strtype = 's';
name1.num     = [1 Inf];
% ---------------------------------------------------------------------
% param Values
% ---------------------------------------------------------------------
param         = cfg_entry;
param.tag     = 'param';
param.name    = 'Values';
param.help    = {'Enter a vector of values, one for each occurence of the event.'};
param.strtype = 'e';
param.num     = [Inf 1];
% ---------------------------------------------------------------------
% poly Polynomial Expansion
% ---------------------------------------------------------------------
poly         = cfg_menu;
poly.tag     = 'poly';
poly.name    = 'Polynomial Expansion';
poly.help    = {'For example, 1st order modulation would model the stick functions and a linear change of the stick function heights over different values of the parameter. Higher order modulation will introduce further columns that contain the stick functions scaled by parameter squared, cubed etc.'};
poly.labels = {
               '1st order'
               '2nd order'
               '3rd order'
               '4th order'
               '5th order'
               '6th order'
}';
poly.values = {1 2 3 4 5 6};
% ---------------------------------------------------------------------
% pmod Parameter
% ---------------------------------------------------------------------
pmod         = cfg_branch;
pmod.tag     = 'pmod';
pmod.name    = 'Parameter';
pmod.val     = {name1 param poly };
pmod.help    = {
                'Model interractions with user specified parameters. This allows nonlinear effects relating to some other measure to be modelled in the design matrix.'
                ''
                'Interactions or response modulations can enter at two levels.  Firstly the stick function itself can be modulated by some parametric variate (this can be time or some trial-specific variate like reaction time) modeling the interaction between the trial and the variate or, secondly interactions among the trials themselves can be modeled using a Volterra series formulation that accommodates interactions over time (and therefore within and between trial types).'
}';
% ---------------------------------------------------------------------
% generic Parametric Modulations
% ---------------------------------------------------------------------
generic2         = cfg_repeat;
generic2.tag     = 'generic';
generic2.name    = 'Parametric Modulations';
generic2.help    = {'The stick function itself can be modulated by some parametric variate (this can be time or some trial-specific variate like reaction time) modeling the interaction between the trial and the variate. The events can be modulated by zero or more parameters.'};
generic2.values  = {pmod };
generic2.num     = [0 Inf];
% ---------------------------------------------------------------------
% cond Condition
% ---------------------------------------------------------------------
cond         = cfg_branch;
cond.tag     = 'cond';
cond.name    = 'Condition';
cond.val     = {name onset duration tmod generic2 };
cond.check   = @cond_check;
cond.help    = {'An array of input functions is contructed, specifying occurrence events or epochs (or both). These are convolved with a basis set at a later stage to give regressors that enter into the design matrix. Interactions of evoked responses with some parameter (time or a specified variate) enter at this stage as additional columns in the design matrix with each trial multiplied by the [expansion of the] trial-specific parameter. The 0th order expansion is simply the main effect in the first column.'};
% ---------------------------------------------------------------------
% generic Conditions
% ---------------------------------------------------------------------
generic1         = cfg_repeat;
generic1.tag     = 'generic';
generic1.name    = 'Conditions';
generic1.help    = {'You are allowed to combine both event- and epoch-related responses in the same model and/or regressor. Any number of condition (event or epoch) types can be specified.  Epoch and event-related responses are modeled in exactly the same way by specifying their onsets [in terms of onset times] and their durations.  Events are specified with a duration of 0.  If you enter a single number for the durations it will be assumed that all trials conform to this duration.For factorial designs, one can later associate these experimental conditions with the appropriate levels of experimental factors. '};
generic1.values  = {cond };
generic1.num     = [0 Inf];
% ---------------------------------------------------------------------
% multi Multiple conditions
% ---------------------------------------------------------------------
multi         = cfg_files;
multi.tag     = 'multi';
multi.name    = 'Multiple conditions';
multi.val{1} = {''};
multi.help    = {
                 'Select the *.mat file containing details of your multiple experimental conditions. '
                 ''
                 'If you have multiple conditions then entering the details a condition at a time is very inefficient. This option can be used to load all the required information in one go. You will first need to create a *.mat file containing the relevant information. '
                 ''
                 'This *.mat file must include the following cell arrays (each 1 x n): names, onsets and durations. eg. names=cell(1,5), onsets=cell(1,5), durations=cell(1,5), then names{2}=''SSent-DSpeak'', onsets{2}=[3 5 19 222], durations{2}=[0 0 0 0], contain the required details of the second condition. These cell arrays may be made available by your stimulus delivery program, eg. COGENT. The duration vectors can contain a single entry if the durations are identical for all events.'
                 ''
                 'Time and Parametric effects can also be included. For time modulation include a cell array (1 x n) called tmod. It should have a have a single number in each cell. Unused cells may contain either a 0 or be left empty. The number specifies the order of time modulation from 0 = No Time Modulation to 6 = 6th Order Time Modulation. eg. tmod{3} = 1, modulates the 3rd condition by a linear time effect.'
                 ''
                 'For parametric modulation include a structure array, which is up to 1 x n in size, called pmod. n must be less than or equal to the number of cells in the names/onsets/durations cell arrays. The structure array pmod must have the fields: name, param and poly.  Each of these fields is in turn a cell array to allow the inclusion of one or more parametric effects per column of the design. The field name must be a cell array containing strings. The field param is a cell array containing a vector of parameters. Remember each parameter must be the same length as its corresponding onsets vector. The field poly is a cell array (for consistency) with each cell containing a single number specifying the order of the polynomial expansion from 1 to 6.'
                 ''
                 'Note that each condition is assigned its corresponding entry in the structure array (condition 1 parametric modulators are in pmod(1), condition 2 parametric modulators are in pmod(2), etc. Within a condition multiple parametric modulators are accessed via each fields cell arrays. So for condition 1, parametric modulator 1 would be defined in  pmod(1).name{1}, pmod(1).param{1}, and pmod(1).poly{1}. A second parametric modulator for condition 1 would be defined as pmod(1).name{2}, pmod(1).param{2} and pmod(1).poly{2}. If there was also a parametric modulator for condition 2, then remember the first modulator for that condition is in cell array 1: pmod(2).name{1}, pmod(2).param{1}, and pmod(2).poly{1}. If some, but not all conditions are parametrically modulated, then the non-modulated indices in the pmod structure can be left blank. For example, if conditions 1 and 3 but not condition 2 are modulated, then specify pmod(1) and pmod(3). Similarly, if conditions 1 and 2 are modulated but there are 3 conditions overall, it is only necessary for pmod to be a 1 x 2 structure array.'
                 ''
                 'EXAMPLE:'
                 'Make an empty pmod structure: '
                 '  pmod = struct(''name'',{''''},''param'',{},''poly'',{});'
                 'Specify one parametric regressor for the first condition: '
                 '  pmod(1).name{1}  = ''regressor1'';'
                 '  pmod(1).param{1} = [1 2 4 5 6];'
                 '  pmod(1).poly{1}  = 1;'
                 'Specify 2 parametric regressors for the second condition: '
                 '  pmod(2).name{1}  = ''regressor2-1'';'
                 '  pmod(2).param{1} = [1 3 5 7]; '
                 '  pmod(2).poly{1}  = 1;'
                 '  pmod(2).name{2}  = ''regressor2-2'';'
                 '  pmod(2).param{2} = [2 4 6 8 10];'
                 '  pmod(2).poly{2}  = 1;'
                 ''
                 'The parametric modulator should be mean corrected if appropriate. Unused structure entries should have all fields left empty.'
}';
multi.filter = 'mat';
multi.ufilter = '.*';
multi.num     = [0 1];
% ---------------------------------------------------------------------
% name Name
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Name';
name.help    = {'Enter name of regressor eg. First movement parameter'};
name.strtype = 's';
name.num     = [1 Inf];
% ---------------------------------------------------------------------
% val Value
% ---------------------------------------------------------------------
val         = cfg_entry;
val.tag     = 'val';
val.name    = 'Value';
val.help    = {'Enter the vector of regressor values'};
val.strtype = 'e';
val.num     = [Inf 1];
% ---------------------------------------------------------------------
% regress Regressor
% ---------------------------------------------------------------------
regress         = cfg_branch;
regress.tag     = 'regress';
regress.name    = 'Regressor';
regress.val     = {name val };
regress.help    = {'regressor'};
% ---------------------------------------------------------------------
% generic Regressors
% ---------------------------------------------------------------------
generic2         = cfg_repeat;
generic2.tag     = 'generic';
generic2.name    = 'Regressors';
generic2.help    = {'Regressors are additional columns included in the design matrix, which may model effects that would not be convolved with the haemodynamic response.  One such example would be the estimated movement parameters, which may confound the data.'};
generic2.values  = {regress };
generic2.num     = [0 Inf];
% ---------------------------------------------------------------------
% multi_reg Multiple regressors
% ---------------------------------------------------------------------
multi_reg         = cfg_files;
multi_reg.tag     = 'multi_reg';
multi_reg.name    = 'Multiple regressors';
multi_reg.val{1} = {''};
multi_reg.help    = {
                     'Select the *.mat/*.txt file containing details of your multiple regressors. '
                     ''
                     'If you have multiple regressors eg. realignment parameters, then entering the details a regressor at a time is very inefficient. This option can be used to load all the required information in one go. '
                     ''
                     'You will first need to create a *.mat file containing a matrix R or a *.txt file containing the regressors. Each column of R will contain a different regressor. When SPM creates the design matrix the regressors will be named R1, R2, R3, ..etc.'
}';
multi_reg.filter = 'mat';
multi_reg.ufilter = '.*';
multi_reg.num     = [0 1];
% ---------------------------------------------------------------------
% hpf High-pass filter
% ---------------------------------------------------------------------
hpf         = cfg_entry;
hpf.tag     = 'hpf';
hpf.name    = 'High-pass filter';
hpf.help    = {'The default high-pass filter cutoff is 128 seconds.Slow signal drifts with a period longer than this will be removed. Use ''explore design'' to ensure this cut-off is not removing too much experimental variance. High-pass filtering is implemented using a residual forming matrix (i.e. it is not a convolution) and is simply to a way to remove confounds without estimating their parameters explicitly.  The constant term is also incorporated into this filter matrix.'};
hpf.strtype = 'e';
hpf.num     = [1 1];
hpf.def     = @(val)spm_get_defaults('stats.fmri.hpf', val{:});
% ---------------------------------------------------------------------
% sess Subject/Session
% ---------------------------------------------------------------------
sess         = cfg_branch;
sess.tag     = 'sess';
sess.name    = 'Subject/Session';
sess.val     = {scans generic1 multi generic2 multi_reg hpf };
sess.check   = @sess_check;
sess.help    = {'The design matrix for fMRI data consists of one or more separable, session-specific partitions.  These partitions are usually either one per subject, or one per fMRI scanning session for that subject.'};
% ---------------------------------------------------------------------
% generic Data & Design
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Data & Design';
generic.help    = {
                   'The design matrix defines the experimental design and the nature of hypothesis testing to be implemented.  The design matrix has one row for each scan and one column for each effect or explanatory variable. (e.g. regressor or stimulus function).  '
                   ''
                   'This allows you to build design matrices with separable session-specific partitions.  Each partition may be the same (in which case it is only necessary to specify it once) or different.  Responses can be either event- or epoch related, where the latter model involves prolonged and possibly time-varying responses to state-related changes in experimental conditions.  Event-related response are modelled in terms of responses to instantaneous events.  Mathematically they are both modelled by convolving a series of delta (stick) or box-car functions, encoding the input or stimulus function. with a set of hemodynamic basis functions.'
}';
generic.values  = {sess };
generic.num     = [1 Inf];
% ---------------------------------------------------------------------
% name Name
% ---------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Name';
name.help    = {'Name of factor, eg. ''Repetition'' '};
name.strtype = 's';
name.num     = [1 Inf];
% ---------------------------------------------------------------------
% levels Levels
% ---------------------------------------------------------------------
levels         = cfg_entry;
levels.tag     = 'levels';
levels.name    = 'Levels';
levels.help    = {'Enter number of levels for this factor, eg. 2'};
levels.strtype = 'e';
levels.num     = [Inf 1];
% ---------------------------------------------------------------------
% fact Factor
% ---------------------------------------------------------------------
fact         = cfg_branch;
fact.tag     = 'fact';
fact.name    = 'Factor';
fact.val     = {name levels };
fact.help    = {'Add a new factor to your experimental design'};
% ---------------------------------------------------------------------
% generic Factorial design
% ---------------------------------------------------------------------
generic1         = cfg_repeat;
generic1.tag     = 'generic';
generic1.name    = 'Factorial design';
generic1.help    = {
                    'If you have a factorial design then SPM can automatically generate the contrasts necessary to test for the main effects and interactions. '
                    ''
                    'This includes the F-contrasts necessary to test for these effects at the within-subject level (first level) and the simple contrasts necessary to generate the contrast images for a between-subject (second-level) analysis.'
                    ''
                    'To use this option, create as many factors as you need and provide a name and number of levels for each.  SPM assumes that the condition numbers of the first factor change slowest, the second factor next slowest etc. It is best to write down the contingency table for your design to ensure this condition is met. This table relates the levels of each factor to the conditions. '
                    ''
                    'For example, if you have 2-by-3 design  your contingency table has two rows and three columns where the the first factor spans the rows, and the second factor the columns. The numbers of the conditions are 1,2,3 for the first row and 4,5,6 for the second. '
}';
generic1.values  = {fact };
generic1.num     = [0 Inf];
% ---------------------------------------------------------------------
% derivs Model derivatives
% ---------------------------------------------------------------------
derivs         = cfg_menu;
derivs.tag     = 'derivs';
derivs.name    = 'Model derivatives';
derivs.help    = {'Model HRF Derivatives. The canonical HRF combined with time and dispersion derivatives comprise an ''informed'' basis set, as the shape of the canonical response conforms to the hemodynamic response that is commonly observed. The incorporation of the derivate terms allow for variations in subject-to-subject and voxel-to-voxel responses. The time derivative allows the peak response to vary by plus or minus a second and the dispersion derivative allows the width of the response to vary. The informed basis set requires an SPM{F} for inference. T-contrasts over just the canonical are perfectly valid but assume constant delay/dispersion. The informed basis set compares favourably with eg. FIR bases on many data sets. '};
derivs.labels = {
                 'No derivatives'
                 'Time derivatives'
                 'Time and Dispersion derivatives'
}';
derivs.values = {[0 0] [1 0] [1 1]};
derivs.val    = {[0 0]};
% ---------------------------------------------------------------------
% hrf Canonical HRF
% ---------------------------------------------------------------------
hrf         = cfg_branch;
hrf.tag     = 'hrf';
hrf.name    = 'Canonical HRF';
hrf.val     = {derivs };
hrf.help    = {'Canonical Hemodynamic Response Function. This is the default option. Contrasts of these effects have a physical interpretation and represent a parsimonious way of characterising event-related responses. This option is also useful if you wish to look separately at activations and deactivations (this is implemented using a t-contrast with a +1 or -1 entry over the canonical regressor). '};
% ---------------------------------------------------------------------
% length Window length
% ---------------------------------------------------------------------
length         = cfg_entry;
length.tag     = 'length';
length.name    = 'Window length';
length.help    = {'Post-stimulus window length (in seconds)'};
length.strtype = 'e';
length.num     = [1 1];
% ---------------------------------------------------------------------
% order Order
% ---------------------------------------------------------------------
order         = cfg_entry;
order.tag     = 'order';
order.name    = 'Order';
order.help    = {'Number of basis functions'};
order.strtype = 'e';
order.num     = [1 1];
% ---------------------------------------------------------------------
% fourier Fourier Set
% ---------------------------------------------------------------------
fourier         = cfg_branch;
fourier.tag     = 'fourier';
fourier.name    = 'Fourier Set';
fourier.val     = {length order };
fourier.help    = {'Fourier basis functions. This option requires an SPM{F} for inference.'};
% ---------------------------------------------------------------------
% length Window length
% ---------------------------------------------------------------------
length         = cfg_entry;
length.tag     = 'length';
length.name    = 'Window length';
length.help    = {'Post-stimulus window length (in seconds)'};
length.strtype = 'e';
length.num     = [1 1];
% ---------------------------------------------------------------------
% order Order
% ---------------------------------------------------------------------
order         = cfg_entry;
order.tag     = 'order';
order.name    = 'Order';
order.help    = {'Number of basis functions'};
order.strtype = 'e';
order.num     = [1 1];
% ---------------------------------------------------------------------
% fourier_han Fourier Set (Hanning)
% ---------------------------------------------------------------------
fourier_han         = cfg_branch;
fourier_han.tag     = 'fourier_han';
fourier_han.name    = 'Fourier Set (Hanning)';
fourier_han.val     = {length order };
fourier_han.help    = {'Fourier basis functions with Hanning Window - requires SPM{F} for inference.'};
% ---------------------------------------------------------------------
% length Window length
% ---------------------------------------------------------------------
length         = cfg_entry;
length.tag     = 'length';
length.name    = 'Window length';
length.help    = {'Post-stimulus window length (in seconds)'};
length.strtype = 'e';
length.num     = [1 1];
% ---------------------------------------------------------------------
% order Order
% ---------------------------------------------------------------------
order         = cfg_entry;
order.tag     = 'order';
order.name    = 'Order';
order.help    = {'Number of basis functions'};
order.strtype = 'e';
order.num     = [1 1];
% ---------------------------------------------------------------------
% gamma Gamma Functions
% ---------------------------------------------------------------------
gamma         = cfg_branch;
gamma.tag     = 'gamma';
gamma.name    = 'Gamma Functions';
gamma.val     = {length order };
gamma.help    = {'Gamma basis functions - requires SPM{F} for inference.'};
% ---------------------------------------------------------------------
% length Window length
% ---------------------------------------------------------------------
length         = cfg_entry;
length.tag     = 'length';
length.name    = 'Window length';
length.help    = {'Post-stimulus window length (in seconds)'};
length.strtype = 'e';
length.num     = [1 1];
% ---------------------------------------------------------------------
% order Order
% ---------------------------------------------------------------------
order         = cfg_entry;
order.tag     = 'order';
order.name    = 'Order';
order.help    = {'Number of basis functions'};
order.strtype = 'e';
order.num     = [1 1];
% ---------------------------------------------------------------------
% fir Finite Impulse Response
% ---------------------------------------------------------------------
fir         = cfg_branch;
fir.tag     = 'fir';
fir.name    = 'Finite Impulse Response';
fir.val     = {length order };
fir.help    = {'Finite impulse response - requires SPM{F} for inference.'};
% ---------------------------------------------------------------------
% bases Basis Functions
% ---------------------------------------------------------------------
bases         = cfg_choice;
bases.tag     = 'bases';
bases.name    = 'Basis Functions';
bases.val     = {hrf };
bases.help    = {'The most common choice of basis function is the Canonical HRF with or without time and dispersion derivatives. '};
bases.values  = {hrf fourier fourier_han gamma fir };
% ---------------------------------------------------------------------
% volt Model Interactions (Volterra)
% ---------------------------------------------------------------------
volt         = cfg_menu;
volt.tag     = 'volt';
volt.name    = 'Model Interactions (Volterra)';
volt.help    = {
                'Generalized convolution of inputs (U) with basis set (bf).'
                ''
                'For first order expansions the causes are simply convolved (e.g. stick functions) in U.u by the basis functions in bf to create a design matrix X.  For second order expansions new entries appear in ind, bf and name that correspond to the interaction among the orginal causes. The basis functions for these efects are two dimensional and are used to assemble the second order kernel. Second order effects are computed for only the first column of U.u.'
                'Interactions or response modulations can enter at two levels.  Firstly the stick function itself can be modulated by some parametric variate (this can be time or some trial-specific variate like reaction time) modeling the interaction between the trial and the variate or, secondly interactions among the trials themselves can be modeled using a Volterra series formulation that accommodates interactions over time (and therefore within and between trial types).'
}';
volt.labels = {
               'Do not model Interactions'
               'Model Interactions'
}';
volt.values = {1 2};
volt.val    = {1};
% ---------------------------------------------------------------------
% global Global normalisation
% ---------------------------------------------------------------------
xGlobal         = cfg_menu;
xGlobal.tag     = 'global';
xGlobal.name    = 'Global normalisation';
xGlobal.help    = {'Global intensity normalisation'};
xGlobal.labels = {
                  'Scaling'
                  'None'
}';
xGlobal.values = {
                  'Scaling'
                  'None'
}';
xGlobal.val     = {'None'};
% ---------------------------------------------------------------------
% mask Explicit mask
% ---------------------------------------------------------------------
mask         = cfg_files;
mask.tag     = 'mask';
mask.name    = 'Explicit mask';
mask.val{1} = {''};
mask.help    = {'Specify an image for explicitly masking the analysis. A sensible option here is to use a segmention of structural images to specify a within-brain mask. If you select that image as an explicit mask then only those voxels in the brain will be analysed. This both speeds the estimation and restricts SPMs/PPMs to within-brain voxels. Alternatively, if such structural images are unavailble or no masking is required, then leave this field empty.'};
mask.filter = 'image';
mask.ufilter = '.*';
mask.num     = [0 1];
% ---------------------------------------------------------------------
% cvi Serial correlations
% ---------------------------------------------------------------------
cvi         = cfg_menu;
cvi.tag     = 'cvi';
cvi.name    = 'Serial correlations';
cvi.help    = {
               'Serial correlations in fMRI time series due to aliased biorhythms and unmodelled neuronal activity can be accounted for using an autoregressive AR(1) model during Classical (ReML) parameter estimation.  '
               ''
               'This estimate assumes the same correlation structure for each voxel, within each session.  ReML estimates are then used to correct for non-sphericity during inference by adjusting the statistics and degrees of freedom appropriately.  The discrepancy between estimated and actual intrinsic (i.e. prior to filtering) correlations are greatest at low frequencies.  Therefore specification of the high-pass filter is particularly important. '
               ''
               'Serial correlation can be ignored if you choose the ''none'' option. Note that the above options only apply if you later specify that your model will be estimated using the Classical (ReML) approach. If you choose Bayesian estimation these options will be ignored. For Bayesian estimation, the choice of noisemodel (AR model order) is made under the estimation options. '
}';
cvi.labels = {
              'none'
              'AR(1)'
}';
cvi.values = {
              'none'
              'AR(1)'
}';
cvi.def     = @(val)spm_get_defaults('stats.fmri.cvi', val{:});
% ---------------------------------------------------------------------
% fmri_spec fMRI model specification
% ---------------------------------------------------------------------
fmri_spec         = cfg_exbranch;
fmri_spec.tag     = 'fmri_spec';
fmri_spec.name    = 'fMRI model specification';
fmri_spec.val     = {dir timing generic generic1 bases volt xGlobal mask cvi };
fmri_spec.help    = {
                     'Statistical analysis of fMRI data uses a mass-univariate approach based on General Linear Models (GLMs). It comprises the following steps (1) specification of the GLM design matrix, fMRI data files and filtering (2) estimation of GLM paramaters using classical or Bayesian approaches and (3) interrogation of results using contrast vectors to produce Statistical Parametric Maps (SPMs) or Posterior Probability Maps (PPMs).'
                     ''
                     'The design matrix defines the experimental design and the nature of hypothesis testing to be implemented.  The design matrix has one row for each scan and one column for each effect or explanatory variable. (eg. regressor or stimulus function). You can build design matrices with separable session-specific partitions.  Each partition may be the same (in which case it is only necessary to specify it once) or different. '
                     ''
                     'Responses can be either event- or epoch related, the only distinction is the duration of the underlying input or stimulus function. Mathematically they are both modeled by convolving a series of delta (stick) or box functions (u), indicating the onset of an event or epoch with a set of basis functions.  These basis functions model the hemodynamic convolution, applied by the brain, to the inputs.  This convolution can be first-order or a generalized convolution modeled to second order (if you specify the Volterra option). The same inputs are used by the Hemodynamic model or Dynamic Causal Models which model the convolution explicitly in terms of hidden state variables. '
                     ''
                     'Basis functions can be used to plot estimated responses to single events once the parameters (i.e. basis function coefficients) have been estimated.  The importance of basis functions is that they provide a graceful transition between simple fixed response models (like the box-car) and finite impulse response (FIR) models, where there is one basis function for each scan following an event or epoch onset.  The nice thing about basis functions, compared to FIR models, is that data sampling and stimulus presentation does not have to be synchronized thereby allowing a uniform and unbiased sampling of peri-stimulus time.'
                     ''
                     'Event-related designs may be stochastic or deterministic.  Stochastic designs involve one of a number of trial-types occurring with a specified probability at successive intervals in time.  These probabilities can be fixed (stationary designs) or time-dependent (modulated or non-stationary designs).  The most efficient designs obtain when the probabilities of every trial type are equal. A critical issue in stochastic designs is whether to include null events If you wish to estimate the evoked response to a specific event type (as opposed to differential responses) then a null event must be included (even if it is not modeled explicitly).'
                     ''
                     'In SPM, analysis of data from multiple subjects typically proceeds in two stages using models at two ''levels''. The ''first level'' models are used to implement a within-subject analysis. Typically there will be as many first level models as there are subjects. Analysis proceeds as described using the ''Specify first level'' and ''Estimate'' options. The results of these analyses can then be presented as ''case studies''. More often, however, one wishes to make inferences about the population from which the subjects were drawn. This is an example of a ''Random-Effects (RFX) analysis'' (or, more properly, a mixed-effects analysis). In SPM, RFX analysis is implemented using the ''summary-statistic'' approach where contrast images from each subject are used as summary measures of subject responses. These are then entered as data into a ''second level'' model. '
}';
fmri_spec.prog = @spm_run_fmri_spec;
fmri_spec.vout = @vout_stats;
fmri_spec.modality = {'FMRI'};
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function t = cond_check(job)
t   = {};
if (numel(job.onset) ~= numel(job.duration)) && (numel(job.duration)~=1),
    t = {sprintf('"%s": Number of event onsets (%d) does not match the number of durations (%d).',...
        job.name, numel(job.onset),numel(job.duration))};
end;
for i=1:numel(job.pmod),
    if numel(job.onset) ~= numel(job.pmod(i).param),
        t = {t{:}, sprintf('"%s" & "%s":Number of event onsets (%d) does not equal the number of parameters (%d).',...
            job.name, job.pmod(i).name, numel(job.onset),numel(job.pmod(i).param))};
    end;
end;
return;
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function t = sess_check(sess)
t = {};
for i=1:numel(sess.regress),
    if numel(sess.scans) ~= numel(sess.regress(i).val),
        t = {t{:}, sprintf('Num scans (%d) ~= Num regress[%d] (%d).',numel(sess.scans),i,numel(sess.regress(i).val))};
    end;
end;
return;
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function dep = vout_stats(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'SPM.mat File';
dep(1).src_output = substruct('.','spmmat');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
%dep(2)            = cfg_dep;
%dep(2).sname      = 'SPM Variable';
%dep(2).src_output = substruct('.','spmvar');
%dep(2).tgt_spec   = cfg_findspec({{'strtype','e'}});
