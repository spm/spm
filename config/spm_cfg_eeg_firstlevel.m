function convmodel = spm_cfg_eeg_firstlevel
% SPM Configuration file for M/EEG convolution modelling
%_______________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_eeg_firstlevel.m 6818 2016-06-21 09:42:45Z peter $

rev = '$Rev: 6818 $';

% ---------------------------------------------------------------------
% units Units for design
% ---------------------------------------------------------------------
units         = cfg_menu;
units.tag     = 'units';
units.name    = 'Units for design';
units.help    = {'The onsets of events can be specified in either samples or seconds.'};
units.labels = {
                'Samples'
                'Seconds'
}';
units.values = {
                'samples'
                'secs'
}';
% ---------------------------------------------------------------------
% Time window
% ---------------------------------------------------------------------
timewin         = cfg_entry;
timewin.tag     = 'timewin';
timewin.name    = 'Time window';
timewin.strtype = 'r';
timewin.num     = [1 2];
timewin.help    = {'Start and end of epoch [ms].'};
% ---------------------------------------------------------------------
% Microtime resolution
% ---------------------------------------------------------------------
utime         = cfg_entry;
utime.tag     = 'utime';
utime.name    = 'Microtime resolution';
utime.help    = {
                  'The microtime resolution, t, is the number of time-bins per input sample used when building regressors. '
                  'Can be modified to make the output up- or downsampled with respect to the input'
}';
utime.strtype = 'r';
utime.num     = [1 1];
utime.val     = {1};

% ---------------------------------------------------------------------
% timing Timing parameters
% ---------------------------------------------------------------------
timing         = cfg_branch;
timing.tag     = 'timing';
timing.name    = 'Timing parameters';
timing.val     = {timewin units utime};

% ---------------------------------------------------------------------
% D M/EEG datasets
% ---------------------------------------------------------------------

D = cfg_files;
D.tag = 'D';
D.name = 'M/EEG dataset';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the M/EEG mat file'};

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
onset.strtype = 'r';
onset.num     = [Inf 1];
% ---------------------------------------------------------------------
% duration Durations
% ---------------------------------------------------------------------
duration         = cfg_entry;
duration.tag     = 'duration';
duration.name    = 'Durations';
duration.help    = {'Specify the event durations. Epoch and event-related responses are modeled in exactly the same way but by specifying their different durations.  Events are specified with a duration of 0.  If you enter a single number for the durations it will be assumed that all trials conform to this duration. If you have multiple different durations, then the number must match the number of onset times.'};
duration.strtype = 'r';
duration.num     = [Inf 1];
%

% ---------------------------------------------------------------------
% manual Specification
% ---------------------------------------------------------------------
manual         = cfg_branch;
manual.tag     = 'manual';
manual.name    = 'Specify manually';
manual.val     = {onset duration};

eventtype         = cfg_entry;
eventtype.tag     = 'eventtype';
eventtype.name    = 'Event type';
eventtype.strtype = 's';

eventvalue         = cfg_entry;
eventvalue.tag     = 'eventvalue';
eventvalue.name    = 'Event value';
eventvalue.strtype = 'e';

trlshift         = cfg_entry;
trlshift.tag     = 'trlshift';
trlshift.name    = 'Shift';
trlshift.strtype = 'r';
trlshift.num     = [1 1];
trlshift.val     = {0};
trlshift.help    = {'shift the triggers by a fixed amount (ms)',... 
                   '(e.g. projector delay).'};

event      = cfg_branch;
event.tag  = 'event';
event.name = 'Event';
event.val  = {eventtype eventvalue trlshift};

% ---------------------------------------------------------------------
% 
% ---------------------------------------------------------------------
fromdata         = cfg_repeat;
fromdata.tag     = 'fromdata';
fromdata.name    = 'Take from dataset';
fromdata.values  = {event};
fromdata.num     = [1 Inf];
% ---------------------------------------------------------------------
% bases How to define events
% ---------------------------------------------------------------------
define         = cfg_choice;
define.tag     = 'define';
define.name    = 'How to define events';
define.val     = {fromdata};
define.help    = {'Choose the way to specify events for building regressors.'};
define.values  = {manual fromdata};

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
param.strtype = 'r';
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
                'Model interactions with user specified parameters. This allows nonlinear effects relating to some other measure to be modelled in the design matrix.'
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
% porth Orthogonalise modulations
% ---------------------------------------------------------------------
porth         = cfg_menu;
porth.tag     = 'orth';
porth.name    = 'Orthogonalise modulations';
porth.help    = {'Orthogonalise regressors within trial types.'};
porth.labels  = {'Yes' 'No'};
porth.values  = {1 0};
porth.val     = {1};
% ---------------------------------------------------------------------
% cond Condition
% ---------------------------------------------------------------------
cond         = cfg_branch;
cond.tag     = 'cond';
cond.name    = 'Condition';
cond.val     = {name define tmod generic2 porth};
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
                 'This *.mat file must include the following cell arrays (each 1 x n): names, onsets and durations. eg. names=cell(1,5), onsets=cell(1,5), durations=cell(1,5), then names{2}=''SSent-DSpeak'', onsets{2}=[3 5 19 222], durations{2}=[0 0 0 0], contain the required details of the second condition. These cell arrays may be made available by your stimulus delivery program, eg. COGENT. The duration vectors can contain a single entry if the durations are identical for all events. Optionally, a (1 x n) cell array named orth can also be included, with a 1 or 0 for each condition to indicate whether parameteric modulators should be orthogonalised.'
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
val.strtype = 'r';
val.num     = [Inf 1];
% ---------------------------------------------------------------------
% convregress Convolution Regressor
% ---------------------------------------------------------------------
convregress         = cfg_branch;
convregress.tag     = 'convregress';
convregress.name    = 'Regressor';
convregress.val     = {name val };
convregress.help    = {'Specification for convolution regressor'};
% ---------------------------------------------------------------------
% generic Regressors
% ---------------------------------------------------------------------
generic3         = cfg_repeat;
generic3.tag     = 'generic';
generic3.name    = 'Convolution regressors';
generic3.help    = {'Convolution regressors are continuous variables that are convolved with a basis set to estimate an impulse-response pattern.'};
generic3.values  = {convregress};
generic3.num     = [0 Inf];
% ---------------------------------------------------------------------
% multi_conv_reg Multiple convolution regressors
% ---------------------------------------------------------------------
multi_conv_reg         = cfg_files;
multi_conv_reg.tag     = 'multi_conv_reg';
multi_conv_reg.name    = 'Multiple convolution regressors';
multi_conv_reg.val{1} = {''};
multi_conv_reg.help    = {
                     'Select the *.mat/*.txt file containing details of your multiple regressors. '
                     ''
                     'If you have multiple regressors eg. realignment parameters, then entering the details a regressor at a time is very inefficient. This option can be used to load all the required information in one go. '
                     ''
                     'You will first need to create a *.mat file containing a matrix R or a *.txt file containing the regressors. Each column of R will contain a different regressor. When SPM creates the design matrix the regressors will be named R1, R2, R3, ..etc.'
}';
multi_conv_reg.filter = 'mat';
multi_conv_reg.ufilter = '.*';
multi_conv_reg.num     = [0 1];
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
% Save regressor coefficients
% ---------------------------------------------------------------------
savereg         = cfg_menu;
savereg.tag     = 'savereg';
savereg.name    = 'Save regressor coefficients';
savereg.help    = {'Choose ''yes'' to save the coefficients for regressors as a separate dataset (of spectra for TF data). If you are only using regressors to model out nuisance variables',
    '(e.g. motion) saving might not be necessary'};
               
savereg.labels = {'yes', 'no'};
savereg.values = {true, false};
savereg.val    = {false};
% ---------------------------------------------------------------------
% hpf High-pass filter
% ---------------------------------------------------------------------
hpf         = cfg_entry;
hpf.tag     = 'hpf';
hpf.name    = 'High-pass filter';
hpf.help    = {'The default high-pass filter cutoff is 10 seconds.Slow signal drifts with a period longer than this will be removed. Use ''explore design'' to ensure this cut-off is not removing too much experimental variance. High-pass filtering is implemented using a residual forming matrix (i.e. it is not a convolution) and is simply to a way to remove confounds without estimating their parameters explicitly.  The constant term is also incorporated into this filter matrix.'};
hpf.strtype = 'r';
hpf.num     = [1 1];
hpf.val     = {10};
% ---------------------------------------------------------------------
% sess Subject/Session
% ---------------------------------------------------------------------
sess         = cfg_branch;
sess.tag     = 'sess';
sess.name    = 'Subject/Session';
sess.val     = {D generic1 multi generic3 multi_conv_reg generic2 multi_reg savereg hpf };
sess.help    = {'The design matrix consists of one or more separable, session-specific partitions.  These partitions are usually either one per subject, or one per scanning session for that subject.'};
% ---------------------------------------------------------------------
% order Order
% ---------------------------------------------------------------------
order         = cfg_entry;
order.tag     = 'order';
order.name    = 'Order';
order.help    = {'Number of basis functions'};
order.strtype = 'n';
order.val     = {12};
order.num     = [1 1];
% ---------------------------------------------------------------------
% fourier Fourier Set
% ---------------------------------------------------------------------
fourier         = cfg_branch;
fourier.tag     = 'fourier';
fourier.name    = 'Fourier Set';
fourier.val     = {order};
fourier.help    = {'Fourier basis functions.'};
% ---------------------------------------------------------------------
% order Order
% ---------------------------------------------------------------------
order         = cfg_entry;
order.tag     = 'order';
order.name    = 'Order';
order.help    = {'Number of basis functions'};
order.strtype = 'n';
order.num     = [1 1];
% ---------------------------------------------------------------------
% fourier_han Fourier Set (Hanning)
% ---------------------------------------------------------------------
fourier_han         = cfg_branch;
fourier_han.tag     = 'fourier_han';
fourier_han.name    = 'Fourier Set (Hanning)';
fourier_han.val     = {order};
fourier_han.help    = {'Fourier basis functions with Hanning Window'};

% ---------------------------------------------------------------------
% order Order
% ---------------------------------------------------------------------
order         = cfg_entry;
order.tag     = 'order';
order.name    = 'Order';
order.help    = {'Number of basis functions'};
order.strtype = 'n';
order.num     = [1 1];
% ---------------------------------------------------------------------
% fir Finite Impulse Response
% ---------------------------------------------------------------------
fir         = cfg_branch;
fir.tag     = 'fir';
fir.name    = 'Finite Impulse Response';
fir.val     = {order};
fir.help    = {'Finite impulse response.'};
% ---------------------------------------------------------------------
% bases Basis Functions
% ---------------------------------------------------------------------
bases         = cfg_choice;
bases.tag     = 'bases';
bases.name    = 'Basis Functions';
bases.val     = {fourier};
bases.help    = {'Choose the basis set'};
bases.values  = {fourier fourier_han fir};
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

%--------------------------------------------------------------------------
% prefix
%--------------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the output dataset. Default prefix is ''C''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'C'};

% ---------------------------------------------------------------------
% eeg_design MEEG model specification 
% ---------------------------------------------------------------------
convmodel         = cfg_exbranch;
convmodel.tag     = 'convmodel';
convmodel.name    = 'Convolution modelling';
convmodel.val     = {spm_cfg_eeg_channel_selector timing sess bases volt, prefix};
convmodel.prog = @eeg_run;
convmodel.vout = @vout_eeg;
convmodel.modality = {'EEG'};

%-------------------------------------------------------------------------
function out = eeg_run(job)
S = job;
S.channels  = spm_cfg_eeg_channel_selector(job.channels);
out.D       = spm_eeg_firstlevel(S);
out.Dfname  = {fullfile(out.D)};
%-------------------------------------------------------------------------
function dep = vout_eeg(job)
% return dependencies
dep(1)            = cfg_dep;
dep(1).sname      = 'Estimated M/EEG data';
dep(1).src_output = substruct('.','D');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'Estimated M/EEG datafile';
dep(2).src_output = substruct('.','Dfname');
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});