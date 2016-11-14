function bms = spm_cfg_dcm_bms
% Configuration file for Bayesian Model Selection (DCM)
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Maria Joao Rosa
% $Id: spm_cfg_dcm_bms.m 6929 2016-11-14 13:07:31Z guillaume $

%--------------------------------------------------------------------------
% dir Directory
%--------------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.help    = {['Select the directory where the files containing the '...
               'results from BMS (BMS.mat) will be written.']};
dir.filter  = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];

%--------------------------------------------------------------------------
% dcmmat Models
%--------------------------------------------------------------------------
dcmmat         = cfg_files;
dcmmat.tag     = 'dcmmat';
dcmmat.name    = 'Models';
dcmmat.help    = {['Select the DCM_*.mat file for each model. '...
                   'DCM_*.mat files (models) should be specified '...
                   'in the same order for each subject and session.']};
dcmmat.filter  = 'mat';
dcmmat.ufilter = '^DCM.*\.mat$';
dcmmat.num     = [0 Inf];

%--------------------------------------------------------------------------
% sess_dcm Sessions
%--------------------------------------------------------------------------
sess_dcm      = cfg_branch;
sess_dcm.tag  = 'sess_dcm';
sess_dcm.name = 'Session';
sess_dcm.val  = { dcmmat };

%--------------------------------------------------------------------------
% subj_dcm Subject
%--------------------------------------------------------------------------
subj_dcm         = cfg_repeat;
subj_dcm.tag     = 'subj_dcm';
subj_dcm.name    = 'Subject';
subj_dcm.values  = {sess_dcm };

%--------------------------------------------------------------------------
% dcm Data
%--------------------------------------------------------------------------
dcm         = cfg_repeat;
dcm.tag     = 'dcm';
dcm.name    = 'Data';
dcm.help    = {['Select the DCM_*.mat file for each model, session and '...
               'subject.']}';
dcm.values  = {subj_dcm };
dcm.num     = [0 Inf];

%--------------------------------------------------------------------------
% model_sp Load model space
%--------------------------------------------------------------------------
model_sp         = cfg_files;
model_sp.tag     = 'model_sp';
model_sp.name    = 'Load model space';
model_sp.help    = {['Optional: load .mat file with all subjects, sessions '...
                  'and models. This option is a faster alternative to selecting '...
                  'the DCM.mat files for each subject/model (above in '...
                  '''Data'').']
                  ['This file is created if the ''Data'' option has been used. '...
                  'It is saved in the same directory as BMS.mat and can then be loaded '...
                  'for future BMS/BMA analyses with the same data.']
                  ['The model space file should contain the structure ''subj''. ' ...
                  'This structure should have the field ''sess'' for sessions, '...
                  'then the subfield ''model'' and in ''model'' there should be '...
                  'five subfields: ''fname'' contains the path to the DCM.mat file, '...
                  '''.F'' the Free Energy of that model, '...
                  '''.Ep'' and ''Cp'' the mean and covariance of the parameters estimates. '...
                  'Finally the subfield ''.nonLin'' should be 1 if the model is non-linear and '...
                  '0 otherwise.']
                  ['Example: subj(3).sess(1).model(4).fname contains the path to the DCM.mat '...
                  'file for subject 3, session 1 and model 4. subj(3).sess(1).model(4).F '...
                  'contains the value of the Free Energy for the same model/session/subject.']};
model_sp.filter  = 'mat';
model_sp.ufilter = '.*';
model_sp.val     = {{''}};
model_sp.num     = [0 1];

%--------------------------------------------------------------------------
% load_f Log-evidence matrix
%--------------------------------------------------------------------------
load_f         = cfg_files;
load_f.tag     = 'load_f';
load_f.name    = 'Log-evidence matrix';
load_f.help    = {['Optional: load .mat file with log-evidence values for '...
                  'comparison. This option is a faster alternative to selecting '...
                  'the DCM.mat files for each subject/model (above in '...
                  '''Data'') but it does not allow for Bayesian Model Averaging. '...
                  'To compute BMA the user needs to specify the DCM.mat files '...
                  'or the model space. ']
                  ['This file should contain an F matrix consisting ' ...
                  'of [s x m] log-evidence values, where s is the number '...
                  'of subjects and m the number of models.']};
load_f.filter  = 'mat';
load_f.ufilter = '.*';
load_f.val     = {{''}};
load_f.num     = [0 1];

%--------------------------------------------------------------------------
% method Inference Method
%--------------------------------------------------------------------------
method         = cfg_menu;
method.tag     = 'method';
method.name    = 'Inference method';
method.help    = {['Specify inference method: random effects '...
                   '(2nd-level, RFX) or fixed effects (1st-level, FFX) analysis. '...
                   'RFX uses Gibbs sampling.']};
method.labels  = {
                  'Fixed effects (FFX)'
                  'Random effects (RFX)'
}';
method.values  = {
                  'FFX'
                  'RFX'
}'; 

% %--------------------------------------------------------------------------
% % priors Priors
% %--------------------------------------------------------------------------
% priors         = cfg_menu;
% priors.tag     = 'priors';
% priors.name    = 'Priors';
% priors.help    = {['Specify priors for family-level inference (RFX only).
% '...
%                    'Options: ''Family'' sets alpha0=1 for each family '...
%                    'while ''Model'' sets alpha0=1 for each model (not '...
%                    'advised).']};
% priors.labels  = {
%                   'Model'
%                   'Family'
% }';
% priors.values  = {
%                   'M-unity'
%                   'F-unity'
% }';
% priors.val      = {'F-unity'};


%--------------------------------------------------------------------------
% family_file Family file
%--------------------------------------------------------------------------
family_file         = cfg_files;
family_file.tag     = 'family_file';
family_file.name    = 'Load family';
family_file.help    = {['Load family.mat file. This file should contain the '...
                        'structure ''family'' with fields ''names'' and '...
                        '''partition''. Example: family.names = {''F1'', '...
                        '''F2''} and family.partition = [1 2 2 1 1]. '...
                        ' This structure specifies two families with names '...
                        '''F1'' and ''F2'' and assigns model 1, 4 and 5 to '...
                        'the first family and models 2 and 3 to the second '...
                        'family.']};
family_file.val     = {{''}};
family_file.filter  = 'mat';
family_file.ufilter = '.*';
family_file.num     = [0 1];

%--------------------------------------------------------------------------
% family_name Family name
%--------------------------------------------------------------------------
family_name         = cfg_entry;
family_name.tag     = 'family_name';
family_name.name    = 'Name';
family_name.help    = {'Specify name for family.'};
family_name.strtype = 's';
family_name.num     = [0 Inf];

%--------------------------------------------------------------------------
% family_models family_models
%--------------------------------------------------------------------------
family_models         = cfg_entry;
family_models.tag     = 'family_models';
family_models.name    = 'Models';
family_models.help    = {['Specify models belonging to this family. '...
                          'Example: write ''2 6'' if the second and sixth model '...
                          'belong to this family.']};
family_models.strtype = 'n';
family_models.num     = [Inf 1];

%--------------------------------------------------------------------------
% family Family
%--------------------------------------------------------------------------
family         = cfg_branch;
family.tag     = 'family';
family.name    = 'Family';
family.val     = {family_name family_models };
family.help    = {'Specify family name and models.'};

%--------------------------------------------------------------------------
% select_family Specify family
%--------------------------------------------------------------------------
select_family         = cfg_repeat;
select_family.tag     = 'select_family';
select_family.name    = 'Construct family';
select_family.values  = {family };
select_family.help    = {'Create family. Specify family name and models.'};

%--------------------------------------------------------------------------
% family_level Specify families
%--------------------------------------------------------------------------
family_level         = cfg_choice;
family_level.tag     = 'family_level';
family_level.name    = 'Family inference';
family_level.help    = {['Optional field to perform family level inference.'...
                         'Options: load family.mat '...
                         'or specify family names and models using '...
                         'the interface.']};
family_level.val     = {family_file };
family_level.values  = {family_file select_family };

%--------------------------------------------------------------------------
% bma_part Choose family
%--------------------------------------------------------------------------
bma_part         = cfg_entry;
bma_part.tag     = 'bma_part';
bma_part.name    = 'Enter family';
bma_part.help    = {['Specify family (integer). E.g. ''2'' for the second '...
                    'family to use in BMA. ']};
bma_part.strtype = 'n';
bma_part.num     = [0 Inf];

%--------------------------------------------------------------------------
% bma_no no
%--------------------------------------------------------------------------
bma_all         = cfg_const;
bma_all.tag     = 'bma_all';
bma_all.name    = 'All families';
bma_all.val     = {'famwin'};
bma_all.help    = {'Use all families for Bayesian Model Averaging (BMA).'}';

%--------------------------------------------------------------------------
% bma_no no
%--------------------------------------------------------------------------
bma_famwin         = cfg_const;
bma_famwin.tag     = 'bma_famwin';
bma_famwin.name    = 'Winning family';
bma_famwin.val     = {'famwin'};
bma_famwin.help    = {'Use winning family for Bayesian Model Averaging (BMA).'}';

%--------------------------------------------------------------------------
% bma_no no
%--------------------------------------------------------------------------
bma_no         = cfg_const;
bma_no.tag     = 'bma_no';
bma_no.name    = 'Do not compute';
bma_no.val     = {0};
bma_no.help    = {'Do not compute Bayesian Model Averaging (BMA).'}';

%--------------------------------------------------------------------------
% bma_yes BMA set
%--------------------------------------------------------------------------
bma_yes         = cfg_choice;
bma_yes.tag     = 'bma_yes';
bma_yes.name    = 'Choose family';
bma_yes.help    = {['Specify family for Bayesian Model Averaging (BMA). '...
                    'Options: ''winning family'', ''enter family'' or '...
                    '''all families''.']};
bma_yes.val     = {bma_famwin };
bma_yes.values  = {bma_famwin bma_all bma_part };

%--------------------------------------------------------------------------
% bma BMA
%--------------------------------------------------------------------------
bma         = cfg_choice;
bma.tag     = 'bma';
bma.name    = 'BMA';
bma.help    = {'Optional field to compute Bayesian Model Averaging (BMA).'};
bma.val     = {bma_no };
bma.values  = {bma_no bma_yes };

%--------------------------------------------------------------------------
% verify_id Verify data ID
%--------------------------------------------------------------------------
verify_id         = cfg_menu;
verify_id.tag     = 'verify_id';
verify_id.name    = 'Verify data identity';
verify_id.help    = {['Verify whether the model comparison is valid '...
                   'i.e. whether the models have been fitted to the same data.']};
verify_id.labels  = {
                  'Yes'
                  'No'
}';
verify_id.values  = {
                  1
                  0
}'; 
verify_id.val     = {1};

%--------------------------------------------------------------------------
% bmsmat BMS.mat
%--------------------------------------------------------------------------
bmsmat         = cfg_files;
bmsmat.tag     = 'bmsmat';
bmsmat.name    = 'BMS.mat';
bmsmat.help    = {['Specify the BMS.mat file obtained from previous BMS '...
                   'analysis (optional). Leave field empty to work on '...
                   'serial mode.']};
bmsmat.filter  = 'mat';
bmsmat.ufilter = '^BMS\.mat$';
bmsmat.val     = {{''}};
bmsmat.num     = [0 1];

%--------------------------------------------------------------------------
% inference: Model Inference
%--------------------------------------------------------------------------
bms_dcm      = cfg_exbranch;
bms_dcm.tag  = 'inference';
bms_dcm.name = 'Model Inference';
bms_dcm.val  = {dir dcm model_sp load_f method family_level bma verify_id};
bms_dcm.help = {['Bayesian Model Selection for Dynamic Causal Modelling '...
    '(DCM) for fMRI or MEEG.']...
    ''...
    ['Input: DCM files (.mat) for each model, session and subject. '...
    'Note that there must be identical numbers of models for all each '...
    'sessions, and identical numbers of sessions for all subjects. ']...
    ''...
    ['Output: For the fixed effects analysis, the log-evidence for each '...
    'model (relative to the worst model) is plotted in the graphics '...
    'window, as well as the posterior probability for each model. In '...
    'addition, the corresponding values are saved in the directory '...
    'specified (BMS.mat). For the random effects analysis, the '...
    'expected posterior probability and exceedance probability of each '...
    'model (i.e. the probability that this model is more likely than '...
    'any other model) are plotted in the graphics window, and the '...
    'corresponding values are saved in the directory specified. If '...
    'there are multiple sessions per subject, the random effects '...
    'analysis operates on the subject-specific sums of log evidences '...
    'across sessions.']};
bms_dcm.prog = @spm_run_dcm_bms;
bms_dcm.vout = @vout;

%--------------------------------------------------------------------------
% results: Visualise BMS results
%--------------------------------------------------------------------------
bms_dcm_vis      = cfg_exbranch;
bms_dcm_vis.tag  = 'results';
bms_dcm_vis.name = 'Review results';
bms_dcm_vis.val  = { bmsmat };
bms_dcm_vis.help = {['Bayesian Model Selection for DCM (Results). '...
                    'Show results from BMS for DCM.']};
bms_dcm_vis.prog = @spm_run_dcm_bms_vis;

%--------------------------------------------------------------------------
% bms Bayesian Model Selection
%--------------------------------------------------------------------------
bms         = cfg_choice;
bms.tag     = 'bms';
bms.name    = 'Bayesian Model Selection';
bms.help    = {['Bayesian Model Selection for group studies (fixed '...
               'effects and random effects analysis).']};
bms.values  = { bms_dcm bms_dcm_vis };


%==========================================================================
function dep = vout(varargin)
% Output file names will be saved in a struct with field .bmsmat
dep(1)            = cfg_dep;
dep(1).sname      = 'BMS.mat File';
dep(1).src_output = substruct('.','bmsmat');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
