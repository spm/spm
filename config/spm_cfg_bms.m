function bms = spm_cfg_bms
% Configuration file for BMS interface.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Maria Joao Rosa
% $Id: spm_cfg_bms.m 3955 2010-06-29 17:26:29Z maria $

% ---------------------------------------------------------------------
% dir Directory
% ---------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.help    = {['Select the directory where the files containing the '...
               'results from BMS (BMS.mat) will be written.']};
dir.filter  = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];

% ---------------------------------------------------------------------
% mod_dcm Models (.mat)
% ---------------------------------------------------------------------
mod_dcm         = cfg_files;
mod_dcm.tag     = 'mod_dcm';
mod_dcm.name    = 'Models';
mod_dcm.help    = {['Select the DCM file (.mat) for each model. '...
                    'DCM.mat files (models) should be specified '...
                    'in the same order for each subject and session.']};
mod_dcm.filter  = 'mat';
mod_dcm.ufilter = '.*';
mod_dcm.num     = [0 Inf];

% ---------------------------------------------------------------------
% sess_dcm Sessions
% ---------------------------------------------------------------------
sess_dcm      = cfg_branch;
sess_dcm.tag  = 'sess_dcm';
sess_dcm.name = 'Session';
sess_dcm.val  = {mod_dcm };

% ---------------------------------------------------------------------
% subj_dcm Subject
% ---------------------------------------------------------------------
subj_dcm         = cfg_repeat;
subj_dcm.tag     = 'subj_dcm';
subj_dcm.name    = 'Subject';
subj_dcm.values  = {sess_dcm };

% ---------------------------------------------------------------------
% dcm Data
% ---------------------------------------------------------------------
dcm         = cfg_repeat;
dcm.tag     = 'dcm';
dcm.name    = 'Data';
dcm.help    = {['Select DCM file (.mat) for each model, session and '...
               'subject.']}';
dcm.values  = {subj_dcm };
dcm.num     = [0 Inf];

% ---------------------------------------------------------------------
% mod_map Models (.img)
% ---------------------------------------------------------------------
mod_map         = cfg_files;
mod_map.tag     = 'mod_map';
mod_map.name    = 'Models';
mod_map.help    = {['Specify the log. evidence map (.img) for each model. '...
                    'Log-evidence maps should be specified '...
                    'in the same order for each subject and session.']};
mod_map.filter  = 'image';
mod_map.ufilter = '.*';
mod_map.num     = [1 Inf];

% ---------------------------------------------------------------------
% sess_map Sessions (Maps)
% ---------------------------------------------------------------------
sess_map      = cfg_branch;
sess_map.tag  = 'sess_map';
sess_map.name = 'Session';
sess_map.val  = {mod_map };

% ---------------------------------------------------------------------
% subj_dcm Subject (Maps)
% ---------------------------------------------------------------------
subj_map         = cfg_repeat;
subj_map.tag     = 'subj_map';
subj_map.name    = 'Subject';
subj_map.values  = {sess_map };

% ---------------------------------------------------------------------
% map Data
% ---------------------------------------------------------------------
map         = cfg_repeat;
map.tag     = 'map';
map.name    = 'Data';
map.help    = {['Select the log. evidence maps (.img) for each '...
               'model, session and subject.']}';
map.values  = {subj_map };
map.num     = [1 Inf];

% ---------------------------------------------------------------------
% mod_name Name
% ---------------------------------------------------------------------
mod_name         = cfg_entry;
mod_name.tag     = 'mod_name';
mod_name.name    = 'Name';
mod_name.help    = {'Specify name for each model (optional).'};
mod_name.strtype = 's';
mod_name.num     = [0 Inf];
mod_name.val     = {''};

% ---------------------------------------------------------------------
% name_mod Name models
% ---------------------------------------------------------------------
name_mod         = cfg_repeat;
name_mod.tag     = 'name_mod';
name_mod.name    = 'Name models';
name_mod.help    = {'Specify name for each model (optional).'}';
name_mod.values  = {mod_name };
name_mod.num     = [0 Inf];

% ---------------------------------------------------------------------
% model_sp Load model space
% ---------------------------------------------------------------------
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

% ---------------------------------------------------------------------
% load_f Log-evidence matrix
% ---------------------------------------------------------------------
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

% ---------------------------------------------------------------------
% method Inference Method
% ---------------------------------------------------------------------
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

% ---------------------------------------------------------------------
% method_maps Inference Method (maps)
% ---------------------------------------------------------------------
method_maps         = cfg_menu;
method_maps.tag     = 'method_maps';
method_maps.name    = 'Inference method';
method_maps.help    = {['Specify inference method: random effects '...
                   '(2nd-level, RFX) or fixed effects (1st-level, FFX) analysis. '...
                   'RFX uses a Variational Bayes approach.']};
method_maps.labels  = {
                  'Fixed effects (FFX)'
                  'Random effects (RFX)'
}';
method_maps.values  = {
                  'FFX'
                  'RFX'
}'; 

% % ---------------------------------------------------------------------
% % priors Priors
% % ---------------------------------------------------------------------
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


% ---------------------------------------------------------------------
% family_file Family file
% ---------------------------------------------------------------------
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
family_file.val{1}  = {''};
family_file.filter  = 'mat';
family_file.ufilter = '.*';
family_file.num     = [0 1];

% ---------------------------------------------------------------------
% family_name Family name
% ---------------------------------------------------------------------
family_name         = cfg_entry;
family_name.tag     = 'family_name';
family_name.name    = 'Name';
family_name.help    = {'Specify name for family.'};
family_name.strtype = 's';
family_name.num     = [0 Inf];

% ---------------------------------------------------------------------
% family_models family_models
% ---------------------------------------------------------------------
family_models         = cfg_entry;
family_models.tag     = 'family_models';
family_models.name    = 'Models';
family_models.help    = {['Specify models belonging to this family. '...
                          'Example: write ''2 6'' if the second and sixth model '...
                          'belong to this family.']};
family_models.strtype = 'e';
family_models.num     = [Inf 1];

% ---------------------------------------------------------------------
% family Family
% ---------------------------------------------------------------------
family         = cfg_branch;
family.tag     = 'family';
family.name    = 'Family';
family.val     = {family_name family_models };
family.help    = {'Specify family name and models.'};

% ---------------------------------------------------------------------
% select_family Specify family
% ---------------------------------------------------------------------
select_family         = cfg_repeat;
select_family.tag     = 'select_family';
select_family.name    = 'Construct family';
select_family.values  = {family };
select_family.help    = {'Create family. Specify family name and models.'};

% ---------------------------------------------------------------------
% family_level Specify families
% ---------------------------------------------------------------------
family_level         = cfg_choice;
family_level.tag     = 'family_level';
family_level.name    = 'Family inference';
family_level.help    = {['Optional field to perform family level inference.'...
                         'Options: load family.mat '...
                         'or specify family names and models using '...
                         'the interface.']};
family_level.val     = {family_file };
family_level.values  = {family_file select_family };

% ---------------------------------------------------------------------
% bma_part Choose family
% ---------------------------------------------------------------------
bma_part         = cfg_entry;
bma_part.tag     = 'bma_part';
bma_part.name    = 'Enter family';
bma_part.help    = {['Specify family (integer). E.g. ''2'' for the second '...
                    'family to use in BMA. ']};
bma_part.strtype = 'e';
bma_part.num     = [0 Inf];

% ---------------------------------------------------------------------
% bma_no no
% ---------------------------------------------------------------------
bma_all         = cfg_const;
bma_all.tag     = 'bma_all';
bma_all.name    = 'All families';
bma_all.val     = {'famwin'};
bma_all.help    = {'Use all families for Bayesian Model Averaging (BMA).'}';

% ---------------------------------------------------------------------
% bma_no no
% ---------------------------------------------------------------------
bma_famwin         = cfg_const;
bma_famwin.tag     = 'bma_famwin';
bma_famwin.name    = 'Winning family';
bma_famwin.val     = {'fanwin'};
bma_famwin.help    = {'Use winning family for Bayesian Model Averaging (BMA).'}';

% ---------------------------------------------------------------------
% bma_no no
% ---------------------------------------------------------------------
bma_no         = cfg_const;
bma_no.tag     = 'bma_no';
bma_no.name    = 'Do not compute';
bma_no.val     = {0};
bma_no.help    = {'Do not compute Bayesian Model Averaging (BMA).'}';

% ---------------------------------------------------------------------
% bma_yes BMA set
% ---------------------------------------------------------------------
bma_yes         = cfg_choice;
bma_yes.tag     = 'bma_yes';
bma_yes.name    = 'Choose family';
bma_yes.help    = {['Specify family for Bayesian Model Averaging (BMA). '...
                    'Options: ''winning family'', ''enter family'' or '...
                    '''all families''.']};
bma_yes.val     = {bma_famwin };
bma_yes.values  = {bma_famwin bma_all bma_part };

% ---------------------------------------------------------------------
% bma BMA
% ---------------------------------------------------------------------
bma         = cfg_choice;
bma.tag     = 'bma';
bma.name    = 'BMA';
bma.help    = {'Optional field to compute Bayesian Model Averaging (BMA).'};
bma.val     = {bma_no };
bma.values  = {bma_no bma_yes };

% ---------------------------------------------------------------------
% verify_id Verify data ID
% ---------------------------------------------------------------------
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
verify_id.val     = {0};

% ---------------------------------------------------------------------
% out_file Output files
% ---------------------------------------------------------------------
out_file         = cfg_menu;
out_file.tag     = 'out_file';
out_file.name    = 'Output files (RFX)';
out_file.help    = {['Specify which output files to save (only valid for'...
                     'RFX analyses). ']...
                     ''...
                    ['Default option (and faster option): '...
                     'PPM = xppm.img (Posterior Probability Maps) '...
                     'for each model. ']...
                     ''...
                    ['Second option: PPM + EPM = xppm.img + '...
                     'epm.img (Exceedance Probability '...
                     'Maps + Exceedance Probability Maps) for each model.']...
                     ''...
                    ['Third option: PPM + EPM + Alpha = xppm.img + '...
                     'epm.img + alpha.img (PPM, EPM and Map of Dirichlet '...
                     'Parameters) for each model.']};
out_file.labels  = {
                   'PPM'
                   'PPM + EPM'
                   'PPM + EPM + Alpha'
                   
}';
out_file.values  = {
                  0
                  1
                  2
}';
out_file.val     = {0};

% ---------------------------------------------------------------------
% mask Mask Image
% ---------------------------------------------------------------------
mask         = cfg_files;
mask.tag     = 'mask';
mask.name    = 'Mask Image';
mask.help    = {['Specify an image for explicitly masking the analysis. '...
                '(optional). '...
                'A sensible option here is to use a segmention of '...
                'structural images to specify a within-brain mask. '...
                'If you select that image as an explicit mask then only '...
                'those voxels in the brain will be analysed. This both '...
                'speeds the inference process and restricts BMS to '...
                'within-brain voxels. Alternatively, if such structural '...
                'images are unavailble or no masking is required, then '...
                'leave this field empty.']};
mask.filter  = 'image';
mask.ufilter = '.*';
mask.val     = {{''}};
mask.num     = [0 1];

% ---------------------------------------------------------------------
% nsamp Number of samples
% ---------------------------------------------------------------------
nsamp         = cfg_entry;
nsamp.tag     = 'nsamp';
nsamp.name    = 'Number of samples';
nsamp.help    = {['Number of samples used to compute exceedance '...
                  'probabilities (default: 1e6). '...
                  'To make computations faster reduce the number of '...
                  'samples when number of models is bigger than 3.']};                 
nsamp.strtype = 's';
nsamp.num     = [1 Inf];
nsamp.val     = {'1e6'};

% ---------------------------------------------------------------------
% file BMS.mat
% ---------------------------------------------------------------------
file         = cfg_files;
file.tag     = 'file';
file.name    = 'BMS.mat';
file.help    = {['Specify the BMS (.mat) file obtained from previous BMS '...
               'analysis (optional). Leave field empty to work on '...
               'serial mode.']};
file.filter  = 'mat';
file.ufilter = '.*';
file.val     = {{''}};
file.num     = [0 1];

% ---------------------------------------------------------------------
% img Map to display
% ---------------------------------------------------------------------
img         = cfg_files;
img.tag     = 'img';
img.name    = 'Map to display';
img.help    = {['Specify map (.img) obtained from BMS Maps '...
               '(optional). Leave field empty to work on serial mode.']};
img.filter  = 'image';
img.ufilter = '.*';
img.val     = {{''}};
img.num     = [0 1];

% ---------------------------------------------------------------------
% thres Probability Threshold
% ---------------------------------------------------------------------
thres         = cfg_entry;
thres.tag     = 'thres';
thres.name    = 'Probability threshold';
thres.help    = {['Specify the probability threshold to apply to the '...
                 'image (optional). Leave field empty to work on '...
                 'serial mode.']};                 
thres.strtype = 'e';
thres.num     = [0 Inf];
thres.val     = {[]};

% ---------------------------------------------------------------------
% k Extent threshold
% ---------------------------------------------------------------------
k         = cfg_entry;
k.tag     = 'k';
k.name    = 'Extent threshold';
k.help    = {['Specify extent threshold (minimum number of voxels '...
                 'per cluster).']};                 
k.strtype = 'e';
k.num     = [0 Inf];
k.val     = {[]};

% ---------------------------------------------------------------------
% scale Map Scale
% ---------------------------------------------------------------------
scale         = cfg_menu;
scale.tag     = 'scale';
scale.name    = 'Map scale';
scale.help    = {['Specify scale to display maps (optional). Default: '...
                 'empty field to work on serial mode. Other options: '...
                 '''None'' will display image with original scale and '...
                 '''Log-odds'' will display image in a log-odds '...
                 ' scale (in this case .img should be a '...
                 'probability map).']};
scale.labels  = {
                  'Empty'
                  'None'
                  'Log-odds'
}';
scale.values  = {
                  []
                  0
                  1
}';
scale.val     = {[]};

% ---------------------------------------------------------------------
% bms_dcm BMS: DCM, output is bar plot
% ---------------------------------------------------------------------
bms_dcm      = cfg_exbranch;
bms_dcm.tag  = 'bms_dcm';
bms_dcm.name = 'BMS: DCM';
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
bms_dcm.prog = @spm_run_bms_dcm;
bms_dcm.vout = @vout;

% ---------------------------------------------------------------------
% bms_dcm_vis: DCM (visualise results)
% ---------------------------------------------------------------------
bms_dcm_vis      = cfg_exbranch;
bms_dcm_vis.tag  = 'bms_dcm_vis';
bms_dcm_vis.name = 'BMS: DCM (Results)';
bms_dcm_vis.val  = {file };
bms_dcm_vis.help = {['Bayesian Model Selection for DCM (Results). '...
                    'Show results from BMS for DCM.']};
bms_dcm_vis.prog = @spm_run_bms_dcm_vis;

% ---------------------------------------------------------------------
% bms_map_inf BMS: Maps (Inference), output is BMS map 
% ---------------------------------------------------------------------
bms_map_inf      = cfg_exbranch;
bms_map_inf.tag  = 'bms_map_inf';
bms_map_inf.name = 'BMS: Maps';
bms_map_inf.val  = {dir map name_mod method_maps out_file mask nsamp };
bms_map_inf.help = {'Bayesian Model Selection for Log-Evidence Maps. '...
    ''...
    ['Input: log-evidence maps (.img) for each model, session and '...
    'subject. Note that there must be identical numbers of models for '...
    'all sessions, and identical numbers of sessions for all '...
    'subjects.']...
    ''...
    ['Output: For the fixed effects analysis, posterior probability maps '...
    '(.img) are created for each model. '...
    'For the random effects analysis, expected posterior probability '...
    'and exceedance probability (i.e. the probability that this model '...
    'is more likely than any other model) maps are created for each '...
    'model. If there are multiple sessions per subject, the random '...
    'effects analysis operates on the subject-specific sums of log '...
    'evidences across sessions. In addition, a BMS.mat file will be save '...
    'in the specified directory for both methods']};
bms_map_inf.prog = @spm_run_bms_map;
bms_map_inf.vout = @vout;

% ---------------------------------------------------------------------
% bms_map_vis BMS: Maps (Results), visualisation of BMS Maps results
% ---------------------------------------------------------------------
bms_map_vis      = cfg_exbranch;
bms_map_vis.tag  = 'bms_map_vis';
bms_map_vis.name = 'BMS: Maps (Results)';
bms_map_vis.val  = {file img thres k scale};
bms_map_vis.help = {['Bayesian Model Selection Maps (Results). '...
                    'Show results from BMS Maps (Inference).']};
bms_map_vis.prog = @spm_run_bms_vis;

% ---------------------------------------------------------------------
% bms Bayesian Model Selection
% ---------------------------------------------------------------------
bms         = cfg_choice;
bms.tag     = 'bms';
bms.name    = 'Bayesian Model Selection';
bms.help    = {['Bayesian Model Selection for group studies (fixed '...
               'effects and random effects analysis).']};
bms.values  = {bms_dcm bms_dcm_vis bms_map_inf bms_map_vis };

%------------------------------------------------------------------------
function dep = vout(varargin)
% Output file names will be saved in a struct with field .files
dep(1)            = cfg_dep;
dep(1).sname      = 'BMS Maps';
dep(1).src_output = substruct('.','files');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});


