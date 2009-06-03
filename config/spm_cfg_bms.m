function bms = spm_cfg_bms
% Configuration file for BMS interface.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Maria Joao Rosa
% $Id: spm_cfg_bms.m 3177 2009-06-03 08:47:41Z vladimir $

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
mod_dcm.help    = {'Select the DCM file (.mat) for each model.'};
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
mod_map.help    = {'Specify the log. evidence map (.img) for each model.'};
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
% file BMS.mat
% ---------------------------------------------------------------------
load_f         = cfg_files;
load_f.tag     = 'load_f';
load_f.name    = 'Log-evidence Matrix';
load_f.help    = {['Load .mat file with log-evidence values for '...
                  'comparison (optional). The file should contain '...
                  'an F matrix consisting of [s x m] log-evidence '...
                  'values, where s is the number of subjects and m '...
                  'the number of models.']};
load_f.filter  = 'mat';
load_f.ufilter = '.*';
load_f.val     = {{''}};
load_f.num     = [0 1];

% ---------------------------------------------------------------------
% method Inference Method
% ---------------------------------------------------------------------
method         = cfg_menu;
method.tag     = 'method';
method.name    = 'Inference Method';
method.help    = {['Specify inference method: random effects '...
                   '(2nd-level) or fixed effects (1st-level) analysis.']};
method.labels  = {
                  'Fixed effects (FFX)'
                  'Random effects (RFX)'
}';
method.values  = {
                  'FFX'
                  'RFX'
}'; 

% ---------------------------------------------------------------------
% verify_id Verify data ID
% ---------------------------------------------------------------------
verify_id         = cfg_menu;
verify_id.tag     = 'verify_id';
verify_id.name    = 'Verify data identity';
verify_id.help    = {['Verify whether the model comparison is valid '...
                   'i.e. whether the models have been fitted to the same data.']};
verify_id.labels  = {
                  'yes'
                  'no'
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
                     'Alpha + PPM = alpha.img (Dirichlet '...
                     'parameters) + xppm.img (Posterior Probability '...
                     'Maps) for each model.. ']...
                     ''...
                    ['Second option: Alpha + PPM + EPM = alpha.img + '...
                     'xppm.img + epm.img (Exceedance Probability '...
                     ' Maps) for each model.']};
out_file.labels  = {
                   'Alpha + PPM'
                   'Alpha + PPM + EPM'
}';
out_file.values  = {
                  0
                  1
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
bms_dcm.val  = {dir dcm load_f method verify_id};
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
bms_map_inf.val  = {dir map method out_file mask nsamp };
bms_map_inf.help = {'Bayesian Model Selection for Log-Evidence Maps. '...
    ''...
    ['Input: log-evidence maps (.img) for each model, session and '...
    'subject. Note that there must be identical numbers of models for '...
    'all each sessions, and identical numbers of sessions for all '...
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
bms_map_vis.val  = {file img thres scale };
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


