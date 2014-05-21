function bms = spm_cfg_bms_map
% Configuration file for BMS interface
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Maria Joao Rosa
% $Id: spm_cfg_bms_map.m 6004 2014-05-21 14:24:14Z guillaume $

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
% mod_map Log evidence maps
%--------------------------------------------------------------------------
mod_map         = cfg_files;
mod_map.tag     = 'mod_map';
mod_map.name    = 'Models';
mod_map.help    = {['Specify the log. evidence map for each model. '...
                    'Log-evidence maps should be specified '...
                    'in the same order for each subject and session.']};
mod_map.filter  = 'image';
mod_map.ufilter = '.*';
mod_map.num     = [1 Inf];

%--------------------------------------------------------------------------
% sess_map Sessions (Maps)
%--------------------------------------------------------------------------
sess_map      = cfg_branch;
sess_map.tag  = 'sess_map';
sess_map.name = 'Session';
sess_map.val  = {mod_map };

%--------------------------------------------------------------------------
% subj_dcm Subject (Maps)
%--------------------------------------------------------------------------
subj_map         = cfg_repeat;
subj_map.tag     = 'subj_map';
subj_map.name    = 'Subject';
subj_map.values  = {sess_map };

%--------------------------------------------------------------------------
% map Data
%--------------------------------------------------------------------------
map         = cfg_repeat;
map.tag     = 'map';
map.name    = 'Data';
map.help    = {['Select the log. evidence maps for each '...
               'model, session and subject.']}';
map.values  = {subj_map };
map.num     = [1 Inf];

%--------------------------------------------------------------------------
% mod_name Name
%--------------------------------------------------------------------------
mod_name         = cfg_entry;
mod_name.tag     = 'mod_name';
mod_name.name    = 'Name';
mod_name.help    = {'Specify name for each model (optional).'};
mod_name.strtype = 's';
mod_name.num     = [0 Inf];
mod_name.val     = {''};

%--------------------------------------------------------------------------
% name_mod Name models
%--------------------------------------------------------------------------
name_mod         = cfg_repeat;
name_mod.tag     = 'name_mod';
name_mod.name    = 'Name models';
name_mod.help    = {'Specify name for each model (optional).'}';
name_mod.values  = {mod_name };
name_mod.num     = [0 Inf];

%--------------------------------------------------------------------------
% method_maps Inference Method (maps)
%--------------------------------------------------------------------------
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
% out_file Output files
%--------------------------------------------------------------------------
out_file         = cfg_menu;
out_file.tag     = 'out_file';
out_file.name    = 'Output files (RFX)';
out_file.help    = {['Specify which output files to save (only valid for'...
                     'RFX analyses). ']...
                     ''...
                    ['Default option (and faster option): '...
                     'PPM = xppm.<ext> (Expected Posterior Probability Maps) '...
                     'for each model ie. posterior mean.']...
                     ''...
                    ['Second option: PPM + EPM = xppm.<ext> + '...
                     'epm.<ext> (Expected Posterior Probability '...
                     'Maps + Exceedance Probability Maps) for each model.']...
                     ''...
                    ['Third option: PPM + EPM + Alpha = xppm.<ext> + '...
                     'epm.<ext> + alpha.<ext> (PPM, EPM and Map of Dirichlet '...
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

%--------------------------------------------------------------------------
% mask Mask Image
%--------------------------------------------------------------------------
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

%--------------------------------------------------------------------------
% nsamp Number of samples
%--------------------------------------------------------------------------
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

%--------------------------------------------------------------------------
% file BMS.mat
%--------------------------------------------------------------------------
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

%--------------------------------------------------------------------------
% img Map to display
%--------------------------------------------------------------------------
img         = cfg_files;
img.tag     = 'img';
img.name    = 'Map to display';
img.help    = {['Specify map obtained from BMS Maps '...
               '(optional). Leave field empty to work on serial mode.']};
img.filter  = 'image';
img.ufilter = '.*';
img.val     = {{''}};
img.num     = [0 1];

%--------------------------------------------------------------------------
% thres Probability Threshold
%--------------------------------------------------------------------------
thres         = cfg_entry;
thres.tag     = 'thres';
thres.name    = 'Probability threshold';
thres.help    = {['Specify the probability threshold to apply to the '...
                 'image (optional). Leave field empty to work on '...
                 'serial mode.']};                 
thres.strtype = 'r';
thres.num     = [0 Inf];
thres.val     = {[]};

%--------------------------------------------------------------------------
% k Extent threshold
%--------------------------------------------------------------------------
k         = cfg_entry;
k.tag     = 'k';
k.name    = 'Extent threshold';
k.help    = {['Specify extent threshold (minimum number of voxels '...
                 'per cluster).']};                 
k.strtype = 'w';
k.num     = [0 Inf];
k.val     = {[]};

%--------------------------------------------------------------------------
% scale Map Scale
%--------------------------------------------------------------------------
scale         = cfg_menu;
scale.tag     = 'scale';
scale.name    = 'Map scale';
scale.help    = {['Specify scale to display maps (optional). Default: '...
                 'empty field to work on serial mode. Other options: '...
                 '''None'' will display image with original scale and '...
                 '''Log-odds'' will display image in a log-odds '...
                 ' scale (in this case image should be a '...
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

%--------------------------------------------------------------------------
% bms_map_inf BMS: Maps (Inference), output is BMS map 
%--------------------------------------------------------------------------
bms_map_inf      = cfg_exbranch;
bms_map_inf.tag  = 'inference';
bms_map_inf.name = 'BMS: Maps (Inference)';
bms_map_inf.val  = {dir map name_mod method_maps out_file mask nsamp };
bms_map_inf.help = {'Bayesian Model Selection for Log-Evidence Maps. '...
    ''...
    ['Input: log-evidence maps for each model, session and '...
    'subject. Note that there must be identical numbers of models for '...
    'all sessions, and identical numbers of sessions for all '...
    'subjects.']...
    ''...
    ['Output: For the fixed effects analysis, posterior probability maps '...
    'are created for each model. '...
    'For the random effects analysis, expected posterior probability '...
    'and exceedance probability (i.e. the probability that this model '...
    'is more likely than any other model) maps are created for each '...
    'model. If there are multiple sessions per subject, the random '...
    'effects analysis operates on the subject-specific sums of log '...
    'evidences across sessions. In addition, a BMS.mat file will be save '...
    'in the specified directory for both methods']};
bms_map_inf.prog = @spm_run_bms_map;
bms_map_inf.vout = @vout;

%--------------------------------------------------------------------------
% bms_map_vis BMS: Maps (Results), visualisation of BMS Maps results
%--------------------------------------------------------------------------
bms_map_vis      = cfg_exbranch;
bms_map_vis.tag  = 'results';
bms_map_vis.name = 'BMS: Maps (Results)';
bms_map_vis.val  = {file img thres k scale};
bms_map_vis.help = {['Bayesian Model Selection Maps (Results). '...
                    'Show results from BMS Maps (Inference).']};
bms_map_vis.prog = @spm_run_bms_vis;

%--------------------------------------------------------------------------
% bms Bayesian Model Selection
%--------------------------------------------------------------------------
bms         = cfg_choice;
bms.tag     = 'bms_map';
bms.name    = 'Bayesian Model Selection';
bms.help    = {['Bayesian Model Selection for group studies (fixed '...
               'effects and random effects analysis).']};
bms.values  = { bms_map_inf bms_map_vis };


%==========================================================================
function dep = vout(varargin)
% Output file names will be saved in a struct with field .files
dep(1)            = cfg_dep;
dep(1).sname      = 'BMS.mat File';
dep(1).src_output = substruct('.','files');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
