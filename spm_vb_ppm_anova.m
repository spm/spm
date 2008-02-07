function spm_vb_ppm_anova(SPM)
% Bayesian ANOVA using model comparison
% FORMAT spm_vb_ppm_anova(SPM)
%
% SPM    -  Data structure corresponding to a full model (ie. one
%           containing all experimental conditions).
%           
% This function creates images of differences in log evidence
% which characterise the average effect, main effects and interactions
% in a factorial design. 
%
% The factorial design is specified in SPM.factor. For a one-way ANOVA 
% the images 
%
%   avg_effect.img
%   main_effect.img
%
% are produced. For a two-way ANOVA the following images are produced
%
%   avg_effect.img
%   main_effect_'factor1'.img
%   main_effect_'factor2'.img
%   interaction.img
%
% These images can then be thresholded. For example a threshold of 4.6 
% corresponds to a posterior effect probability of [exp(4.6)] = 0.999. 
% See paper VB4 for more details.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_vb_ppm_anova.m 1143 2008-02-07 19:33:33Z spm $

    
disp('Warning: spm_vb_ppm_anova only works for single session data.');

model = spm_vb_models(SPM,SPM.factor);

analysis_dir = pwd;
for m=1:length(model)-1,
    model_subdir = ['model_',int2str(m)];
    mkdir(analysis_dir,model_subdir);
    SPM.swd = fullfile(analysis_dir,model_subdir);
    SPM.Sess(1).U = model(m).U;
    SPM.Sess(1).U = spm_get_ons(SPM,1);
    SPM = spm_fMRI_design(SPM,0); % 0 = don't save SPM.mat
    SPM.PPM.update_F = 1; % Compute evidence for each model
    SPM.PPM.compute_det_D = 1; 
    spm_spm_vb(SPM);
end

% Compute differences in contributions to log-evidence images
% to assess main effects and interactions
nf = length(SPM.factor);
if nf==1
    % For a single factor
    
    % Average effect
    image1 = fullfile(analysis_dir, 'model_1','LogEv.img');
    image2 = fullfile(analysis_dir, 'model_2','LogEv.img');
    imout  = fullfile(analysis_dir, 'avg_effect.img');
    img_subtract(image1,image2,imout);
    
    % Main effect of factor
    image1 = fullfile(analysis_dir, 'model_2','LogEv.img');
    image2 = fullfile(analysis_dir, 'LogEv.img');
    imout  = fullfile(analysis_dir, 'main_effect.img');
    img_subtract(image1,image2,imout);
    
elseif nf==2
    % For two factors
    
    % Average effect
    image1 = fullfile(analysis_dir, 'model_1','LogEv.img');
    image2 = fullfile(analysis_dir, 'model_2','LogEv.img');
    imout  = fullfile(analysis_dir, 'avg_effect.img');
    img_subtract(image1,image2,imout);
    
    % Main effect of factor 1
    image1 = fullfile(analysis_dir, 'model_2','LogEv.img');
    image2 = fullfile(analysis_dir, 'model_3','LogEv.img');
    imout  = fullfile(analysis_dir, ['main_effect_',SPM.factor(1).name,'.img']);
    img_subtract(image1,image2,imout);
    
    % Main effect of factor 2
    image1 = fullfile(analysis_dir, 'model_2','LogEv.img');
    image2 = fullfile(analysis_dir, 'model_4','LogEv.img');
    imout  = fullfile(analysis_dir, ['main_effect_',SPM.factor(2).name,'.img']);
    img_subtract(image1,image2,imout);
    
    % Interaction
    image1 = fullfile(analysis_dir, 'model_5','LogEv.img');
    image2 = fullfile(analysis_dir, 'LogEv.img');
    imout  = fullfile(analysis_dir, 'interaction.img');
    img_subtract(image1,image2,imout);
    
end
    
%-----------------------------------------------------------------------
function img_subtract(image1,image2,image_out)
% Subtract image 1 from image 2 and write to image out
% Note: parameters are names of files

Vi = spm_vol(strvcat(image1,image2));
Vo = struct(...
    'fname',   image_out,...
    'dim',     [Vi(1).dim(1:3)],...
    'dt',      [spm_type('float32') spm_platform('bigend')],...
    'mat',     Vi(1).mat,...
    'descrip', 'Difference in Log Evidence');
f = 'i2-i1';
flags = {0,0,1};
Vo = spm_imcalc(Vi,Vo,f,flags);
