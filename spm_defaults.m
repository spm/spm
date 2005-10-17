function spm_defaults
% Sets the defaults which are used by SPM
%
% FORMAT spm_defaults
%_______________________________________________________________________
%
% This file is intended to be customised for the site.
% Individual users can make copies which can be stored in their own
% matlab subdirectories. If ~/matlab is ahead of the SPM directory
% in the MATLABPATH, then the users own personal defaults are used.
%
% Care must be taken when modifying this file
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner, Andrew Holmes
% $Id: spm_defaults.m 256 2005-10-17 18:57:24Z guillaume $


global defaults

% Misc
%=======================================================================
defaults.grid     = 0.4;
defaults.cmdline  = 0;

% File format specific
%=======================================================================
defaults.analyze.multivol = 0;
defaults.analyze.flip     = 1; % <<= Very important.  Relates to L/R

% Stats defaults
%=======================================================================
defaults.stats.maxmem   = 2^20;
defaults.stats.maxres   = 64;
defaults.stats.fmri.ufp = 0.001;
defaults.stats.pet.ufp  = 0.05;
defaults.stats.eeg.ufp  = 1;

% Mask defaults
%=======================================================================
defaults.mask.thresh    = 0.8;

% Realignment defaults
%=======================================================================
defaults.realign.estimate.quality = 0.9;
defaults.realign.estimate.weight = 0;
defaults.realign.estimate.interp = 2;
defaults.realign.estimate.wrap   = [0 0 0];
defaults.realign.estimate.sep    = 4;
defaults.realign.estimate.fwhm   = 5;
defaults.realign.estimate.rtm    = 1;
defaults.realign.write.mask      = 1;
defaults.realign.write.interp    = 4;
defaults.realign.write.wrap      = [0 0 0];

% Unwarp defaults
%=======================================================================
defaults.unwarp.estimate.fwhm    = 4;
defaults.unwarp.estimate.basfcn  = [12 12];
defaults.unwarp.estimate.regorder= 1;
defaults.unwarp.estimate.regwgt  = 1e5;
defaults.unwarp.estimate.soe     = 1;
defaults.unwarp.estimate.rem     = 1;
defaults.unwarp.estimate.jm      = 0;
defaults.unwarp.estimate.noi     = 5;
defaults.unwarp.estimate.expround= 'Average';
%
% Unwarp uses defaults.realign.write 
% defaults for writing.
%

% Coregistration defaults
%=======================================================================
defaults.coreg.estimate.cost_fun = 'nmi';
defaults.coreg.estimate.sep      = [4 2];
defaults.coreg.estimate.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
defaults.coreg.estimate.fwhm     = [7 7];
defaults.coreg.write.interp      = 1;
defaults.coreg.write.wrap        = [0 0 0];
defaults.coreg.write.mask        = 0;

% Spatial Normalisation defaults
%=======================================================================
defaults.normalise.estimate.smosrc  = 8;
defaults.normalise.estimate.smoref  = 0;
defaults.normalise.estimate.regtype = 'mni';
defaults.normalise.estimate.weight  = '';
defaults.normalise.estimate.cutoff  = 25;
defaults.normalise.estimate.nits    = 16;
defaults.normalise.estimate.reg     = 1;
defaults.normalise.estimate.wtsrc   = 0;
defaults.normalise.write.preserve   = 0;
defaults.normalise.write.bb         = [[-78 -112 -50];[78 76 85]];
defaults.normalise.write.vox        = [2 2 2];
defaults.normalise.write.interp     = 1;
defaults.normalise.write.wrap       = [0 0 0];

% Segmentation defaults
%=======================================================================
defaults.segment.estimate.priors = char(...
	fullfile(spm('Dir'),'apriori','grey.nii'),...
	fullfile(spm('Dir'),'apriori','white.nii'),...
	fullfile(spm('Dir'),'apriori','csf.nii'));
defaults.segment.estimate.reg    = 0.01;
defaults.segment.estimate.cutoff = 30;
defaults.segment.estimate.samp   = 3;
defaults.segment.estimate.bb     =  [[-88 88]' [-122 86]' [-60 95]'];
defaults.segment.estimate.affreg.smosrc = 8;
defaults.segment.estimate.affreg.regtype = 'mni';
%defaults.segment.estimate.affreg.weight = fullfile(spm('Dir'),'apriori','brainmask.nii'); 
defaults.segment.estimate.affreg.weight = '';
defaults.segment.write.cleanup   = 1;
defaults.segment.write.wrt_cor   = 1;
defaults.segment.write.wrt_brV   = 1;

% Bias field estimation defaults (UNUSED)
%=======================================================================
defaults.bias.nbins  = 256;	% Number of histogram bins
defaults.bias.reg    = 0.01;	% Regularisation
defaults.bias.cutoff = 30;	% DCT frequency cutoff (mm)

% VBM Preprocessing defaults
%=======================================================================
defaults.preproc.tpm     = char(...
    fullfile(spm('Dir'),'tpm','grey.nii'),...
    fullfile(spm('Dir'),'tpm','white.nii'),...
    fullfile(spm('Dir'),'tpm','csf.nii')); % Prior probability maps
defaults.preproc.ngaus    = [2 2 2 4];     % Gaussians per class
defaults.preproc.warpreg  = 1;             % Warping Regularisation
defaults.preproc.warpco   = 25;            % Warp Frequency Cutoff
defaults.preproc.biasreg  = 0.0001;        % Bias regularisation
defaults.preproc.biasfwhm = 60;            % Bias FWHM
defaults.preproc.regtype  = 'mni';         % Affine Regularisation
defaults.preproc.samp     = 3;             % Sampling distance

% ImCalc defaults
%=======================================================================
defaults.imcalc.dmtx   = 0;                 % Data Matrix
defaults.imcalc.mask   = 0;                 % Masking
defaults.imcalc.interp = 1;                 % Interpolation
defaults.imcalc.dtype  = spm_type('int16'); % Data Type

% User Interface Defaults
%=======================================================================
defaults.ui.print   = struct('opt',{{'-dpsc2','-append'}},'append',true,'ext','.ps');
defaults.ui.colour1 = [0.8 0.8 1.0];
defaults.ui.colour2 = [1.0 1.0 0.8];
defaults.ui.colour3 = [0.0 0.0 0.0];
defaults.ui.fs      = 14;

% EEG specific
%=======================================================================
defaults.eeg.dtype = 'float';

