
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
% %W% John Ashburner, Andrew Holmes %E%

global defaults

% Misc
%=======================================================================
defaults.grid     = 0.4;
defaults.cmdline  = 0;
defaults.logfile  = '';
defaults.printstr =  [spm_figure('DefPrintCmd'),'spm2.ps'];

% File format specific
%=======================================================================
defaults.analyze.multivol = 0;
defaults.analyze.flip     = 1; % <<= Very important.  Relates to L/R

% Stats defaults
%=======================================================================
defaults.stats.maxmem   = 2^20;
defaults.stats.maxres   = 64;
defaults.stats.fmri.ufp = 0.001;
defaults.stats.fmri.t   = 16;
defaults.stats.fmri.t0  = 1;
defaults.stats.pet.ufp  = 0.05;

% Mask defaults
%=======================================================================
defaults.mask.thresh    = 0.8;

% Realignment defaults
%=======================================================================
defaults.realign.estimate.quality = 0.75;
defaults.realign.estimate.weight = 0;
defaults.realign.estimate.interp = 2;
defaults.realign.estimate.wrap   = [0 0 0];
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
defaults.segment.estimate.priors = str2mat(...
	fullfile(spm('Dir'),'apriori','gray.mnc'),...
	fullfile(spm('Dir'),'apriori','white.mnc'),...
	fullfile(spm('Dir'),'apriori','csf.mnc'));
defaults.segment.estimate.reg    = 0.01;
defaults.segment.estimate.cutoff = 30;
defaults.segment.estimate.samp   = 3;
defaults.segment.estimate.bb     =  [[-88 88]' [-122 86]' [-60 95]'];
defaults.segment.estimate.affreg.smosrc = 8;
defaults.segment.estimate.affreg.regtype = 'mni';
%defaults.segment.estimate.affreg.weight = fullfile(spm('Dir'),'apriori','brainmask.mnc'); 
defaults.segment.estimate.affreg.weight = '';
defaults.segment.write.cleanup   = 1;
defaults.segment.write.wrt_cor   = 1;

% Bias field estimation defaults
%=======================================================================
defaults.bias.nbins  = 256;	% Number of histogram bins
defaults.bias.reg    = 0.01;	% Regularisation
defaults.bias.cutoff = 30;	% DCT frequency cutoff (mm)

