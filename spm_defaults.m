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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner, Andrew Holmes
% $Id: spm_defaults.m 2821 2009-03-03 19:54:19Z guillaume $


global defaults

% Misc
%=======================================================================
defaults.grid     = 0.4;
defaults.cmdline  = 0;
defaults.renderer = 'zbuffer';

% File format specific
%=======================================================================
% Note that defaults.analyze.flip is no longer used.  Specifying the
% left/right handedness of the voxel indices is now done entirely by
% spm_flip_analyze_images.m


% Stats defaults
%=======================================================================
defaults.stats.maxmem   = 2^20;
defaults.stats.maxres   = 64;
defaults.stats.fmri.ufp = 0.001;
defaults.stats.pet.ufp  = 0.05;
defaults.stats.eeg.ufp  = 1;
defaults.stats.topoFDR  = 0;

% fMRI design defaults
%=======================================================================
defaults.stats.fmri.fmri_t    = 16;
defaults.stats.fmri.fmri_t0   = 1;
defaults.stats.fmri.cond.tmod = 0;
defaults.stats.fmri.hpf       = 128;
defaults.stats.fmri.hrf.derivs = [0 0];
defaults.stats.fmri.volt      = 1;
defaults.stats.fmri.global    = 'None';
defaults.stats.fmri.cvi       = 'AR(1)';

% Factorial design defaults
%=======================================================================
defaults.stats.fact.dept     = 0;
defaults.stats.fact.variance = 1;
defaults.stats.fact.t2.gmsca = 0;
defaults.stats.fact.ancova   = 0;
defaults.stats.fact.mcov.iCC = 1;
defaults.stats.fact.iCFI     = 1;
defaults.stats.fact.iCC      = 1;
defaults.stats.fact.athresh  = 100;
defaults.stats.fact.rthresh  = .8;
defaults.stats.fact.imask    = 1;
defaults.stats.fact.gmsca    = 50;
defaults.stats.fact.glonorm  = 1;
defaults.stats.fact.mreg_int = 1;

% Model estimation defaults
%=======================================================================
defaults.stats.est.signal = 'UGL';
defaults.stats.est.ARP    = 3;

% Contrast manager batch defaults
%=======================================================================
defaults.stats.con.delete = 0;

% Results report batch defaults
%=======================================================================
defaults.stats.results.threshtype = 'FWE'; % Threshold type
defaults.stats.results.thresh     = 0.05;  % Threshold value
defaults.stats.results.extent     = 0;     % Spatial extent
defaults.stats.results.maskthresh = 0.05;  % (Uncorrected) Threshold for masking
defaults.stats.results.print      = true;  % Print report to file

% Mask defaults
%=======================================================================
defaults.mask.thresh    = 0.8;

% Filename prefix defaults
%=======================================================================
defaults.slicetiming.prefix     = 'a';
defaults.realign.write.prefix   = 'r';
defaults.coreg.write.prefix     = 'r';
defaults.unwarp.write.prefix    = 'u';
defaults.normalise.write.prefix = 'w';
defaults.smooth.prefix          = 's';

% Realignment defaults
%=======================================================================
defaults.realign.estimate.quality = 0.9;
defaults.realign.estimate.weight = {''};
defaults.realign.estimate.interp = 2;
defaults.realign.estimate.wrap   = [0 0 0];
defaults.realign.estimate.sep    = 4;
defaults.realign.estimate.fwhm   = 5;
defaults.realign.estimate.rtm    = 1;
defaults.realign.write.mask      = 1;
defaults.realign.write.interp    = 4;
defaults.realign.write.wrap      = [0 0 0];
defaults.realign.write.which     = [2 1];

% Unwarp defaults
%=======================================================================
defaults.unwarp.estimate.rtm     = 0;
defaults.unwarp.estimate.fwhm    = 4;
defaults.unwarp.estimate.basfcn  = [12 12];
defaults.unwarp.estimate.regorder= 1;
defaults.unwarp.estimate.regwgt  = 1e5;
defaults.unwarp.estimate.foe     = [4 5];
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
defaults.bias.nbins  = 256; % Number of histogram bins
defaults.bias.reg    = 0.01;    % Regularisation
defaults.bias.cutoff = 30;  % DCT frequency cutoff (mm)

% VBM Preprocessing defaults
%=======================================================================
defaults.preproc.tpm     = cellstr(char(...
    fullfile(spm('Dir'),'tpm','grey.nii'),...
    fullfile(spm('Dir'),'tpm','white.nii'),...
    fullfile(spm('Dir'),'tpm','csf.nii'))); % Prior probability maps
defaults.preproc.ngaus    = [2 2 2 4]';     % Gaussians per class
defaults.preproc.warpreg  = 1;             % Warping Regularisation
defaults.preproc.warpco   = 25;            % Warp Frequency Cutoff
defaults.preproc.biasreg  = 0.0001;        % Bias regularisation
defaults.preproc.biasfwhm = 60;            % Bias FWHM
defaults.preproc.regtype  = 'mni';         % Affine Regularisation
defaults.preproc.samp     = 3;             % Sampling distance
defaults.preproc.output.GM  = [0 0 1];
defaults.preproc.output.WM  = [0 0 1];
defaults.preproc.output.CSF = [0 0 0];
defaults.preproc.output.biascor = 1;
defaults.preproc.output.cleanup = 0;

% Smooth defaults
%=======================================================================
defaults.smooth.fwhm  = [8 8 8];
defaults.smooth.dtype = 0;

% ImCalc defaults
%=======================================================================
defaults.imcalc.output = 'output.img';      % Output Filename
defaults.imcalc.dmtx   = 0;                 % Data Matrix
defaults.imcalc.mask   = 0;                 % Masking
defaults.imcalc.interp = 1;                 % Interpolation
defaults.imcalc.dtype  = spm_type('int16'); % Data Type

% Deformation utility defaults
%=======================================================================
defaults.defs.vox = [NaN NaN NaN];
defaults.defs.bb  = [NaN NaN NaN
                    NaN NaN NaN];
defaults.defs.K   = 6;

% DICOM Import defaults
%=======================================================================
defaults.dicom.root    = 'flat';
defaults.dicom.format  = 'img';
defaults.dicom.icedims = 0;

% MINC Import defaults
%=======================================================================
defaults.minc.dtype    = spm_type('int16');
defaults.minc.ext      = '.img';

% ECAT Import defaults
%=======================================================================
defaults.ecat.ext      = '.img';

% User Interface Defaults
%=======================================================================
defaults.ui.print   = struct('opt',{{'-dpsc2','-append'}},'append',true,'ext','.ps');
defaults.ui.colour1 = [0.8 0.8 1.0];
defaults.ui.colour2 = [0.73 0.78 0.96];
defaults.ui.colour3 = [0.0 0.0 0.0];
defaults.ui.fs      = 14;

% EEG specific
%=======================================================================
defaults.eeg.dtype = 'float';

