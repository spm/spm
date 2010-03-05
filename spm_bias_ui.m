function spm_bias_ui(P)
% Non-uniformity correct images.
%
% The objective function is related to minimising the entropy of
% the image histogram, but is modified slightly.
% This fixes the problem with the SPM99 non-uniformity correction
% algorithm, which tends to try to reduce the image intensities. As
% the field was constrainded to have an average value of one, then
% this caused the field to bend upwards in regions not included in
% computations of image non-uniformity.
%
%_______________________________________________________________________
% Ref:
% J Ashburner. 2002. "Another MRI Bias Correction Approach" [abstract].
% Presented at the 8th International Conference on Functional Mapping of
% the Human Brain, June 2-6, 2002, Sendai, Japan. Available on CD-Rom
% in NeuroImage, Vol. 16, No. 2.
%
%_______________________________________________________________________
%
%                        The Prompts Explained
%_______________________________________________________________________
%
% 'Scans to correct' - self explanatory
%
%_______________________________________________________________________
%
%                           Defaults Options
%_______________________________________________________________________
%[   things in square brackets indicate corresponding defaults field   ]
%
% 'Number of histogram bins?'
% The probability density of the image intensity is represented by a
% histogram. The optimum number of bins depends on the number of voxels
% in the image.  More voxels allows a more detailed representation.
% Another factor is any potential aliasing effect due to there being a
% discrete number of different intensities in the image.  Fewer bins
% should be used in this case.
% [defaults.bias.nbins]
%
% 'Regularisation?'
% The importance of smoothness for the estimated bias field. Without
% any regularisation, the algorithm will attempt to correct for
% different grey levels arising from different tissue types, rather than
% just correcting bias artifact.
% Bias correction uses a Bayesian framework (again) to model intensity
% inhomogeneities in the image(s).  The variance associated with each
% tissue class is assumed to be multiplicative (with the
% inhomogeneities).  The low frequency intensity variability is
% modelled by a linear combination of three dimensional DCT basis
% functions (again), using a fast algorithm (again) to generate the
% curvature matrix.  The regularization is based upon minimizing the
% integral of square of the fourth derivatives of the modulation field
% (the integral of the squares of the first and second derivs give the
% membrane and bending energies respectively).
% [defaults.bias.reg]
%
% 'Cutoff?'
% Cutoff of DCT bases.  Only DCT bases of periods longer than the
% cutoff are used to describe the warps. The number used will
% depend on the cutoff and the field of view of the image.
% [defaults.bias.cutoff]
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_bias_ui.m 3756 2010-03-05 18:43:37Z guillaume $


global defaults

if nargin==1 && strcmpi(P,'defaults');
    defaults.bias = edit_defaults(defaults.bias);
    return;
end;
bias_ui(defaults.bias);
return;
%=======================================================================

%=======================================================================
function bias_ui(flags)
% User interface for nonuniformity correction
spm('FnBanner',mfilename,'$Rev: 3756 $');
[Finter,unused,CmdLine] = spm('FnUIsetup','Flatten');
spm_help('!ContextHelp',mfilename);
PP = spm_select(Inf, 'image', 'Scans to correct');
spm('Pointer','Watch');
for i=1:size(PP,1),
    spm('FigName',['Flatten: working on scan ' num2str(i)],Finter,CmdLine);
    drawnow;
    P              = deblank(PP(i,:));
    T              = spm_bias_estimate(P,flags);
    [pth,nm,xt,vr] = spm_fileparts(P);
    S              = fullfile(pth,['bias_' nm '.mat']);
    %S             = ['bias_' nm '.mat'];
    spm_bias_apply(P,S);
end;
if 0,
fg = spm_figure('FindWin','Interactive');
if ~isempty(fg), spm_figure('Clear',fg); end;
end
spm('FigName','Flatten: done',Finter,CmdLine);
spm('Pointer');
return;
%=======================================================================

%=======================================================================
function flags = edit_defaults(flags)

nb  = [32 64 128 256 512 1024 2048];
tmp = find(nb == flags.nbins);
if isempty(tmp), tmp = 6; end;
flags.nbins = spm_input('Number of histogram bins?','+1','m',...
        ['  32 bins |  64 bins| 128 bins| 256 bins| 512 bins|1024 bins|2048 bins'],...
         nb, tmp);

rg  = [0 0.00001 0.0001 0.001 0.01 0.1 1.0 10];
tmp = find(rg == flags.reg);
if isempty(tmp), tmp = 4; end;
flags.reg = spm_input('Regularisation?','+1','m',...
        ['no regularisation (0)|extremely light regularisation (0.00001)|'...
     'very light regularisation (0.0001)|light regularisation (0.001)|',...
     'medium regularisation (0.01)|heavy regularisation (0.1)|'...
     'very heavy regularisation (1)|extremely heavy regularisation (10)'],...
         rg, tmp);

co  = [20 25 30 35 40 45 50 60 70 80 90 100];
tmp = find(co == flags.cutoff);
if isempty(tmp), tmp = 4; end;
flags.cutoff =  spm_input('Cutoff?','+1','m',...
    [' 20mm cutoff| 25mm cutoff| 30mm cutoff| 35mm cutoff| 40mm cutoff|'...
     ' 45mm cutoff| 50mm cutoff| 60mm cutoff| 70mm cutoff| 80mm cutoff|'...
     ' 90mm cutoff|100mm cutoff'],...
    co, tmp);

return;
%=======================================================================
