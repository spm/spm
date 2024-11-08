%+-------------------------------------------------------------------------
%| An SPM GUI-based toolbox to perform pTFCE
%| (probabilistic Threshold-free Cluster Enhancement)
%+-------------------------------------------------------------------------
%|
%| The threshold-free cluster enhancement (TFCE) approach integrates cluster
%| information into voxel-wise statistical inference to enhance detectability
%| of neuroimaging signal. Despite the significantly increased sensitivity,
%| the application of TFCE is limited by several factors: (i) generalization
%| to data structures, like brain network connectivity data is not trivial,
%| (ii) TFCE values are in an arbitrary unit, therefore, P-values can only be
%| obtained by a computationally demanding permutation-test.
%| Here, we introduce a probabilistic approach for TFCE (pTFCE), that gives
%| a simple general framework for topology-based belief boosting.
%| The core of pTFCE is a conditional probability, calculated based on Bayes'
%| rule, from the probability of voxel intensity and the threshold-wise
%| likelihood function of the measured cluster size. We provide an estimation
%| of these distributions based on Gaussian Random Field (GRF) theory.
%| The conditional probabilities are then aggregated across cluster-forming
%| thresholds by a novel incremental aggregation method. Our approach is
%| validated on simulated and real fMRI data.
%| pTFCE is shown to be more robust to various ground truth shapes and
%| provides a stricter control over cluster "leaking" than TFCE and, in the
%| most realistic cases, further improves its sensitivity. Correction for
%| multiple comparison can be trivially performed on the enhanced P-values,
%| without the need for permutation testing, thus pTFCE is well-suitable for
%| the improvement of statistical inference in any neuroimaging workflow.
%|
%| This matlab implementation is a port of the pTFCE R-package (v0.0.4) and
%| validated to that.
%+-----------------------------------------------------------------------+
%| PLEASE CITE                                                           |
%+-----------------------------------------------------------------------+
%| T. Spisák, Z. Spisák, M. Zunhammer, U. Bingel, S. Smith, T. Nichols, T.
%| Kincses, Probabilistic TFCE: a generalized combination of cluster size
%| and voxel intensity to increase statistical power, Neuroimage 185:12-26,
%| 2019.
%+-----------------------------------------------------------------------+
%| WEBSITE                                                               |
%| https://github.com/spisakt/pTFCE                                      |
%+-----------------------------------------------------------------------+
%| Version 0.1.3                                                     |
%+-----------------------------------------------------------------------+
%| University Hospital Essen                                             |
%| Department of Neurology                                               |
%| Bingel Lab                                                            |
%+-----------------------------------------------------------------------+
%|                                                                       |
%| Please acknowledge in publications.                                   |
%|                                                                       |
%| Kindly,                                                               |
%|                                                                       |
%| Tamas Spisak                                                          |
%| tamas.spisak@uk-essen.de                                              |
%+-----------------------------------------------------------------------+



function varargout = gui_pTFCE(varargin)

SCCSid  = '0.5';

global BCH; %- used as a flag to know if we are in batch mode or not.

% First read SPM.mat
[spmmatfile, sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
swd = spm_file(spmmatfile,'fpath');

try
    load(fullfile(swd,'SPM.mat'));
catch
    error(['Cannot read ' fullfile(swd,'SPM.mat')]);
end
SPM.swd = swd;
        
if isempty(SPM) 
    spm('alert','Invalid SPM.mat!','pTFCE',[],0);
	varargout = {[]};
	return;
end

%-Change directory so that relative filenames are valid
%--------------------------------------------------------------------------
cd(SPM.swd);

%-Check the model has been estimated
%--------------------------------------------------------------------------
try
    SPM.xVol.S;
catch
    spm('alert*',{'This model has not been estimated.','',...
        fullfile(swd,'SPM.mat')}, mfilename, [], ~spm('CmdLine'));
    SPM = []; xSPM = [];
    return
end

xX   = SPM.xX;                      %-Design definition structure
XYZ  = SPM.xVol.XYZ;                %-XYZ coordinates
S    = SPM.xVol.S;                  %-search Volume {voxels}
R    = SPM.xVol.R(4);               %-search Volume {resels}
M    = SPM.xVol.M(1:3,1:3);         %-voxels to mm matrix
VOX  = sqrt(diag(M'*M))';           %-voxel dimensions

%-Get contrasts
%--------------------------------------------------------------------------
try, xCon = SPM.xCon; catch, xCon = {}; end

try
    Ic        = xSPM.Ic;
catch
    [Ic,xCon] = spm_conman(SPM,'T&F',1,...
                           '    Select contrasts...',' for conjunction',0);
end
if isempty(xCon)
    % figure out whether new contrasts were defined, but not selected
    % do this by comparing length of SPM.xCon to xCon, remember added
    % indices to run spm_contrasts on them as well
    try
        noxCon = numel(SPM.xCon);
    catch
        noxCon = 0;
    end
    IcAdd = (noxCon+1):numel(xCon);
else
    IcAdd = [];
end

nc        = length(Ic);  % Number of contrasts

if nc > 1
    error(['Only one contrast can be selected!']);
end

% clone contrast
Iptfce=length(SPM.xCon)+1;
SPM.xCon(Iptfce)=SPM.xCon(Ic);

% extract statistical iamge from the contrast

STAT = SPM.xCon(Iptfce).STAT;
DOF = SPM.xX.erdf;
fname = SPM.xCon(Iptfce).Vspm.fname;
maskfname = SPM.VM.fname;
descrip = SPM.xCon(Iptfce).Vspm.descrip;

% load image
if ~spm_existfile(fname)
        error('File not found: %s',fname);
end
fprintf('Image to be enhanced: %s\n', fname);

V = spm_vol(fname);
M = spm_vol(maskfname);

img = spm_read_vols(V);
mask = spm_read_vols(M);

% convert it to Z-score (from T or F)

if strcmp(STAT,'T')
    imgZ = img_t2z(img, DOF, 1);
elseif strcmp(STAT,'F')
    %error('F-tests not yet supported! Please contact the authors');
    % Note that this is a hack!!
    % ToDo: find a better way to get DOF for F-tests
    % must be very simple...
    %x=strsplit('SPM{F_[1.0,11.0]} - contrast 7: f-test', '[')
    %xx=strsplit(x{2}, ',')
    %DOF1=str2num(xx{1})
    %xxx=strsplit(xx{2}, ']')
    %DOF2=str2num(xxx{1})
    DOF2=DOF
    DOF1=SPM.xCon(Iptfce).eidf
    imgZ = img_f2z(img, DOF1, DOF2, 1);
else
    error('Statistical parameter type not supported: %s',STAT);
end

% developer: write out unenhanced Z-score image for validation
%V_ptfce = V;
%V_ptfce.fname = 'test_unenhnaced_Z.nii';
%spm_write_vol(V_ptfce,imgZ);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%developer: set R and S here explicitly for validation
% these are the values for the SPM tutorial dataset
% as calcuklated by the R-package
%R=658.4739
%S=51560
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% apply pTFCE
[pTFCE_Z, ptfce_p] = pTFCE(imgZ, mask, R, S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% developer: write out Z-score image for validation
V_ptfce = V;
V_ptfce.fname = 'test_pTFCE_p.nii';
spm_write_vol(V_ptfce,ptfce_p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert back to original statistics (T or F)

if strcmp(STAT,'T')
    img_out = tinv(1-ptfce_p, DOF);
elseif strcmp(STAT,'F')
    img_out = finv(1-ptfce_p, DOF1, DOF2);
end

% save volume
V_ptfce=V;
V_ptfce.descrip=strcat('pTFCE_', V.descrip);
[path, name, ext]=fileparts(V.fname);
if path
    V_ptfce.fname = strcat(path, filesep, 'pTFCE_', name, ext);
else
    V_ptfce.fname = strcat('pTFCE_', name, ext);
end

spm_write_vol(V_ptfce,img_out);

% update and save SPM.mat

SPM.xCon(Iptfce).Vspm.fname = V_ptfce.fname;
SPM.xCon(Iptfce).Vspm.descrip = V_ptfce.descrip;
SPM.xCon(Iptfce).name = strcat(SPM.xCon(Iptfce).name, ' (pTFCE)');


save('SPM.mat','SPM', spm_get_defaults('mat.format'));
fprintf('%s%30s\n','Saving SPM.mat','...done')   

