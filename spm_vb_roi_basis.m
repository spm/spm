function [F,pm] = spm_vb_roi_basis (VOI_fnames,SPM,bases,model)
% Compare Hemodynamic Basis sets for a cluster of interest
% FORMAT [F,pm] = spm_vb_roi_basis (VOI_fnames,SPM,bases,model)
%
% VOI_fname     VOI filenames eg. VOI_fnames{1}='Test_VOI.mat'
%
% SPM           SPM data structure (this must be loaded in from an 
%               SPM.mat file). If this field is not specified this function
%               wil prompt you for the name of an SPM.mat file
%
% bases         Specifies which basis sets to compare:
%
%               'all'   - the 7 default types (see help spm_get_bf)
%               'fir'   - Finite Impulse Response with variable number of bins
%               'fh'    - Fourier + Hanning window with variable number of bins
%               'user'  - user specified models set by model variable 
%                         (see below). This allows a user-specified set of
%                         models to be compared.
%
%               The default option is 'all'
%
% model         For ith basis set specify
%
%               model(i).name - see help spm_get_bf
%               model(i).sname - short name to be used in results histogram 
%               model(i).order - number of basis functions/number of bins
%               model(i).length - overall window length in seconds
%
%               for i=1..number of models
%
%               This variable only needs to be specified if the bases option
%               is set to 'user'.
%
%               Typical function usages: 
%
%               [F,pm]=spm_vb_roi_basis('Test_VOI.mat');
%               [F,pm]=spm_vb_roi_basis('Test_VOI.mat',SPM);
%               [F,pm]=spm_vb_roi_basis('Test_VOI.mat',SPM,'fir');
%               [F,pm]=spm_vb_roi_basis('Test_VOI.mat',SPM,'user',model);
%
% F             model evidences 
% pm            posterior model probability
%
% See W. Penny et al. (2007). Bayesian Model Comparison of Spatially 
% Regularised General Linear Models. Human Brain Mapping.
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_vb_roi_basis.m 6079 2014-06-30 18:25:37Z spm $


%-Check input parameters
%--------------------------------------------------------------------------
if ~nargin
    [VOI_fnames, sts] = spm_select([1 Inf],'^VOI.*\.mat$','select VOIs');
    if ~sts, F = []; pm = []; return; end
end
VOI_fnames = cellstr(VOI_fnames);
if nargin < 2
    [Pf, sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
    if ~sts, return; end
    swd = spm_file(Pf,'fpath');
    load(fullfile(swd,'SPM.mat'))
end
if nargin < 3, bases = 'all'; end

%-Check parameters
%--------------------------------------------------------------------------
if strcmpi(bases,'user') && nargin ~= 4
    error('Models not specified.');
end

nsess = length(SPM.Sess);
nVOIs = length(VOI_fnames);
if nsess~=nVOIs
    error('SPM design specified for %d sessions and %d VOI files specified.',nsess,nVOIs);
end

%-Load VOIs
%--------------------------------------------------------------------------
for p=1:nVOIs
    load(VOI_fnames{p});
    data_sess(p).Y=xY.y;
end
xyz  = xY.XYZmm;
N    = size(xyz,2);
xyz  = SPM.xVol.M\[xyz;ones(1,size(xyz,2))];
xyz  = xyz(1:3,:);
vxyz = spm_vb_neighbors(xyz',1);

% Set number of AR coefficients
%--------------------------------------------------------------------------
try 
    SPM.PPM.AR_P;
catch
    SPM.PPM.AR_P = 3;
end

% Specify type of prior for regression coefficients
%--------------------------------------------------------------------------
try
    SPM.PPM.priors.W;
catch
    if N==1
        SPM.PPM.priors.W = 'Voxel - Shrinkage';
    else
        SPM.PPM.priors.W = 'Spatial - GMRF';
    end
end

% Specify type of prior for AR coefficients
%--------------------------------------------------------------------------
try
    SPM.PPM.priors.A;
catch
    if N==1
        SPM.PPM.priors.A = 'Voxel - Shrinkage';
    else
        SPM.PPM.priors.A = 'Spatial - GMRF';
    end
end

% Get matrices that will remove low-frequency drifts 
% if high pass filters have been specified
%--------------------------------------------------------------------------
for s=1:nsess
    sess_nScan = length(SPM.xX.K(s).row);
    if size(SPM.xX.K(s).X0,2) > 0
        X0 = SPM.xX.K(s).X0;
        hpf(s).R0 = eye(sess_nScan)-X0*pinv(X0);
    else
        hpf(s).R0 = eye(sess_nScan);
    end    
end

% Set optimisation parameters
%--------------------------------------------------------------------------
try
    SPM.PPM.maxits;
catch
    SPM.PPM.maxits = 16;
end
try
    SPM.PPM.tol;
catch
    SPM.PPM.tol = 0.00001;
end

% Specify basis functions
%--------------------------------------------------------------------------
window_length = 32; % length in seconds

switch bases
    case 'fir'
        window_length = 20; % length in seconds
        bins = [1:1:10];
        for i=1:length(bins),
            model(i).name   = 'Finite Impulse Response';
            model(i).sname  = int2str(bins(i));
            model(i).order  = bins(i); 
            model(i).length = window_length; 
        end
    case 'fh'
        window_length = 20; % length in seconds
        bins = [1:1:10];
        for i=1:length(bins),
            model(i).name   = 'Fourier set (Hanning)';
            model(i).sname  = int2str(bins(i));
            model(i).order  = bins(i); 
            model(i).length = window_length; 
        end
    case 'all'
        model(1).name   = 'hrf';
        model(1).sname  = 'Inf-1';
        model(1).order  = 1;
        model(1).length = window_length; 
        
        model(2).name   = 'hrf (with time derivative)';
        model(2).sname  = 'Inf-2';
        model(2).order  = 2;
        model(2).length = window_length;
        
        model(3).name   = 'hrf (with time and dispersion derivatives)';
        model(3).sname  = 'Inf-3';
        model(3).order  = 3;  
        model(3).length = window_length;                
        
        model(4).name   = 'Fourier set';
        model(4).sname  = 'F';
        model(4).order  = 5;
        model(4).length = window_length;
        
        model(5).name   = 'Fourier set (Hanning)';
        model(5).sname  = 'FH';
        model(5).order  = 5;
        model(5).length = window_length;
        
        model(6).name   = 'Gamma functions';
        model(6).sname  = 'Gamm3';
        model(6).order  = 3;  
        model(6).length = window_length;                
        
        model(7).name   = 'Finite Impulse Response';
        model(7).sname  = 'FIR';
        model(7).order  = 10; 
        model(7).length = window_length;        
    case 'user'
        disp('Using user-specified models.');
    otherwise
        error('Unknown bases option.');
end

% Fit models
%==========================================================================
for m=1:length(model)
    fprintf('Model %d (%s)\n',m,model(m).name);                             %-#
    
    % Switch basis functions
    %----------------------------------------------------------------------
    SPM.xBF.name   = model(m).name;
    SPM.xBF.order  = model(m).order;
    SPM.xBF.length = model(m).length;
    SPM.xBF        = spm_get_bf(SPM.xBF);
    
    SPM = spm_fmri_design(SPM,0);  % 0 to override writing of SPM.mat file
    
    Ev = 0;
    % Loop over sessions
    %======================================================================
    for s=1:nsess
        fprintf('Session %d\n',s);                                      %-#
        
        % Design matrix
        %------------------------------------------------------------------
        X = SPM.xX.X(SPM.Sess(s).row,SPM.Sess(s).col);
        X = [X ones(length(SPM.Sess(s).row),1)]; % Add one constant 
        
        slice = spm_vb_init_volume (X,SPM.PPM.AR_P);
        slice = spm_vb_set_priors(slice,SPM.PPM.priors,vxyz);
        
        slice.maxits        = SPM.PPM.maxits;
        slice.tol           = SPM.PPM.tol;
        slice.compute_det_D = 1;
        slice.verbose       = 1;
        slice.update_w      = 1;
        slice.update_lambda = 1;
        slice.update_F      = 1;
    
        % Filter data to remove low frequencies
        %------------------------------------------------------------------
        R0Y = hpf(s).R0 * data_sess(s).Y;
        
        % Variational Bayes for GLM-AR modelling
        %------------------------------------------------------------------
        slice = spm_vb_glmar(R0Y,slice);
        
        % Add up evidence from different sessions
        %------------------------------------------------------------------
        Ev = Ev + slice.F;
    end
    
    F(m) = Ev;
end

% Display results
%==========================================================================
Fm    = F - mean(F);
Fnorm = F - min(F);
pm    = exp(Fm)./sum(exp(Fm));
if any(isnan(pm))
    [pF,pi] = max(F);
    pm      = zeros(size(pi));
    pm(pi)  = 1;
end

figure;
bar(pm);
set(gca,'XTickLabel',{model.sname});
ylabel('Posterior model probability');
set(gca,'FontSize',12);

figure;
bar(Fnorm);
set(gca,'XTickLabel',{model.sname});
ylabel('Log Evidence');
set(gca,'FontSize',12);
