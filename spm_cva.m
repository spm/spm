function [CVA] = spm_cva(xSPM,SPM,hReg,CVA)
% VOI extraction of adjusted data and CVA
% FORMAT [CVA] = spm_cva(xSPM,SPM,hReg,CVA)
%
% xSPM   - structure containing specific SPM details
%     xSPM.Ic  - indice of contrast (in SPM.xCon)
% SPM    - structure containing generic analysis details
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% CVA.contrast =  contrast name
% CVA.name     =  CVA name
% CVA.c        =  contrast weights
% CVA.X        =  contrast subspace
% CVA.Y        =  whitened and adjusted data
% CVA.X0       =  null space of contrast
% 
% CVA.XYZ      =  locations of voxels (mm)
% CVA.xyz      =  seed voxel location (mm)
% CVA.VOX      =  dimension of voxels (mm)
% 
% CVA.V        =  canonical vectors  (data)
% CVA.v        =  canonical variates (data)
% CVA.W        =  canonical vectors  (design)
% CVA.w        =  canonical variates (design)
% CVA.C        =  canonical contrast (design)
% 
% CVA.chi      =  Chi-squared statistics testing D >= i
% CVA.df       =  d.f.
% CVA.p        =  p-values
%
% also saved in CVA_*.mat in the SPM working directory
%
%__________________________________________________________________________
% 
% This routine allows one to make inferences about effects that are
% distributed in a multivariate fashion or pattern over voxels. It uses
% conventional canonical variates (CVA) analysis (also know as canonical
% correlation analysis, ManCova and linear discriminant analysis).  CVA is
% a complement to MVB, in that the predictor variables remain the design
% matrix and the response variable is the imaging data in the usual way.
% However, the multivariate aspect of this model allows one to test for
% designed effects that are distributed over voxels and thereby increase
% the sensitivity of the analysis.
% 
% Because there is only one test, there is no multiple comparison problem.
% The results are shown in term of the maximum intensity projection of the
% (positive) canonical image or vector and the canonical variates based on
% (maximally) correlated mixtures of the explanatory variables and data.
% 
% CVA uses the generalised eigenvalue solution to the treatment and
% residual sum of squares and products of a general linear model. The
% eigenvalues (i.e., canonical values), after transformation, have a
% chi-squared distribution and allow one to test the null hypothesis that
% the mapping is D or more dimensional. This inference is shown as a bar
% plot of p-values.  The first p-value is formally identical to that
% obtained using Wilks' Lambda and tests for the significance of any
% mapping.
% 
% This routine uses the current contrast to define the subspace of interest
% and treats the remaining design as uninteresting. Conventional results
% for the canonical values are used after the data (and design matrix) have
% been whitened; using the appropriate ReML estimate of non-sphericity.
% 
% CVA can be used for decoding because the model employed by CVA does not
% care about the direction of the mapping (hence canonical correlation
% analysis). However, one cannot test for mappings between nonlinear
% mixtures of regional activity and some experimental variable (this is
% what the MVB was introduced for).
% 
% References:
% 
% Characterizing dynamic brain responses with fMRI: a multivariate
% approach. Friston KJ, Frith CD, Frackowiak RS, Turner R. NeuroImage. 1995
% Jun;2(2):166-72.
%
% A multivariate analysis of evoked responses in EEG and MEG data. Friston
% KJ, Stephan KM, Heather JD, Frith CD, Ioannides AA, Liu LC, Rugg MD,
% Vieth J, Keber H, Hunter K, Frackowiak RS. NeuroImage. 1996 Jun;
% 3(3):167-174.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_cva.m 4352 2011-06-13 17:27:46Z ged $
 

% get figure handles
%--------------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
spm_results_ui('Clear');
spm_input('!DeleteInputObj');
header = get(Finter,'Name');
set(Finter,'Name','Canonical Variates analysis')

% review old analysis or proceed with a new one
%--------------------------------------------------------------------------
if nargin < 4
    action = spm_input('Canonical Variates Analysis','!+1','b', ...
        {'New Analysis','Results'},{'A','R'},1);
    if strcmpi(action,'r')
        CVA = spm_cva_results;
        return
    end
end
 
%-Contrast and VOI specification
%==========================================================================
 
% get contrast
%--------------------------------------------------------------------------
con    = SPM.xCon(xSPM.Ic).name;
c      = SPM.xCon(xSPM.Ic).c;
c      = full(c);
 
% get VOI name
%--------------------------------------------------------------------------
try
    name = CVA.name;
catch
    name = spm_input('name','-8','s',con);
end
name   = strrep(name,' ','_');
name   = ['CVA_' name '.mat'];
 
% get current location {mm}
%--------------------------------------------------------------------------
try
    xyzmm = CVA.xY.xyz;
catch
    xyzmm = spm_results_ui('GetCoords');
end
    
% specify search volume
%--------------------------------------------------------------------------
try
    xY = CVA.xY;
    CVA = rmfield(CVA,'xY');
catch
    xY = [];
end
xY.xyz = xyzmm;

Q      = ones(1,size(SPM.xVol.XYZ, 2));
XYZmm  = SPM.xVol.M(1:3,:)*[SPM.xVol.XYZ; Q];

[xY, XYZ, j] = spm_ROI(xY, XYZmm);
 
% voxels defined
%--------------------------------------------------------------------------
spm('Pointer','Watch')
 
%-Extract required data from results files
%==========================================================================
 
% get explanatory variables (data)
%--------------------------------------------------------------------------
Y    = spm_get_data(SPM.xY.VY,SPM.xVol.XYZ(:,j));
 
if isempty(Y)
    spm('alert*',{'No voxels in this VOI';'Please use a larger volume'},...
        'Canonical Variates analysis');
    return
end
 
% remove serial correlations and get design (note X := W*X)
%--------------------------------------------------------------------------
Y   = SPM.xX.W*Y;
X   = SPM.xX.xKXs.X;
 
% and get null-space of contrast
%--------------------------------------------------------------------------
X0  = X - X*c*pinv(c);
try,  X0 = [X0 SPM.xX.K.X0]; end          % add drift terms
try,  X0 = [X0 SPM.xGX.gSF]; end          % add global estimate
X   = full(X*c);
X0  = spm_svd(X0);
 

%-Dimension reduction (if necessary)
%==========================================================================
U   = spm_mvb_U(Y,'compact',X0,XYZ);
U   = spm_svd(U);
Y   = Y*U;

 
%-Canonical Variates Analysis
%==========================================================================

% remove null space of contrast
%--------------------------------------------------------------------------
Y     = Y - X0*(X0'*Y);
X     = X - X0*(X0'*X);
P     = pinv(X);
 
% degrees of freedom
%--------------------------------------------------------------------------
[n,m] = size(Y);
b     = rank(X);
h     = min(b,m);
f     = n - b - size(X0,2);
 
 
% generalised eigensolution for treatment and residual sum of squares
%--------------------------------------------------------------------------
T     = X*(P*Y);
SST   = T'*T;
SSR   = Y - T;
SSR   = SSR'*SSR;
[v,d] = eig(SSR\SST);
[q,r] = sort(-real(diag(d)));
r     = r(1:h);
d     = real(d(r,r));
v     = real(v(:,r));
V     = U*v;                          % canonical vectors  (data)
v     = Y*v;                          % canonical variates (data)
W     = P*v;                          % canonical vectors  (design)
w     = X*W;                          % canonical variates (design)
C     = c*W;                          % canonical contrast (design)
 
% inference on dimensionality - p(i) test of D >= i; Wilks' Lambda := p(1)
%--------------------------------------------------------------------------
cval  = log(diag(d) + 1);
for i = 1:h
  chi(i) = (f - (m - b + 1)/2)*sum(cval(i:h));
  df(i)  = (m - i + 1)*(b - i + 1);
  p(i)   = 1 - spm_Xcdf(chi(i),df(i));
end
 
% prevent overflow
%--------------------------------------------------------------------------
p     = max(p,exp(-16));

%-Save results
%==========================================================================
M     = SPM.xVol.M(1:3,1:3); %-voxels to mm matrix
VOX   = sqrt(diag(M'*M))';   %-voxel dimensions 

% assemble results
%--------------------------------------------------------------------------
CVA.contrast = con;          % contrast name
CVA.name     = name;         % CVA name
CVA.c        = c;            % contrast weights
CVA.X        = X;            % contrast subspace
CVA.Y        = Y;            % whitened and adjusted data
CVA.X0       = X0;           % null space of contrast
 
CVA.XYZ      = XYZ;          % locations of voxels (mm)
CVA.xyz      = xyzmm;        % seed voxel location (mm)
CVA.VOX      = VOX;          % dimension of voxels (mm)
 
CVA.V        = V;            % canonical vectors  (data)
CVA.v        = v;            % canonical variates (data)
CVA.W        = W;            % canonical vectors  (design)
CVA.w        = w;            % canonical variates (design)
CVA.C        = C;            % canonical contrast (design)
 
CVA.chi      = chi;          % Chi-squared statistics testing D >= i
CVA.df       = df;           % d.f.
CVA.p        = p;            % p-values
 
 
% save
%--------------------------------------------------------------------------
if spm_check_version('matlab','7') >= 0
    save(fullfile(SPM.swd,name),'-V6','CVA')
else
    save(fullfile(SPM.swd,name),'CVA')
end
assignin('base','CVA',CVA)
 
% display results
%--------------------------------------------------------------------------
spm_cva_results(CVA);

% reset title
%--------------------------------------------------------------------------
set(Finter,'Name',header)
spm('Pointer','Arrow')
