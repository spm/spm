function [CVA] = spm_cva(xSPM,SPM,hReg)
% VOI extraction of adjusted data and CVA
% FORMAT [CVA] = spm_cva(xSPM,SPM,hReg);
%
% xSPM   - structure containing specific SPM details
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
% obtained using Wilk’s Lambda and tests for the significance of any
% mapping.
% 
% This routine uses the current contrast to define the subspace of interest
% and treats the remaining design as uninteresting. Conventional results
% for the canonical values are used after the data (and design matrix) have
% been whitened; using the appropriate ReML estimate of non-sphericity.
% 
% CVA can be used to for decoding because the model employed by CVA design
% not care about the direction of the mapping (hence canonical correlation
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
% Vieth J, Keber H, Hunter K, Frackowiak RS. NeuroImage. 1996 Jun;3(3 Pt
% 1):167-74.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
 
% Karl Friston
% $Id: spm_cva.m 2662 2009-01-28 20:23:11Z karl $
 

% review old analysis or proceed with a new one
%--------------------------------------------------------------------------
switch questdlg('new canonical variaiates analaysis?');
    case {'No'}
        CVA = spm_cva_results;
        return
    case {'Cancel'}
        return
end

% get figure handles
%--------------------------------------------------------------------------
Fcva   = spm_figure('GetWin','MVB');
Finter = spm_figure('GetWin','Interactive');
header = get(Finter,'Name');
set(Finter,'Name','Canonical Variates analysis')
 
 
% contrast and VOI specification
%==========================================================================
 
%-Get contrast
%--------------------------------------------------------------------------
con    = SPM.xCon(xSPM.Ic).name;
c      = SPM.xCon(xSPM.Ic).c;
c      = full(c);
 
%-Get VOI name
%--------------------------------------------------------------------------
name   = ['CVA_' spm_input('name','-8','s',con)];
 
%-Get current location {mm}
%--------------------------------------------------------------------------
xyzmm  = spm_results_ui('GetCoords');
 
%-Specify search volume
%--------------------------------------------------------------------------
SPACE  = spm_input('Search volume...','!+1','b',...
                  {'Sphere','Box','Image'},['S','B','I']);
Q      = ones(1,size(SPM.xVol.XYZ, 2));
XYZmm  = SPM.xVol.M*[SPM.xVol.XYZ; Q];
XYZmm  = XYZmm(1:3,:);
 
switch SPACE
 
    case 'S' %-Sphere
    %----------------------------------------------------------------------
    D     = spm_input('radius of VOI {mm}','!+1');
    str   = sprintf('%0.1fmm sphere',D);
    j     = find(sum((XYZmm - xyzmm*Q).^2) <= D^2);
 
    case 'B' %-Box
    %----------------------------------------------------------------------
    D     = spm_input('box dimensions [k l m] {mm}','!+1');
    if length(D) < 3
        D = D(1)*[1 1 1];
    end
    str   = sprintf('%0.1f x %0.1f x %0.1f mm box',D(1),D(2),D(3));
    j     = find(all(abs(XYZmm - xyzmm*Q) <= D(:)*Q/2));
 
    case 'I' %-Mask Image
    %----------------------------------------------------------------------
    Msk   = spm_select(1,'image','Image defining search volume');
    D     = spm_vol(Msk);
    str   = sprintf('image mask: %s',spm_str_manip(Msk,'a30'));
    XYZ   = D.mat \ [XYZmm; Q];
    j     = find(spm_sample_vol(D, XYZ(1,:), XYZ(2,:), XYZ(3,:),0) > 0);
 
end
 
 
% voxels defined
%--------------------------------------------------------------------------
spm('Pointer','Watch')
 
%-Extract required data from results files
%==========================================================================
 
% get explanatory variables (data)
%--------------------------------------------------------------------------
XYZ  = XYZmm(:,j);
Y    = spm_get_data(SPM.xY.VY,SPM.xVol.XYZ(:,j));
 
if isempty(Y)
    warndlg({'No voxels in this VOI';'Please use a larger volume'})
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
[n,b] = size(X);                       
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
 
% inference on dimensionality - p(i) test of D >= i; Wilk’s Lambda := p(1)
%--------------------------------------------------------------------------
cval  = log(diag(d) + 1);
for i = 1:h
  chi(i) = (f - (m - b + 1)/2)*sum(cval(i:h));
  df(i)  = (m - i + 1)*(b - i + 1);
  p(i)   = 1 - spm_Xcdf(chi(i),df(i));
end
 

% save results
%==========================================================================
 
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
CVA.VOX      = xSPM.VOX;     % dimension of voxels (mm)
 
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
save(fullfile(SPM.swd,name),'CVA')
assignin('base','CVA',CVA)
 
% display results
%--------------------------------------------------------------------------
spm_cva_results(CVA);

%-Reset title
%--------------------------------------------------------------------------
set(Finter,'Name',header)
spm('Pointer','Arrow')
