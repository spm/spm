function [SPM] = spm_mfx(SPM,c)
% Convert a 1st-level design specification into a MFX specification
% FORMAT [SPM] = spm_mfx(SPM,c)
% SPM {in} - design and estimation structure after a 1st-level analysis
% c        - contrast used to define 2nd level design matrix. If this is
%            not specified spm_mfx will (1) suggest the ones(n,1) contrast
%            where n is the number of sessions/subjects, (2) call
%            spm_conman to allow this contrast to be modified interactively
%
%            Note: the specification of a contrast that is not ones(n,1) allows, 
%            for example, specified sessions/subjects to be ignored.
%            
% SPM {out} is saved in fullfile(SPM.swd,'mfx','SPM.mat')
%
% spm_mfx takes the SPM.mat of a 1st-level estimation of a repeated-measure
% multi-session study and produces the SPM design specification for a
% full mixed-effects (MFX) analysis. The 1st-level design (X1) must have
% the same number of parameters for each session. These are assumed to
% represent session-specific realisations of 2nd-level effects.
%
% spm_mfx prompts for a 2nd-level design matrix (X2) in the form of an
% F-contrast. This is expanded using the Kronecker tensor product to
% model the effects of each 2nd-level parameter separately. A new
% SPM.mat structure is saved in a subdirectory of the 1st-level results
% directory and can be estimated in the usual way. 2nd-level contrasts
% can then be used to test specific hypotheses at the 2nd-level in terms
% of compounds of 1st-level parameters specified by X2 (e.g. their
% mean).
%
% spm_mfx is a full mixed effects analysis in the sense that it allows
% for unbalanced designs at the 1st-level and different 1st-level error
% covariances. Operationally, ReML estimates of the 1st and 2nd-level
% covariance components are computed by projecting the 2nd-level effects
% down to the 1st-level and partitioning the covariance of the data in
% observation space. The 2nd-level parameter estimates are then computed
% as linear mixtures of the 1st-level estimates, using the appropriate
% non-sphericity. This non-sphericity is a mixture of 1st- and 2nd-level
% components that renders the ensuing 2nd-level estimates ML.
%
% In summary;
%
% ReML estimates of V1 are obtained where
%
%                y    = X1*B1 + X0*B0 + e1
%                B1   = X2*B2 + e2;
%
% giving;        y    = X1*X2*B2 + X0*B0 + X1*e2 + e1
%
% where          V1   = cov(X1*e2 + e1)
%
% V1 is now used to give the covariance components of any 1st-level
% parameter estimators B1h
%
%                B1h  = M1*y
% such that      V2   = cov(B1h) = M1*V1*M1'
%
% is the error covariance for the single level model
%
%                B1h  = X2*B2 + r2
%
% where cov(r2) = cov(B1h) = V2, which can be estimated non-iteratively
% in the usual way to give the ML estimates of B2.
%
% Note that with balanced designs and equal error covariances over
% sessions, at the 1st level there is no need to compute multiple
% covariance components because, at the 2nd-level, they are exactly the
% same (i.e. M1*X1*cov(e2)*X1*M1 has the same form as M1*cov(e1)*M1).
%
% The ReML hyperparameters are estimated using the covariance of y over
% voxels. This means that the relative amounts of within and
% between-session variance are assumed to be fixed over voxels but can
% vary in their overall expression. The voxels used for this pooling are
% those that show 1st-level responses.
%
% See spm_reml.m
%
%__________________________________________________________________________
% Copyright (C) 2002-2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_mfx.m 6380 2015-03-17 16:13:50Z guillaume $

SVNid = '$Rev: 6380 $';

%-Say hello
%--------------------------------------------------------------------------
SPMid  = spm('FnBanner',mfilename,SVNid);
Finter = spm('FigName','MFX specification...'); spm('Pointer','Arrow')

%-Get SPM.mat if necessary
%--------------------------------------------------------------------------
if ~nargin
    load(spm_select(1,'^SPM\.mat$','Select SPM.mat'));
end
try
    swd = SPM.swd;
catch
    error('This model has not been estimated.');
end

%-Change to SPM.swd
%--------------------------------------------------------------------------
try
    cd(swd);
end
 
%-Check this is a repeated measures design
%--------------------------------------------------------------------------
if isfield(SPM,'Sess')
    n     = length(SPM.Sess);            % number of sessions
    for i = 1:n
        if length(SPM.Sess(i).col) ~= length(SPM.Sess(1).col)
            error('This is not a repeated measures design.')
        end
    end
else
    % This is not first level fMRI so ask for number of subjects
    n     = spm_input('Number of subjects', '+1', 'r', [], 1);
end


%-Build MFX design specification
%==========================================================================

%-Response variable (parameter estimates of B1 from the 1st-level)
%--------------------------------------------------------------------------
iX1      = [SPM.xX.iH SPM.xX.iC];
iX0      = [SPM.xX.iB SPM.xX.iG];
nY       = length(iX1);              % number of B1
nP       = nY/n;                     % number of B2

%-Change relative filenames to full
%--------------------------------------------------------------------------
S.xY.P   = spm_file(char(SPM.Vbeta(iX1).fname),'path',swd);
S.xY.VY  = spm_vol(S.xY.P);

%-Design matrices
%==========================================================================

% 1st-level (X1) and confounds X0 (including those in filter structure)
%--------------------------------------------------------------------------
X1    = SPM.xX.X(:,iX1);
X0    = SPM.xX.X(:,iX0);
K0    = sparse(0,0);
try
    for i = 1:n
        K0 = blkdiag(K0, SPM.xX.K(i).X0);
    end
end
X0    = [X0 K0];

%-Ensure X1 is orthogonalized w.r.t. X0
%--------------------------------------------------------------------------
X1    = X1 - X0*inv(X0'*X0)*(X0'*X1);
X1    = sparse(X1);


% 2nd-level (X2) design (as an F-contrast to ensure estimability)
% expanded to cover all parameters (kron(X2,eye(nP)))
%--------------------------------------------------------------------------
sX        = spm_sp('set',eye(n));
if nargin < 2
    c         = ones(n,1);
    SPM2.xCon = spm_FcUtil('Set','one-sample t-test','F','c',c, sX);
    SPM2.xX.xKXs  = sX;
    for    i  = 1:n
        SPM2.xX.name{i} = sprintf('Session %i',i);
    end
    tmpdir = tempname;
    try
        mkdir(tmpdir);
        cd(tmpdir);
        [I,xCon]  = spm_conman(SPM2,'F',1,'2nd-level contrast','',1);
    end
    try, cd(swd); end
    try, spm_unlink(fullfile(tmpdir,'SPM.mat')); end
    try, rmdir(tmpdir); end
else
    I = 1;
    xCon(I).c = c;
    xCon(I).name = '2nd level';
end
X2        = xCon(I).c;
X2        = kron(X2,speye(nP,nP));
nC        = size(X2,2);


%-Construct 2nd-level SPM specification (S)
%==========================================================================
fprintf('%-40s: %30s\n','Mixed-Effect Model','...ReML estimation');     %-#
spm('FigName','Stats: MFX-ReML',Finter); spm('Pointer','Watch')

xsDes.Design = '2nd-level MFX analysis';
xsDes.Name   = xCon(I).name;
S.xsDes      = xsDes;       % description


%-Names for nC contrasts in X2 and nP parameters
%--------------------------------------------------------------------------
name  = {};
for i = 1:nC
for j = 1:nP
    name{end + 1} = sprintf('contrast %i parameter %i',i,j);
end
end

sF    = {'parameter','session','',''};
I     = [];
for i = 1:n
for j = 1:nP
    I(end + 1,:)  = [j i 1 1];
end
end

%-Set fields
%--------------------------------------------------------------------------
S.xX.X    = X2;
S.xX.name = name;
S.xX.iH   = [];
S.xX.iC   = [1:size(X2,2)];
S.xX.iB   = [];
S.xX.iG   = [];
S.xX.I    = I;
S.xX.sF   = sF;


%-Mixed covariance components
%==========================================================================

% 1st-level covariance components 
%--------------------------------------------------------------------------
try
    Q = SPM.xVi.Vi;
catch
    Q = {SPM.xVi.V};
end


% 2nd-level covariance components (projected to first level)
%--------------------------------------------------------------------------
for i = 1:nP

    % unequal variances
    %----------------------------------------------------------------------
    s          = zeros(nP,nP);
    s(i,i)     = 1;
    Q{end + 1} = X1*kron(speye(n,n),s)*X1';

    % correlations
    %----------------------------------------------------------------------
    for  j = (i + 1):nP
        s          = zeros(nP,nP);
        s(i,j)     = 1;
        s(j,i)     = 1;
        Q{end + 1} = X1*kron(speye(n,n),s)*X1';
    end
end

% 1st-level non-sphericity - ReML estimates, restricted to the Null
% space of 'fixed' effects X1*X2 and X0
%--------------------------------------------------------------------------
[V1,h]    = spm_reml(SPM.xVi.Cy,[X1*X2 X0],Q);

% 2nd-level non-sphericity (including original whitening and filtering)
%--------------------------------------------------------------------------
W         = SPM.xX.W;
try
    K     = SPM.xX.K;
catch
    K     = 1;
end
pX1       = SPM.xX.pKX;
M1        = pX1(iX1,:)*spm_filter(K,W);

%-Project 1st-level covariance components to 2nd-level and save in xVi
%--------------------------------------------------------------------------
V2        = M1*V1*M1';
V2        = V2*length(V2)/trace(V2);
for i = 1:length(Q);
    Vi{i} = M1*Q{i}*M1';
end

S.xVi.V   = sparse(V2);
S.xVi.Vi  = Vi;
S.xVi.h   = h;

%-Smoothness and volume information
%--------------------------------------------------------------------------
S.xVol    = SPM.xVol;

%-Change to SPM.swd/mfx and save analysis parameters in SPM.mat file
%--------------------------------------------------------------------------
SPM       = S;
SPM.swd   = fullfile(swd,'mfx');
[st, me]  = mkdir(SPM.swd);
if st
    save(fullfile(SPM.swd,'SPM.mat'), 'SPM', spm_get_defaults('mat.format'));
else
    error('Could not save SPM.mat in mfx: %s', me)
end

%==========================================================================
%- E N D: Cleanup GUI
%==========================================================================
spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#
fprintf('...you may now estimate this mixed-effects model\n\n')         %-#
