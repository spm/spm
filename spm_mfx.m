function [SPM] = spm_mfx(SPM)
% converts a 1st-level design specification into a MFX specification
% FORMAT [SPM] = spm_mfx(SPM)
%___________________________________________________________________________
% %W% Karl Friston %E%
SCCSid = '%I%';

%-Say hello
%-----------------------------------------------------------------------
SPMid  = spm('FnBanner',mfilename,SCCSid);
Finter = spm('FigName','MFX specification...'); spm('Pointer','Arrow')

%-Get SPM.mat if necessary
%-----------------------------------------------------------------------
if nargin ==0
	load(spm_get(1,'SPM.mat','Select SPM.mat'));
end
swd   = SPM.swd;


%-Check this is a repeated measures design
%-----------------------------------------------------------------------
n     = length(SPM.Sess);			% number of sessions
for i = 1:n
	if length(SPM.Sess(i).col) ~= length(SPM.Sess(1).col)
		error('this is not a repeated measures design')
		return
	end
end

% Build MFX design specification
%-----------------------------------------------------------------------

% response variable (parameter estimates from 1st-level)
%=======================================================================
nY       = length(SPM.xX.iC);
nP       = nY/n;
S.xY.VY  = SPM.Vbeta(SPM.xX.iC);
for i = 1:nY
	S.xY.VY(i).fname = fullfile(swd,S.xY.VY(i).fname);
end

% designs
%=======================================================================

% 1st-level (X1) and confounds X0 (inlcuding those in filter structure);
%-----------------------------------------------------------------------
iX1   = SPM.xX.iC;
iX0   = SPM.xX.iB;
X1    = SPM.xX.X(:,iX1);
X0    = SPM.xX.X(:,iX0);
K0    = sparse(0,0);
try
	for i = 1:n
		K0 = blkdiag(K0, SPM.xX.K(i).X0);
	end
end
X0    = [X0 K0];
X1    = X1 - X0*inv(X0'*X0)*(X0'*X1);
X1    = sparse(X1);


% 2nd-level (X2) design (as an F-contrast to ensure estimability)
%-----------------------------------------------------------------------
sX       = spm_sp('set',eye(n));
c        = ones(n,1);
xCon     = spm_FcUtil('Set','one-sample t-test','F','c',c, sX);
xX.xKXs  = sX;
for    i = 1:n
	xX.name{i}  = sprintf('Session %i',i);
end
[I,xCon] = spm_conman(xX,xCon,'F',1,'2nd-level contrast','',1)
X2       = xCon(I).c;
X2       = kron(X2,speye(nP,nP));
nC       = size(X2,2);

% names
%-----------------------------------------------------------------------
name  = {};
for i = 1:nC
for j = 1:nP
	name{end + 1} = sprintf('contrast %i parameter %i',i,j);
end
end

% set fields
%-----------------------------------------------------------------------
S.xX.X    = X2;
S.xX.name = name;


% mixed covariance components
%=======================================================================

% 1st-level covariance components 
%-----------------------------------------------------------------------
Q     = SPM.xVi.Vi;

% 2nd-level covariance components (projected to first level)
%-----------------------------------------------------------------------
for i = 1:nP

	% unequal variances
	%---------------------------------------------------------------
	s          = zeros(nP,nP);
	s(i,i)     = 1;
	Q{end + 1} = X1*kron(speye(n,n),s)*X1';

	% correlations
	%---------------------------------------------------------------
	for  j = (i + 1):nP
		s          = zeros(nP,nP);
		s(i,j)     = 1;
		s(j,i)     = 1;
		Q{end + 1} = X1*kron(speye(n,n),s)*X1';
	end
end

% 1st-level non-sphericity
%-----------------------------------------------------------------------
V1       = spm_reml(SPM.xVi.Cy,[X1*X2 X0],Q);

% 2nd-level non-sphericity
%-----------------------------------------------------------------------
W        = SPM.xX.W;
K        = SPM.xX.K;
pX1      = SPM.xX.pKX;
V2       = pX1*spm_filter(K,spm_filter(K,W*V1*W')')*pX1';
S.xVi.V  = sparse(V2(iX1,iX1));


% smoothness and volume infomation
%-----------------------------------------------------------------------
S.xVol   = SPM.xVol;

% implicit masking
%-----------------------------------------------------------------------
S.xM.TH  = -ones(nY,1)/0
S.xM.I   = 1;


%-Change to SPM.swd/mfx and save analysis parameters in SPM.mat file
%-----------------------------------------------------------------------
mkdir mfx;
cd mfx
SPM.swd = pwd;
save SPM SPM

%=======================================================================
%- E N D: Cleanup GUI
%=======================================================================
fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')               %-#
spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                     %-#
fprintf('...you may now estimate this mixed-effects model\n\n')      %-#
