function spm_spm_Bayes
% Conditional parameter estimation of a General Linear Model
% FORMAT spm_spm_Bayes
%_______________________________________________________________________
% spm_spm_Bayes returns to voxels identified by spm_spm (OLS parameter
% estimation) to get conditional parameter estimates and ReML hyper-
% parameter estimates.  These estimates use prior covariances, on the
% parameters, from emprical Bayes.  These PEB prior variances come from
% the hierarchical model that obtains by considering voxels as providing
% a second level.  Put simply, the variance in parameters, over voxels,
% is used as a prior variance from the point of view of any one voxel.
% The error covariance hyperparameters are re-estimated in the light of
% these priors.  The approach adopted is essentially a fully Bayesian
% analysis at each voxel, using emprical Bayesian prior variance
% estimators over voxels.
%
% Each separable partition (i.e. Session) is assigned its own
% hyperparameter but within session covariance components are lumped
% together, using their relative expectations over voxels.  This makes
% things much more computationally efficient and avoids inefficient
% voxel-specific multiple hyperparameter estimates.
%
% spm_spm_Bayes saves:
%
%                           ----------------
%
%
%	PPM.l           = session-specific hyper-parameter means
%	PPM.C(v,v)      = conditional covariances of parameters;
%	PPM.dC{i}(v,v)  = dC/dl;
%	PPM.ddC{i}(v,v) = ddC/dldl
%
% The derivatives are used to compute the conditional variance of various
% contrasts in spm_getSPM, using a second order Taylor expansion about the
% hyperparameter means.
%
%
%                           ----------------
%
% Cbeta_????.{img,hdr}                     - conditional parameter images
% These are 16-bit (float) images of the conditional estimates. The image
% files are numbered according to the corresponding column of the
% design matrix. Voxels outside the analysis mask (mask.img) are given
% value NaN.
%
%                           ----------------
%
% CHp_????.{img,hdr}              - error covariance hyperparamter images
% This is a 32-bit (double) image of the ReML error variance estimate.
% for each separable partition (Session).  Voxels outside the analysis 
% mask are given value NaN.
%
%_______________________________________________________________________
% %W% Karl Friston %E%

%-Say hello
%-----------------------------------------------------------------------
Finter = spm('FigName','Stats: Bayesian estimation...');

%-Select SPM.mat & change directory
%-----------------------------------------------------------------------
swd    = spm_str_manip(spm_get(1,'SPM.mat','Select SPM.mat'),'H');
load(fullfile(swd,'SPM.mat'))
spm('Pointer','Watch')
cd(swd)


%-Parameters
%-----------------------------------------------------------------------
global MAXMEM 
if length(MAXMEM)	
	maxMem = MAXMEM;	
else
	maxMem = 2^20;	%-Max data chunking size, in bytes
end
M      = VY(1,1).mat;
DIM    = VY(1,1).dim(1:3)';
N      = 3 - sum(DIM == 1);


%=======================================================================
% - A N A L Y S I S   P R E L I M I N A R I E S
%=======================================================================

%-Initialise output images
%=======================================================================
fprintf('%-40s: %30s','Output images','...initialising')             %-#

%-Image dimensions
%-----------------------------------------------------------------------
xdim  = DIM(1); ydim = DIM(2); zdim = DIM(3);

%-Intialise oonditional estimate image files
%-----------------------------------------------------------------------
clear Vbeta
[nScan nBeta]  = size(xX.X);
Vbeta(1:nBeta) = deal(struct(...
			'fname',	[],...
			'dim',		[DIM',spm_type('float')],...
			'mat',		M,...
			'pinfo',	[1 0 0]',...
			'descrip',	''));
for i = 1:nBeta
	Vbeta(i).fname   = sprintf('Cbeta_%04d.img',i);
	Vbeta(i).descrip = sprintf('Cond. beta (%04d) - %s',i,xX.Xnames{i});
	spm_unlink(Vbeta(i).fname)
	Vbeta(i)         = spm_create_image(Vbeta(i));
end

%-Intialise ReML hyperparameter image files
%-----------------------------------------------------------------------
if exist('Sess')
	nHp = length(Sess);
else
	nHp = 1;
end
VHp(1:nHp)  = deal(struct(...
			'fname',	[],...
			'dim',		[DIM',spm_type('double')],...
			'mat',		M,...
			'pinfo',	[1 0 0]',...
			'descrip',	''));
for i = 1:nHp
	VHp(i).fname   = sprintf('Hp_%04d.img',i);
	VHp(i).descrip = sprintf('Hyperparameter (%04d)',i);
	spm_unlink(VHp(i).fname)
	VHp(i)         = spm_create_image(VHp(i));
end

fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...initialised')        %-#


%=======================================================================
% - E M P I R I C A L  B A Y E S  F O R  P R I O R  V A R I A N C E
%=======================================================================
fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...estimatng priors')   %-#

% get row u{i} and column {vi}/v0{i} indices for separable designs
%----------------------------------------------------------------------
if exist('Sess')
	s     = length(Sess);
	for i = 1:s
		 u{i} = Sess{i}.row;
		 v{i} = Sess{i}.col;
		v0{i} = xX.iB(i);

	end
else
            s = 1;
	 u{1} = [1:nScan];
	 v{1} = [xX.iH xX.iC];
	v0{1} = [xX.iB xX.iG]
end

% cycle over sperarable paritions
%-----------------------------------------------------------------------
for i = 1:s

	% Get design X and confounds X0
	%---------------------------------------------------------------
	fprintf('%-30s- %i\n','  ReML Session',i);
	X      = xX.X(u{i}, v{i});
	X0     = xX.X(u{i},v0{i});
	[m n]  = size(X);

	% add confound in 'filter'
	%---------------------------------------------------------------
	if iscell(xX.K)
		X0 = full([X0 xX.K{i}.KH]);
	end

	% orthogonalize X w.r.t. X0
	%---------------------------------------------------------------
	X      = X - X0*(pinv(X0)*X);

	% covariance components induced by parameter variation over voxels
	%---------------------------------------------------------------
	for  j = 1:n
		Q{j} = X*sparse(j,j,1,n,n)*X';
	end

	% covariance components induced by error non-sphericity {Q}
	%---------------------------------------------------------------
	V      = {};
	if iscell(xX.xVi.Vi)
		for j = 1:length(xX.xVi.Vi)
			q    = xX.xVi.Vi{j}(u{i},u{i});
			if any(find(q));
				V{end + 1} = q;
			end
		end
	else
		V     = {speye(m,m)};
	end


	% ReML ovariance component estimation
	%---------------------------------------------------------------
	[C h W] = spm_reml(CY,X0,{Q{:} V{:}});


	% lump error covariance components togther to form V
	%---------------------------------------------------------------
	if length(V) > 1
		q     = sparse(m,m);
		for j = 1:length(V)
			q = q + V{j}*h(n + j);
		end
		V     = {q*m/trace(q)};
	end

	%-design structure for this partition using prior variances sP(i)
	% treat confounds as fixed (i.e. infinite prior variance)
	%---------------------------------------------------------------
	n0      = size(X0,2);
	Cb      = blkdiag(diag(h(1:n)),speye(n0,n0)*1e8);
	P{1}.X  = [X X0];
	P{1}.C  = V;
	P{2}.X  = sparse(size(P{1}.X,2),1);
	P{2}.C  = Cb;

	sP(i).P = P;
	sP(i).u = u{:};
	sP(i).v = [v{:} v0{:}];
	sP(i).h = h;

end



%=======================================================================
% - F I T   M O D E L   &   W R I T E   P A R A M E T E R    I M A G E S
%=======================================================================

%-Cycle to avoid memory problems (plane by plane)
%=======================================================================
spm_progress_bar('Init',100,'Bayesian estimation','');

%-Find a suitable block size for loop over planes (2D or 3D data)
%-----------------------------------------------------------------------
blksz = ceil(maxMem/8/nScan);
SHp   = 0;
for z = 1:zdim

    % current plane-specific parameters
    %-------------------------------------------------------------------
    U       = find(XYZ(3,:) == z);
    nbch    = ceil(length(U)/blksz);
    CrBl    = zeros(nBeta,length(U));	%-conditional parameter estimates
    CrHp    = zeros(nHp,  length(U));	%-ReML hyperparameter estimates
    for bch = 1:nbch			%-loop over bunches of lines (planks)

	%-construct list of voxels in this block
	%---------------------------------------------------------------
	I     = [1:blksz] + (bch - 1)*blksz;
	I     = I(I <= length(U));
	xyz   = XYZ(:,U(I));
	nVox  = size(xyz,2);

	%-Get response variable
	%---------------------------------------------------------------
	Y     = zeros(nScan,nVox);
	for i = 1:nScan
		Y(i,:) = spm_sample_vol(VY(i,1),xyz(1,:),xyz(2,:),xyz(3,:),0);
	end

	%-Conditional estimates (per partition, per voxel)
	%---------------------------------------------------------------
	beta  = zeros(nBeta,nVox);
	Hp    = zeros(nHp,  nVox);
	for j = 1:length(sP)
		P     = sP(j).P;
		u     = sP(j).u;
		v     = sP(j).v;
		for i = 1:nVox
			[C P]      = spm_PEB(Y(u,i),P);
			beta(v,i)  = C{2}.E(1:length(v));
			Hp(j,i)    = C{1}.h;
		end
	end

	%-Save for current plane in memory as we go along
	%---------------------------------------------------------------
	CrBl(:,I) = beta;
	CrHp(:,I) = Hp;
	SHp       = SHp + sum(Hp,2);

    end   % (bch)


    %-write out plane data to image files
    %===================================================================

    %-Write conditional beta images
    %-------------------------------------------------------------------
    for i = 1:nBeta
	tmp       = sparse(XYZ(1,U),XYZ(2,U),CrBl(i,:),xdim,ydim);
   	tmp       = tmp + NaN*(~tmp);
	Vbeta(i)  = spm_write_plane(Vbeta(i),tmp,z);
    end

    %-Write hyperparameter images
    %-------------------------------------------------------------------
    for i = 1:nHp
	tmp       = sparse(XYZ(1,U),XYZ(2,U),CrHp(i,:),xdim,ydim);
   	tmp       = tmp + NaN*(~tmp);
	VHp(i)    = spm_write_plane(VHp(i),tmp,z);
    end


    %-Report progress
    %-------------------------------------------------------------------
    spm_progress_bar('Set',100*(z - 1)/zdim);


end % (for z = 1:zdim)
fprintf('\n')                                                        %-#
spm_progress_bar('Clear')

%=======================================================================
% - P O S T   E S T I M A T I O N
%=======================================================================

% Taylor expansion for conditional covariance
%-----------------------------------------------------------------------
fprintf('%-40s: %30s\n','Non-sphericity','...REML estimation') %-#

% expansion point (mean hyperparameters)
%-----------------------------------------------------------------------
l     = SHp/S;

% change on conditional coavriance w.r.t. hyperparameters
%-----------------------------------------------------------------------
n     = size(xX.X,2);
PPM.l = l;
for i = 1:s
	PPM.dC{i}  = sparse(n,n);
	PPM.ddC{i} = sparse(n,n);
end
for i = 1:s

	P     = sP(j).P;
	u     = sP(j).u;
	v     = sP(j).v;

	% derivatives
	%---------------------------------------------------------------
	d     = P{1}.X'*inv(P{1}.C{1})*P{1}.X;
	Cby   = inv(d/l(i) + inv(P{2}.C));
	d     = d*Cby;
	dC    = Cby*d/(l(i)^2);
	ddC   = 2*(dC/(l(i)^2) - Cby/(l(i)^3))*d;

	% place in output structure
	%---------------------------------------------------------------
	j               = 1:length(v);
	PPM.Cb(v,v)     = P{2}.C(j,j);
	PPM.Cby(v,v)    = Cby(j,j);
	PPM.dC{i}(v,v)  = dC(j,j);
	PPM.ddC{i}(v,v) = ddC(j,j);
end


%-"close" written image files, updating scalefactor information
%=======================================================================
fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...closing image files')  %-#
for i=1:nBeta, Vbeta(i) = spm_create_image(Vbeta(i)); end
for i=1:nHp  ,   VHp(i) = spm_create_image(VHp(i));   end

%-Retain relative filenames
%-----------------------------------------------------------------------
fprintf('%s%30s',sprintf('\b')*ones(1,30),'...tidying file handles') %-#
Vbeta  = {Vbeta.fname}';
VHp    = {  VHp.fname}';

fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')               %-#


%-Save remaining results files and analysis parameters
%=======================================================================
fprintf('%-40s: %30s','Saving results','...writing')                 %-#

%-Save analysis parameters in SPM.mat file
%-----------------------------------------------------------------------
PPMvars = {	'Vbeta','VHp',...		%-Filenames
		'PPM'};				%-Conditional covariances
		

save('PPM',PPMvars{:})

fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')               %-#


%=======================================================================
%- E N D: Cleanup GUI
%=======================================================================
spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                     %-#
fprintf('...use the results section for assessment\n\n')             %-#
