function [SPM,VOL,xX,xSDM] = spm_getSPM
% computs a specified and thresholded following parameter estimation
% FORMAT [SPM,VOL,xX,xSDM] = spm_getSPM;
%
% SPM    - structure containing SPM, distribution & filtering detals
% .swd   - SPM working directory - directory containing current SPM.mat
% .c     - contrast(s) - in cell array
% .Z     - minimum of n Statistics {filtered on u and k}
% .n     - number of conjoint tests        
% .STAT  - distribution {Z, T, X or F}     
% .df    - degrees of freedom [df{interest}, df{residual}]
% .u     - height threshold
% .k     - extent threshold {resels}
% .XYZ   - location of voxels {voxel coords}
% .XYZmm - location of voxels {mm}
% .QQ    - indices of volxes in Y.mad file
%
%
% VOL    - structure containing details of volume analysed
% .S     - search Volume {voxels}
% .R     - search Volume {resels}
% .FWHM  - smoothness {voxels}     
% .M     - voxels - > mm matrix
% .iM    - mm -> voxels matrix
% .VOX   - voxel dimensions {mm}
% .DIM   - image dimensions {voxels}
% 
%
% xX     - Design Matrix structure
%        - (see spm_spm.m for structure)
%
% xSDM   - structure containing contents of SPM.mat file
%          ( see spm_spm.m for contents...
%          ( ...except that xX & XYZ are removed
%
%_______________________________________________________________________
%
% spm_getSPM prompts for an SPM and applies thresholds {u & k}
% to a point list of voxel values (specified with their locations {XYZ})
% This allows the SPM be displayed and characterized in terms of regionally 
% significant effects by subsequent routines.
% 
% Either single contrasts can be examined or conjunctions of different
% contrasts.  In the latter case a new SPM is created, that reflects
% the miminum of all specified effects.  A conjunction is therefore
% the conjoint expression of two or more effects.  The contrasts are
% successively enforced to be orthogonal so the order of non-orthogonal
% contrast specification is important.
%
% Masking simply eliminates voxels from the current contrast if they
% do not survive an uncorrected p value (based on height) in one or
% more further contrasts.  No account is taken of this masking in the
% statistical inference pertaining to the masked contrast.
% 
% The SPM is subject to thresholding on the basis of height (u) and the
% number of resels comprising its clusters {k}. The height threshold is
% specified as above in terms of an [un]corrected p value or
% statistic.  Clusters can also be thresholded on the basis of their
% spatial extent.  (SPM99 expresses all cluster volumes in terms of
% resels.  These are effective resolution elements and comprise several
% voxels, depending on the smoothness estimator.) If you want to see
% all voxels simply enter 0.  In this instance the 'set-level'
% inference can be considered an 'omnibus test' based on the number of
% clusters that obtain.
% 
% see spm_list.m and spm_P.m for further details
%
%_______________________________________________________________________
% %W% Karl Friston, Andrew Holmes %E%


%-GUI setup
%-----------------------------------------------------------------------
spm_help('!ContextHelp',mfilename)

%-Select SPM.mat
%-----------------------------------------------------------------------
swd  = spm_str_manip(spm_get(1,'SPM.mat','Select SPM.mat for analysis'),'H');

%-What sort of SPM
%-----------------------------------------------------------------------
STAT = spm_input('...which SPM?','+1','b','SPM{T}|SPM{F}',['T','F'],1);


%-Preliminaries...
%=======================================================================

%-Get Stats data from SPM.mat
%-----------------------------------------------------------------------
xSDM = load(fullfile(swd,'SPM.mat'));
xX   = xSDM.xX;
XYZ  = xSDM.XYZ;
S    = xSDM.S;
R    = xSDM.R;

%-Build index from XYZ into corresponding Y.mad locations
%-----------------------------------------------------------------------
str = fullfile(swd,'Yidx.mat');
if exist(str), load(str), else, Yidx = []; end
QQ         = zeros(1,S);
QQ(Yidx)   = 1:length(Yidx);

%-Canonicalise relative pathnames in xSDM by prepending swd
% (SPM result images kept with relative pathnames so SPM is robust to
% (user moving files around, but use full pathnames here so code is robust
% (to changing directory whilst using the results section
%-----------------------------------------------------------------------
nbeta = length(xSDM.Vbeta);
for i=1:nbeta, xSDM.Vbeta(i).fname = fullfile(swd,xSDM.Vbeta(i).fname); end
xSDM.VResMS.fname = fullfile(swd,xSDM.VResMS.fname);


%=======================================================================
% - C O N T R A S T S ,  S P M   C O M P U T A T I O N ,   M A S K I N G
%=======================================================================
%-**** Could ask for contrast(s) here, and if matrix ask if F/T,
%-**** interpreting T contrast matrix as multiple contrasts for conjunctions

if STAT == 'T'

	%-Get contrast[s] as column vectors
	%---------------------------------------------------------------
	c = spm_input('contrast(s)','+1','x','',Inf,xX.xKXs)';
	n = size(c,2);

	%-Get any contrast(s) for mask
	%---------------------------------------------------------------
	if spm_input('mask with other contrast(s)','+1','y/n',[1,0],2)
		%-Get contrast[s] for mask
		%-------------------------------------------------------
		cm = spm_input('contrast(s)','+1','x','',Inf,xX.xKXs)';
		m  = size(cm,2);
		%-Threshold for mask (uncorrected)
		%-------------------------------------------------------
		u = spm_input('threshold for mask','+1','e',0.05);
		if u <= 1; u = spm_u(u,[1,xX.erdf],'T'); end
	else
		m  = 0;
		cm = [];
	end

	%-Enforce orthogonality of multiple contrasts for conjunction
	% (Orthogonality within subspace spanned by contrasts)
	%---------------------------------------------------------------
	for i = 2:n
		j      = find(c(:,i));
		tmp    = xX.Bcov(j,j)*c(j,1:i-1);
		c(j,i) = c(j,i) - tmp*(pinv(tmp)*c(j,i));
	end

	%-NB:	If contrasts are linearly dependant, will have a zero
	%	column in the contrast matrix.

	%-Compute SPM{t}(s) from parameter images and ResMS image
	%---------------------------------------------------------------
	Z   = zeros(n+m,S);
	tmp = [c,cm];

	%-Accumulate weighted sums of parameter estimates
	for j = 1:nbeta
		Z = Z + tmp(j,:)' * spm_sample_vol(xSDM.Vbeta(j),...
				XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
	end

	%-Normalise contrast(s) by estimated s.d(s) to gain SPM{t}(s)
	ResMS = spm_sample_vol(xSDM.VResMS,XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
	for i = 1:n+m
		Z(i,:) = Z(i,:)./(ResMS*(tmp(:,i)'*xX.Bcov*tmp(:,i)));
	end


	%-Apply any masks
	%===============================================================
	if m>0
		Q     = find(min(Z(n+1:n+m,:),[],1) > u);
		Z     = Z(1:n,Q);
		XYZ   = XYZ(:,Q);
		QQ    = QQ(Q);
	end

	%-Canonicalise variables for results section usage...
	%---------------------------------------------------------------
	edf = [1,xX.erdf];
	c   = num2cell(c, 1);
%	cm  = num2cell(cm,1);


else

	%-Get F-contrast matrix (contrasts as column vectors)
	%---------------------------------------------------------------
	c = spm_input('F-contrast','+1','x','',Inf,xX.xKXs)';
	n = 1;

	RB2           = spm_SpUtil('BetaRC',xX.xKXs,c);
	KX1           = spm_SpUtil('cTestSp',xX.xKXs,c);
	[trMV,trMVMV] = spm_SpUtil('trMV',KX1,xX.V);
	edf           = [trMV^2/trMVMV, xX.erdf];
	
	[u,s,v] = svd(RB2);
	RB      = u*sqrt(s)*v';	%-Residual (in parameter space) forming mtx
	

	%-Compute Extra sum-of-squares from parameter images
	%---------------------------------------------------------------
	VESS = struct(	'fname',	fullfile(swd,'tmpESS.img'),...
			'dim',		[xSDM.Vbeta(1).dim(1:3),16],...
			'mat',		xSDM.Vbeta(1).mat,...
			'pinfo',	[1,0,0]',...
			'descrip',	'spm_results: Extra sum-of-squares');
	VESS = spm_resss(xSDM.Vbeta,VESS,RB);

%	beta = zeros(nbeta,S);					%-****
%	for j = 1:nbeta						%-****
%		beta(j,:) = spm_sample_vol(xSDM.Vbeta(j),...	%-****
%				XYZ(1,:),XYZ(2,:),XYZ(3,:),0);	%-****
%	end							%-****
%	ESSc    = sum((RB*beta).^2)				%-****

	%-Compute SPM{F} from ESS and ResMS images
	%---------------------------------------------------------------
	Z = (spm_sample_vol(VESS,  XYZ(1,:),XYZ(2,:),XYZ(3,:),0) / trMV) ./ ...
	    (spm_sample_vol(xSDM.VResMS,XYZ(1,:),XYZ(2,:),XYZ(3,:),0) / xX.trRV);

	%-Canonicalise contrast storage (i.e. in cell array)
	%---------------------------------------------------------------
	c = {c};

end


%-Implement conjunction as minima of orthogonal contrasts
%-----------------------------------------------------------------------
Z = min(Z,[],1);



%=======================================================================
% - H E I G H T   &   E X T E N T   T H R E S H O L D S
%=======================================================================

%-Get and apply height threshold
%-----------------------------------------------------------------------
if spm_input('corrected height threshold','+1','y/n',[1,0],2)
	u  = spm_input('corrected p value','+0','e',0.05);
	u  = spm_U(u,edf,STAT,xSDM.R,n);
else
	%-NB: p-value imput won't work properly for conjunctions,
	%     since spm_u doesn't cope with minumum fields...
	u  = spm_input(['threshold {',STAT,' or p value}'],'+0','e',0.001);
	if u <= 1; u = spm_u(u,edf,STAT); end
end


%-Get extent threshold [default = 0]
%-----------------------------------------------------------------------
k     = spm_input('& extent threshold {resels}','+1','e',0);


% Eliminate voxels based on height
%-----------------------------------------------------------------------
Q     = find(Z > u);
Z     = Z(:,Q);
XYZ   = XYZ(:,Q);
QQ    = QQ(Q);


% ...return if there are no voxels
%-----------------------------------------------------------------------
if ~length(Q)
	error('No voxels above this threshold {u}')
end

%-Apply extent threshold
%-----------------------------------------------------------------------
A     = spm_clusters(XYZ);
Q     = [];
for i = 1:max(A)
	j = find(A==i);
	if length(j)*R(end)/S >= k; Q = [Q j]; end
end

% eliminate voxels
%-----------------------------------------------------------------------
Z     = Z(:,Q);
XYZ   = XYZ(:,Q);
QQ    = QQ(Q);


% ...and return if there are no clusters
%-----------------------------------------------------------------------
if ~length(Q)
	error('No clusters above this threshold {k}')
end

%=======================================================================
% - E N D
%=======================================================================

%-Compute mm<->voxel matrices
%-----------------------------------------------------------------------
 M  = xSDM.Vbeta(1).mat;
iM  = inv(M);


%-Assemble output structures
%=======================================================================
SPM    = struct('swd',		swd,...
		'c',		{c},...
		'Z',		Z,...
		'n',		n,...
		'STAT',		STAT,...
		'df',		edf,...
		'u',		u,...
		'k',		k,...
		'XYZ',		XYZ,...
		'XYZmm',	M(1:3,:)*[XYZ; ones(1,size(XYZ,2))],...
		'QQ',		QQ);

VOL    = struct('S',		S,...
		'R',		R,...
		'FWHM',		xSDM.FWHM,...
		'M',		M,...
		'iM',		iM,...
		'VOX',		sqrt(sum(M(1:3,1:3).^2))',...
		'DIM',		xSDM.Vbeta(1).dim(1:3)');

%-Remove (large) duplicate fields from xSDM
%-----------------------------------------------------------------------
xSDM = rmfield(xSDM,{'xX','XYZ'});
