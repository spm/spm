function [SPM,VOL,DES] = spm_getSPM
% computs a specified and thresholded following parameter estimation
% FORMAT [SPM,VOL,DES] = spm_getSPM;
%
% SPM    = struct('Z'		minimum of n Statistics {filtered on u and k}
% 		  'n'		number of conjoint tests	
% 		  'STAT'	distribution {Z, T, X or F}	
% 		  'df'		[df{interest} df{error}]
%		  'u'		height threshold
% 		  'k'		extent threshold {resels}
% 
% VOL    = struct('R'		Search Volume {resels}
% 		  'S'		search Volume {voxels}
% 		  'FWHM'	Smoothness {voxels}	
% 		  'DIM'		image dimensions {voxels}
% 		  'VOX'		voxel dimensions {mm}
% 		  'XYZ'		location of voxels (mm}	
% 		  'ORG'		origin x,y,z  = 0 {voxels}	
% 		  'M'		voxels - > mm matrix	
% 		  'QQ'		indices of volxes in ????.mat files	
% 
% DES    = struct('X'		Design Matrix
% 		  'C'		Contrast weights	
% 		  'B'		pinv(X)*V*pinv(X)'/trace(R*V)
%_______________________________________________________________________
%
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
% The SPM is subject to thresholding on the basis of height (u) and
% the number of resels comprising its clusters {k}. The height threshold
% is specified as above in terms of an [un]corrected p value or statistic.
% Clusters can also be thresholded on the basis of their spatial extent.
% SPM98 expresses all cluster volumes in terms of resels.  These are effective
% resolution elements and comprise several voxels, depending on the smoothness
% estimator. 
% If you want to see all voxels simply enter 0.  In this instance the
% 'set-level' inference can be considered an 'omnibus test' based on the
% number of clusters that obtain.
% 
% see spm_list.m and spm_P.m for further details
%
%__________________________________________________________________________
% %W% Karl Friston %E%

%-Defalt Action
%--------------------------------------------------------------------------
global CWD

%-Get figure handles and filenames
%--------------------------------------------------------------------------
Finter  = spm_figure('FindWin','Interactive');
spm_clf(Finter)

% What sort of SPM
%---------------------------------------------------------------------------
STAT    = spm_input('which SPM',1,...
		    'b','SPM{T}|SPM{F}',['T' 'F']);

set(Finter,'Name',['SPM{' STAT '} projections'])

str     = ['SPM.mat'];
tmp     = spm_get(1,str,['select ' str]);
CWD     = strrep(tmp,['/' str],'');

%-Get data
%--------------------------------------------------------------------------
load([CWD,'/SPM'])
load([CWD,'/XYZ'])
QQ         = 1:size(XYZ,2);

% backwards compatibility
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DES        = [H C B G];

% re-compute smoothness of component fields
%---------------------------------------------------------------------------
W          = W*sqrt(spm_lambda(Fdf(2)));
FWHM       = W*sqrt(8*log(2));
R          = spm_resels(FWHM,S);
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


% Get contrasts for SPM{Z} and SPM{t} and apply masks
%===========================================================================
if STAT == 'T'


	%-Get contrast[s]
	%-------------------------------------------------------------------
	m      = size(DES,2);
	str    = sprintf('contrast[s] n x 1 - %i',m);
	c      = spm_input(str,1);
	n      = size(c,1);

	%-Zero pad
	%-------------------------------------------------------------------
	c      = [c zeros(n,m - size(c,2))];

	%-Enforce orthogonality {within subspace spanned by contrasts}
	%-------------------------------------------------------------------
	p      = c;
	c      = c(1,:);
	for  i = 2:n
		j      = find(p(i,:));
		d      = p(i,j);
		q      = BCOV(j,j)*c(:,j)';
		d      = d' - q*(pinv(q)*d');
		c(i,j) = d';
	end

	% get parameter estimates and contruct new partitions
	%-------------------------------------------------------------------
	load([CWD,'/BETA'])
	load([CWD,'/RES'])


	% compute SPM
	%-------------------------------------------------------------------
	Z      = zeros(n,length(RES));
	for  i = 1:n
		Z(i,:) = c(i,:)*BETA./sqrt(RES*(c(i,:)*BCOV*c(i,:)'));
	end
	CON    = c;


	%-Get and apply any masks
	%===================================================================
	if spm_input('mask with other contrast[s]','!+1','b','no|yes',[0 1],1)

		%-Get contrast[s] for mask
		%-----------------------------------------------------------
		str   = sprintf('contrast[s] n x 1 - %i',m);
		c     = spm_input(str,1);
		n     = size(c,1);

		%-Zero pad
		%-----------------------------------------------------------
		c     = [c zeros(n,m - size(c,2))];


		% threshold for mask
		%-----------------------------------------------------------
		u     = spm_input('threshold for mask','!+1','e',0.05);
		if u <= 1; u = spm_u(u,Fdf,'T'); end

		% compute mask
		%-----------------------------------------------------------
		M     = zeros(n,length(RES));
		for i = 1:n
			M(i,:) = c(i,:)*BETA./sqrt(RES*(c(i,:)*BCOV*c(i,:)'));
		end

		% eliminate voxels
		%-----------------------------------------------------------
		Q     = find(min(M,[],1) > u);
		Z     = Z(:,Q);
		XYZ   = XYZ(:,Q);
		QQ    = QQ(Q);
	end
else

	load([CWD,'/SPMF'])
	Z     = SPMF;
	CON   = [];
end

%-Get and apply height threshold
%---------------------------------------------------------------------------
n     = size(Z,1);
Z     = min(Z,[],1);

% get height and extent thresholds
%===========================================================================


%-Get and apply height threshold
%---------------------------------------------------------------------------
if spm_input('corrected height threshold','!+1','b','no|yes',[0 1],1)
	u  = spm_input('corrected p value','!+0','e',0.05);
	u  = spm_U(u,Fdf,STAT,R,n);
else
	%-Get and apply height threshold [default p < 0.001 uncorrected]
	%-------------------------------------------------------------------
	u  = spm_input('threshold {F or p value}','!+0','e',0.001);
	if u <= 1; u = spm_u(u,Fdf,STAT); end
end


%-Get extent threshold [default = 0]
%---------------------------------------------------------------------------
k     = spm_input('& extent threshold {resels}','!+1','e',0);


% eliminate voxels based on height
%---------------------------------------------------------------------------
Q     = find(Z > u);
Z     = Z(:,Q);
XYZ   = XYZ(:,Q);
QQ    = QQ(Q);


%-Return if there are no voxels
%---------------------------------------------------------------------------
if ~length(Q)
	error('No voxels above this threshold {u}')
end

% and apply extent threshold
%---------------------------------------------------------------------------
A     = spm_clusters(XYZ,V([4 5 6]));
Q     = [];
for i = 1:max(A)
	j = find(A == i);
	if length(j)*R(length(R))/S >= k; Q = [Q j]; end
end

% eliminate voxels
%---------------------------------------------------------------------------
Z     = Z(:,Q);
XYZ   = XYZ(:,Q);
QQ    = QQ(Q);


%-Return if there are no clusters
%---------------------------------------------------------------------------
if ~length(Q)
	error('No clusters above this threshold {k}')
end

%-Finished
%---------------------------------------------------------------------------
spm_clf(Finter)
set(Finter,'Name',' ','Pointer','Arrow')


% assemble output structures
%===========================================================================
SPM    = struct('Z',	Z,...
		'n',	n,...
		'STAT',	STAT,...
		'df',	Fdf,...
		'u',	u,...
		'k',	k);

VOL    = struct('R',	R,...
		'FWHM',	FWHM,...
		'S',	S,...
		'DIM',	V(1:3),...
		'VOX',	V(4:6),...
		'ORG',	V(7:9),...
		'M',    [diag(V(4:6)) -V(7:9).*V(4:6); 0 0 0 1],...
		'XYZ',	XYZ,...
		'QQ',	QQ);

DES    = struct('X',	DES,...
		'C',	CON,...
		'B',	BCOV);

