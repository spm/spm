function [Z,XYZ,QQ,u,k] = spm_projections_ui(Action,STAT)
% used to review results of statistical analysis (SPM{Z})
% FORMAT [Z,XYZ,QQ,u,k] = spm_projections_ui(Action,STAT)
%
% Action - 'Display'  - Calls spm_projections
%        - 'Results'  - Just returns output variables
% STAT   - Statistic
%
% Z      - Statistic values after filtering on height and size thresholds
% XYZ    - location in mm
% QQ     - Indices of selected voxels
% u      - selected height threshold
% k      - selected extent threshold {resels}
%_______________________________________________________________________
%
% 
% spm_projections_ui prompts for an SPM and applies thresholds {u & k}
% to a point list of voxel values (specified with their locations {XYZ})
% This allows the SPM be displayed and characterized in terms of regionally 
% significant effects.
% 
% Either single contrasts can be examined or conjunctions of different
% contrasts.  In the latter case a new SPM is created, that reflects
% the miminum of all specified effects.  A conjunction is therefore
% the conjoint expression of two or more effects.
%
% Masking simply eliminates voxels from the current contrast if they
% do not survive an uncorrected p value (based on height) in one or
% more further contrasts.  No account is taken of this masking in the
% statistical inference pertaining to the masked contrast.
% 
% The SPM is subject to thresholding on the basis of height (u) and
% the number of voxels comprising its clusters {k}. The height threshold
% is specified as before in terms of an [un]corrected p value or statistic.
% If you only want to see clusters that survive a corrected p value
% (based on spatial extent) than the corrected p value you enter will
% specify the extent threshold employed.  If however you choose to see
% all clusters you can specify an extent threshold in terms of an
% uncorrected p value (i.e. the probability of getting a cluster that
% size or larger, assuming it exists) or voxels.  If you want to see all
% voxels simply enter 0.  In this instance the 'set-level' inference can
% be considered an 'omnibus test' based on the number of clusters that
% obtain.
% 
% see spm_projections.m and spm_P for further details
%
%__________________________________________________________________________
% %W% Karl Friston %E%

%-Defalt Action
%--------------------------------------------------------------------------
global CWD

%-Get figure handles and filenames
%--------------------------------------------------------------------------
Finter  = spm_figure('FindWin','Interactive');
Fgraph  = spm_figure('FindWin','Graphics');
spm_clf(Finter)
set(Finter,'Name',['SPM{' STAT '} projections'])

str     = ['SPM.mat'];
tmp     = spm_get(1,str,['select ' str]);
CWD     = strrep(tmp,['/' str],'');

%-Get data
%--------------------------------------------------------------------------
load([CWD,'/SPM'])
load([CWD,'/XYZ'])
QQ         = 1:size(XYZ,2);
DES        = [H C B G];

% backwards compatibility
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% re-compute smoothness of component fields
%---------------------------------------------------------------------------
if STAT ~= 'Z'
	W  = W*sqrt(spm_lambda(Fdf(2)));
end
FWHM       = W*sqrt(8*log(2));
R          = spm_resels(FWHM,S);
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


% Get contrasts for SPM{Z} and SPM{t} and apply masks
%===========================================================================
if STAT == 'Z' | STAT == 'T'


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

	% Gaussianization
	%-------------------------------------------------------------------
	if STAT == 'Z'
		Z     = spm_t2z(Z,Fdf(2));
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


% get height and extent thresholds
%===========================================================================


%-Get and apply height threshold
%-----------------------------------------------------------------------
if spm_input('corrected height threshold','!+1','b','no|yes',[0 1],1)
	u  = spm_input('corrected p value','!+0','e',0.05);
	u  = spm_U(u,Fdf,STAT,R,size(Z,1));
else
	%-Get and apply height threshold [default p < 0.001 uncorrected]
	%---------------------------------------------------------------
	u  = spm_input('threshold {F or p value}','!+0','e',0.001);
	if u <= 1; u = spm_u(u,Fdf,STAT); end
end


%-Get extent threshold [default = 0]
%-----------------------------------------------------------------------
k     = spm_input('& extent threshold {resels}','!+1','e',0);


% eliminate voxels based on height
%-----------------------------------------------------------------------
Q     = find(min(Z,[],1) > u);
Z     = Z(:,Q);
XYZ   = XYZ(:,Q);
QQ    = QQ(Q);


%-Return if there are no voxels
%-----------------------------------------------------------------------
if ~length(Q)
	figure(Fgraph); spm_clf
	axis off
	text(0,0.3,spm_str_manip(CWD,'a50'));
	text(0,0.2,'No voxels above this threshold {u}','FontSize',16);
	return
end

% and apply extent threshold
%-----------------------------------------------------------------------
A     = spm_clusters(XYZ,V([4 5 6]));
Q     = [];
for i = 1:max(A)
	j = find(A == i);
	if length(j)*R(length(R))/S >= k; Q = [Q j]; end
end

% eliminate voxels
%-----------------------------------------------------------------------
Z     = Z(:,Q);
XYZ   = XYZ(:,Q);
QQ    = QQ(Q);


%-Return if there are no clusters
%-----------------------------------------------------------------------
if ~length(Q)
	figure(Fgraph); spm_clf
	axis off
	text(0,0.3,spm_str_manip(CWD,'a50'));
	text(0,0.2,'No clusters above this threshold {k}','FontSize',16);
	return
end


%-Switch depending on Action
%=======================================================================

%-Display and characterize SPM{Z}
%-----------------------------------------------------------------------
if      strcmp(lower(Action),lower('Display'))

	set(Finter,'Name','Thankyou','Pointer','watch')
	spm_projections(Z,Fdf,STAT,R,S,FWHM,u,k*max(R)/S,XYZ,V,DES,CON)

% Results
%-----------------------------------------------------------------------
elseif  strcmp(lower(Action),lower('Results'))

	spm_clf(Finter)
	set(Finter,'Name','Thankyou')

end

%-Finished
%-----------------------------------------------------------------------
spm_clf(Finter)
set(Finter,'Name',' ','Pointer','Arrow')

return
