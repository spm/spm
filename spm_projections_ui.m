function [Z,XYZ,QQ,u,k,S,W] = spm_projections_ui(Action)
% used to review results of statistical analysis (SPM{Z})
% FORMAT [Z,XYZ,QQ,u,k,S,W] = spm_projections_ui(Action)
%
% Action - 'Display'  - Calls spm_projections
%        - 'Results'  - Just returns output variables
%        - 'Writing'  - writes filtered SPM{t}
%
% Z      - Z values after filtering on height and size thresholds
% XYZ    - location in mm
% QQ     - Indices of selected voxels
% u      - selected height threshold
% k      - selected extent threshold {voxels}
% S      - search volume {voxels}
% W      - smoothness estimators (of Gaussianized t fields)
%_______________________________________________________________________
%
% 
% spm_projections_ui allows the SPM{Z} created by spm_spm.m to be
% re-displayed and characterized in terms of regionally significant
% effects.
% 
% Either single contrasts can be examined or conjunctions of different
% contrasts.  In the latter case a new SPM{Z} is created, that reflects
% the sum of all specified effects and voxels are eliminated where there
% are significant differences among these effects (at p < 0.05
% uncorrected).  A conjunction is therefore the conjoint expression of
% two or more effects (each specified in terms of a contrast) to a
% similar degree.
%
% Masking simply eliminates voxels from the current contrast if they
% do not survive an uncorrected p value (based on height) in one or
% more further contrasts.  No account is taken of this masking in the
% statistical inference pertaining to the masked contrast.
% 
% The resulting SPM{Z} can be masked by another (corresponding to another
% contrast).  This masking is not taken into account when presenting Z
% scores or corrected staistical inference.
% 
% The SPM{Z} is subject to thresholding on the basis of height (Z) and
% the number of voxels comprising its clusters {k}. The height threshold
% is specified as before in terms of an [un]corrected p value or Z score.
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
%_______________________________________________________________________
% %W% Karl Friston %E%

%-Defalt Action
%-----------------------------------------------------------------------
if nargin == 0, Action = 'Display'; end
global CWD

%-Get figure handles and filenames
%-----------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
Fgraph = spm_figure('FindWin','Graphics');
spm_clf(Finter)
set(Finter,'Name','SPM{Z} projections')

tmp = spm_get(1,'.mat','select SPMt.mat for analysis','SPMt');
CWD = strrep(tmp,'/SPMt.mat','');

%-Get data
%-----------------------------------------------------------------------
load([CWD,'/SPM'])
load([CWD,'/XYZ'])
load([CWD,'/SPMt'])
load([CWD,'/RES']);

QQ = (1:size(XYZ,2))';


%-Get contrast[s]
%-----------------------------------------------------------------------
DES  = [H C B G]; 	% design matrix
i    = 0;
if spm_input('Conjunction analysis',1,'b','no|yes',[0 1],1)

	% conjunction analysis
	%---------------------------------------------------------------
	while any(i < 1 | i > size(CONTRAST,1))
		str = sprintf('contrasts ? 1 - %i',size(CONTRAST,1));
		i   = spm_input(str,1);
		CON = CONTRAST(i,:);
	end

	% correlation between the Z scores under the null hypothesis
	%---------------------------------------------------------------
	J     = CON*BCOV*CON';
	J     = inv(diag(sqrt(diag(J))))'*J*inv(diag(sqrt(diag(J))));

	% transform
	%---------------------------------------------------------------
	[e v] = eig(J);
	T     = e*inv(sqrt(v));
	Z     = SPMt(i,:)'*T;
	o     = inv(T)*ones(length(i),1);
	r     = sum(Z'.^2);
	Z     = Z*o/sqrt(o'*o);

	% eliminate voxels with significant deviation from a conjunction
	%---------------------------------------------------------------
	r     = r - Z'.^2;
	U     = spm_invXcdf((1 - 0.05),length(i) - 1);
	Q     = find(r' < U);
	S     = S - (length(Z) - length(Q)) ;
	Z     = Z(Q);
	RES   = RES(Q);
	XYZ   = XYZ(:,Q);
	SPMt  = SPMt(:,Q);
	QQ    = QQ(Q);

else
	% straightforward contrast
	%---------------------------------------------------------------
	if size(CONTRAST,1) > 1
			while any(i < 1 | i > size(CONTRAST,1))
			str = sprintf('contrast ? 1 - %i',size(CONTRAST,1));
			i   = spm_input(str,1);
			CON = CONTRAST(i(1),:);
		end
	else
		i = 1;
	end
	Z     = SPMt(i,:)';

end


%-Get and apply any masks
%-----------------------------------------------------------------------
I     = i;
i     = 0;
if size(CONTRAST,1) > 1
    if spm_input('mask with other contrast[s]','!+1','b','no|yes',[0 1],1)
	while any(i < 1 | i > size(CONTRAST,1))
		str = sprintf('contrasts ? 1 - %i',size(CONTRAST,1));
		i   = spm_input(str,'!+0');
	end

	% threshold for mask
	%---------------------------------------------------------------
	u     = spm_input('threshold for mask','!+1','e',0.05);
	if u < 1; u = spm_invNcdf(1 - u); end
	if length(i) > 1
		Q   = find(all(SPMt(i,:) > u));
	else
		Q   = find(SPMt(i,:) > u);
	end

	% eliminate voxels
	%---------------------------------------------------------------
	Z     = Z(Q);
	XYZ   = XYZ(:,Q);
	QQ    = QQ(Q);
    end
end


%-Get and apply height threshold
%-----------------------------------------------------------------------
if spm_input('use corrected height threshold','!+1','b','no|yes',[0 1],1)
	u  = spm_input('corrected p value','!+0','e',0.05);
	u  = spm_z(u,W,S);
else
	%-Get and apply height threshold [default p < 0.001 uncorrected]
	%---------------------------------------------------------------
	u  = spm_input('height threshold {Z or p value}','!+0','e',0.001);
	if u < 1; u = spm_invNcdf(1 - u); end
end


% eliminate voxels
%-----------------------------------------------------------------------
Q     = find(Z > u);
Z     = Z(Q);
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

%-Get extent threshold
%-----------------------------------------------------------------------
if spm_input('use corrected extent threshold','!+1','b','no|yes',[0 1],1)
	k  = spm_input('corrected p value','!+0','e',0.05);
	k  = spm_k(k,W,u,S);
else
	%-Get and apply extent threshold [default p < 0.5 uncorrected]
	%---------------------------------------------------------------
	k  = spm_input('extent threshold {k or p value}','!+0','e',0.5);
	if (k < 1) & (k > 0); k = spm_invkcdf(1 - k,u,W); end
end

% and apply
%-----------------------------------------------------------------------
A     = spm_clusters(XYZ,V([4 5 6]));
Q     = [];
for i = 1:max(A)
	j = find(A == i);
	if length(j) >= k
		Q = [Q j];
	end
end

% eliminate voxels
%-----------------------------------------------------------------------
Z     = Z(Q);
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
	spm_projections(Z,XYZ,u,k,V,W,S,DES,CON,df);

% Results
%-----------------------------------------------------------------------
elseif  strcmp(lower(Action),lower('Results'))

	spm_clf(Finter)
	set(Finter,'Name','Thankyou')


% Write SPM{Z}
%-----------------------------------------------------------------------
elseif  strcmp(lower(Action),lower('Writing'))

	FName   = spm_input('Filename ?',4,'s','SPM_filtered');
	str     = sprintf('spm{Z}-filtered: u = %5.3f, k = %d',u,k);

	%-Reconstruct filtered image from XYZ & t
	%---------------------------------------------------------------
	n       = size(XYZ,2);
	rcp     = round(XYZ./meshgrid([1 - 2*FLIP;1;1].*V(4:6),1:n)' + ...
		meshgrid(V(7:9),1:n)');
	Dim     = cumprod([1,V(1:2)']);
	OffSets = meshgrid([0,1,1],1:n)';
	e       = ((rcp - OffSets)'*Dim')';
	Z       = Z.*(Z > 0);
	T       = zeros(1,prod(V(1:3)));
	T(e)    = Z;

	%-Write out to analyze file
	%---------------------------------------------------------------
	spm_hwrite([FName,'.hdr'],V(1:3),V(4:6),1/16,2,0,V(7:9),str);
	fid     = fopen([FName,'.img'],'w');
	fwrite(fid,T*16,spm_type(2));
	fclose(fid);

% Unknwn Action
%-----------------------------------------------------------------------
else
	error('Unknown action string')

end

%-Finished
%-----------------------------------------------------------------------
spm_clf(Finter)
set(Finter,'Name',' ','Pointer','Arrow')

return
