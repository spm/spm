function [Z,XYZ,XA,u,k,S,W] = spm_projectionsF_ui(Action,Fname)
% used to review results of statistical analysis (SPM{Z})
% FORMAT [F,XYZ,XA,u,k,S,W] = spm_projectionsF_ui(Action,Fname)
%
% Action - 'Display'  - Calls spm_projections
%        - 'Results'  - Just returns output variables
%        - 'Writing'  - writes filtered SPM{F} to Fname
%
% F      - F values after filtering on height and size thresholds
% XYZ    - location in mm
% XA     - Adjusted activities ('Results' option only)
% u      - selected height threshold
% k      - selected extent threshold {voxels}
% S      - search volume {voxels}
% W      - smoothness estimators (of component fields)
%_______________________________________________________________________
%
% spm_projectionsF_ui allows the SPM{F} created by spm_spm.m to be
% re-displayed and characterized in terms of regionally significant
% effects. This is based on K. Worsley results for the expected maximum
% value in an F-field and requires the smoothness estimation of the
% original component fields.
%
%
% see spm_projectionsF.m for further details
%
%_______________________________________________________________________
% %W% JBP, Karl Friston %E%


%-Defalt Action
%-----------------------------------------------------------------------
if nargin == 0, Action = 'Display'; end
global CWD

%-Get figure handles and filenames
%-----------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
Fgraph = spm_figure('FindWin','Graphics');
spm_clf(Finter)
set(Finter,'Name','SPM{F} projections')

tmp = spm_get(1,'.mat','select SPMF.mat for analysis','SPMF');
CWD = strrep(tmp,'/SPMF.mat','');
K   = [];
XA  = [];

%-Get data
%-----------------------------------------------------------------------
load([CWD,'/SPM'])
load([CWD,'/XYZ'])
load([CWD,'/SPMF'])
if strcmp(lower(Action),lower('Results'))
	load([CWD,'/XA']); end

%-Design matrix
%-----------------------------------------------------------------------
DES   = [K H C B G];
Z     = SPMF;


% compute smoothness of component fields
%-----------------------------------------------------------------------
[Lc2z Lc2t Lt2z] = spm_lambda(df);
W     = W*sqrt(Lc2z);


%-Get and apply height threshold [default p < 0.001]
%-----------------------------------------------------------------------
u     = spm_input('height threshold for SPM{F}',1,'e',0.001);
if u < 1; u = spm_invFcdf(1 - u,Fdf); end

% eliminate voxels
%-----------------------------------------------------------------------
Q     = find(Z > u);
Z     = Z(Q);
XYZ   = XYZ(:,Q);
if strcmp(lower(Action),lower('Results'))
	XA = XA(:,Q); end

%-Return if there are no voxels
%-----------------------------------------------------------------------
if ~length(Q)
	axis off
	text(0,0.3,spm('DirTrunc',CWD));
	text(0,0.2,'No voxels above this threshold {u}','FontSize',16);
	return
end

%-Get and apply extent threshold [default = E{n}]
%-----------------------------------------------------------------------
k     = spm_input('& extent threshold {voxels}',2,'e',round(prod(W)));
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
if strcmp(lower(Action),lower('Results'))
	XA = XA(:,Q);
end


%-Return if there are no clusters
%-----------------------------------------------------------------------
if ~length(Q)
	axis off
	text(0,0.3,spm('DirTrunc',CWD));
	text(0,0.2,'No clusters above this threshold {k}','FontSize',16);
	return
end

%-Switch depending on Action
%=======================================================================

%-Display and characterize SPM{Z}
%-----------------------------------------------------------------------
if      strcmp(lower(Action),lower('Display'))

	set(Finter,'Name','Thankyou','Pointer','watch')
	spm_projections(Z,XYZ,u,k,V,W,S,DES,[],Fdf);

% Results
%-----------------------------------------------------------------------
elseif  strcmp(lower(Action),lower('Results'))

	spm_clf(Finter)
	set(Finter,'Name','Thankyou')


% Write SPM{Z}
%-----------------------------------------------------------------------
elseif  strcmp(lower(Action),lower('Writing'))

	if nargin < 2
		FName = spm_input('Filename ?',3,'s','SPM_filtered');
	else
		FName = P2;
	end
	str = sprintf('spm{F}-filtered: u = %5.3f, k = %d',u,k);

	%-Reconstruct filtered image from XYZ & t
	%---------------------------------------------------------------
	n       = size(XYZ,2);
	rcp     = round(XYZ./meshgrid([1 - 2*FLIP;1;1].*V(4:6),1:n)' + ...
		meshgrid(V(7:9),1:n)');
	Dim     = cumprod([1,V(1:2)']);
	OffSets = meshgrid([0,1,1],1:n)';
	e       = ((rcp - OffSets)'*Dim')';
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
