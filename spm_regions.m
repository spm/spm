function [Y,xY] = spm_regions(SPM,VOL,xX,xCon,xSDM,hReg)
% VOI time-series extraction of adjusted data (local eigenimage analysis)
% FORMAT [Y xY] = spm_regions(SPM,VOL,xX,xCon,xSDM,hReg);
%
% SPM    - structure containing SPM, distribution & filtering detals
% VOL    - structure containing details of volume analysed
% xX     - Design Matrix structure
% xSDM   - structure containing contents of SPM.mat file
% xCon   - Contrast definitions structure (see spm_FcUtil.m for details)
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% Y      - first eigenvariate of VOI
% xY     - structure with:
%       xY.name         - name of VOI
%       xY.y            - voxel-wise data (filtered and adjusted)
%       xY.u            - first eigenvariate
%       xY.v            - first eigenimage
%       xY.s            - eigenimages
%       xY.XYZmm        - Co-ordinates of voxels used within VOI
%       xY.xyz          - centre of VOI (mm)
%       xY.radius       - radius of VOI (mm)
%       xY.dstr         - description of filtering & adjustment applied
%
% Y and xY are also saved in VOI_*.mat in the SPM working directory
%
% (See spm_getSPM for details on the SPM,VOL, xX & xSDM structures.)
%
%_______________________________________________________________________
%
% spm_regions extracts a representative time course from voxel data in
% Y.mad in terms of the first eigenvariate of the filtered abd adjusted
% data in all suprathreshold voxels saved (in Y.mad) within a spherical VOI
% centered on the nearest Y.mad voxel to the selected location.
%
% If temporal filtering has been specified (fMRI), then the data is
% temporally filtered. Adjustment is with respect to the null space of
% a selected contrast, or can be omitted.
%
% For a VOI of radius 0, the (filtered &) adjusted voxel time-series is
% returned, scaled to have a 2-norm or 1. The actual (filtered &)
% adjusted voxel time series can be extracted from xY.y, and will be
% the same as the (filtered &) adjusted data returned by the plotting
% routine (spm_graph.m) for the same contrast.
%
% See spm_spm.m for further details of how voxels are selected for the
% saving of their raw data in Y.mad.
%_______________________________________________________________________
% %W% Karl Friston %E%


% get figure handles
%-----------------------------------------------------------------------
Finter = spm_figure('GetWin','Interactive');
Fgraph = spm_figure('GetWin','Graphics');
set(Finter,'Name','VOI time-series extraction')


%-Find nearest voxel [Euclidean distance] in point list with data saved
% in Y.mad, and update GUI location
%-----------------------------------------------------------------------
if ~length(SPM.XYZmm)
	spm('alert!','No suprathreshold voxels!',mfilename,0);
	Y = []; xY = [];
	return
elseif exist(fullfile(SPM.swd,'Y.mad')) ~= 2
	spm('alert!','No raw data saved with this analysis!',mfilename,0);
	Y = []; xY = [];
	return
elseif ~any(SPM.QQ)
	spm('alert!','No raw data saved for any suprathreshold location!',...
		mfilename,0);
	Y = []; xY = [];
	return
end


[xyz,i] = spm_XYZreg('NearestXYZ',...
		spm_XYZreg('GetCoords',hReg),SPM.XYZmm(:,find(SPM.QQ)));
spm_XYZreg('SetCoords',xyz,hReg);


%-Get adjustment options, VOI name and radius
%-----------------------------------------------------------------------
spm_input(sprintf('at [%3.0f %3.0f %3.0f]',xyz),1,'d',...
	'VOI time-series extraction')
Ic      = spm_input('Adjust data for (select contrast)',2,'m',...
		{'<don''t adjust>',xCon.name})-1;
Rname   = spm_input('name of region',3,'s');
R       = spm_input('VOI radius (mm)',4,'r',0,1,[0,Inf]);

%-Find suprathreshold voxels within radius, & note those also in Y.mad
%-----------------------------------------------------------------------
d       = [SPM.XYZmm(1,:)-xyz(1); SPM.XYZmm(2,:)-xyz(2); SPM.XYZmm(3,:)-xyz(3)];
Q       = find(sum(d.^2) <= R^2);
q       = find(SPM.QQ(Q));

if any(SPM.QQ(Q)==0)
    spm('alert"',{...
        sprintf('Don''t have raw data for all %d suprathreshold',length(Q)),...
	sprintf('voxels within %.1gmm of [%3.0f %3.0f %3.0f]',R,xyz),' ',...
	sprintf('Proceeding using the %d voxels that are.',length(q)),...
	},mfilename,sqrt(-1));
end


%-Extract required data from results files
%=======================================================================

%-Get (approximate) raw data y from Y.mad file
%-NB: Data in Y.mad file is compressed, and therefore not fully accurate
%-----------------------------------------------------------------------
y       = spm_extract(fullfile(SPM.swd,'Y.mad'),SPM.QQ(Q(q)));
rcp     = VOL.iM(1:3,:)*[SPM.XYZmm(:,Q(q));ones(size(q))];

%-Parameter estimates: beta = xX.pKX*xX.K*y; (load from file for accuracy)
%-----------------------------------------------------------------------
nBeta = length(xSDM.Vbeta);
beta  = ones(nBeta,size(q,2));
for   i = 1:nBeta
	beta(i,:) = spm_sample_vol(xSDM.Vbeta(i),rcp(1,:),rcp(2,:),rcp(3,:),0);
end


%-Computation
%=======================================================================

% filter and remove confounds
%-----------------------------------------------------------------------
y    = spm_filter('apply',xX.K, y);

if Ic
	y = y - spm_FcUtil('Y0',xCon(Ic),xX.xKXs,beta);
	dstr = ['adjusted for ',xCon(Ic).name];
else
	dstr = 'not adjusted';
end

% compute regional response in terms of first eigenvariate
%-----------------------------------------------------------------------
[m n]   = size(y);
if m > n
	[v s v] = svd(spm_atranspa(y));
	s       = diag(s);
	v       = v(:,1);
	u       = y*v/sqrt(s(1));
else
	[u s u] = svd(spm_atranspa(y'));
	s       = diag(s);
	u       = u(:,1);
	v       = y'*u/sqrt(s(1));
end
d       = sign(sum(v));
u       = u*d;
v       = v*d;
Y       = u;


%-Display MIP of VOI weighting and timecourse
% NB: timecourse axes overlap MIP: This obscures the coronal MIP (for MIP96)
%     in 3D mode, but for 2D data half the plane image will be obscured.
%-----------------------------------------------------------------------
spm_results_ui('Clear',Fgraph);
figure(Fgraph);
axes('Position',[0.0500,0.1100,0.5500,0.3388])
spm_mip(v(:,1),SPM.XYZmm(:,Q(q)),VOL.M,VOL.DIM)
title('VOI weighting                                          ')

axes('Position',[0.376,0.130,0.520,0.294])
plot(Y)
set(gca,'YAxisLocation','Right')
title(['1st eigenvariate: ' Rname],'FontSize',16)
str = {	'scan number';' ';sprintf(...
	'%d voxels in sphere of radius %.1gmm at [%3.0f %3.0f %3.0f]',...
	length(q),R,xyz);dstr;sprintf('Variance: %0.2f%%',s(1)*100/sum(s))};
xlabel(str)

% create structure
%-----------------------------------------------------------------------
xY      = struct('name',	Rname,...
		 'y',		y,...
		 'u',		u,...
		 'v',		v,...
		 's',		s,...
		 'XYZmm',	SPM.XYZmm(:,Q(q)),...
		 'xyz',		xyz,...
		 'radius',	R,...
		 'dstr',	dstr);
% save
%-----------------------------------------------------------------------
save(fullfile(SPM.swd,['VOI_',Rname]),'Y','xY')

%-Reset title
%-----------------------------------------------------------------------
spm('FigName',['SPM{',SPM.STAT,'}: Results']);
