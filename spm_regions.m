function [Y,xY] = spm_regions(SPM,VOL,xX,xSDM,hReg)
% VOI time-series extraction of adjusted data (local eigenimage analysis)
% FORMAT [Y xY] = spm_regions(SPM,VOL,xX,xSDM,hReg);
%
% SPM    - structure containing SPM, distribution & filtering detals
% VOL    - structure containing details of volume analysed
% xX     - Design Matrix structure
% xSDM   - structure containing contents of SPM.mat file
% hReg   - handle of MIP register
%
% Y      - first eigenvariate of VOI
% xY     - structure with:
% 	xY.name 	- name of VOI
% 	xY.Y 		- first eigenvariate
% 	xY.y 		- voxel-wise data (filtered and adjusted)
% 	xY.XYZ 		- centre of VOI (mm)
% 	xY.radius 	- radius of VOI (mm)
%
% Y and xY are also saved in VOI_???.mat in the SPM directory
%___________________________________________________________________________
% spm_regions extracts a representative time course from voxel data in
% Y.mad in terms of the first eigenvariate of filtered and adjusted [for
% confounds] data in all suprathreshold voxels saved within a spherical
% VOI centered on the nearest Y.mad voxel to the selected location.
%---------------------------------------------------------------------------
% %W% Karl Friston %E%


% get figure handles
%---------------------------------------------------------------------------
Finter = spm_figure('GetWin','Interactive');
Fgraph = spm_figure('GetWin','Graphics');
set(Finter,'Name','VOI time-series extraction')

%-Find nearest voxel [Euclidean distance] in point list & update GUI
%---------------------------------------------------------------------------
if ~length(SPM.XYZmm)
	spm('alert!','No suprathreshold voxels!',mfilename,0);
	Y = []; xY = [];
	return
end

% get VOI and name
%---------------------------------------------------------------------------
Rname   = spm_input('name of region',1,'s');
R       = spm_input('VOI radius (mm)',2,'e',0);

% get selected location
%---------------------------------------------------------------------------
Q       = find(SPM.QQ);
XYZ     = SPM.XYZmm(:,Q);
[L,i]   = spm_XYZreg('NearestXYZ',spm_XYZreg('GetCoords',hReg),XYZ);
spm_XYZreg('SetCoords',L,hReg);


% find voxels within radius
%---------------------------------------------------------------------------
d       = [XYZ(1,:)-L(1); XYZ(2,:)-L(2); XYZ(3,:)-L(3)];
q       = find(sum(d.^2) <= R^2);
y       = spm_extract('Y.mad',SPM.QQ(Q(q)));
rcp     = VOL.iM(1:3,:)*[XYZ(:,q);ones(size(q))];

%-Parameter estimates: beta = xX.pKX*xX.K*y;
%---------------------------------------------------------------------------
if isstruct(xSDM.F_iX0)
	Ic      = xSDM.F_iX0(1).iX0;
else
	Ic      = xSDM.F_iX0;
end
beta    = ones(length(Ic),size(q,2));
for   i = 1:length(Ic)
	beta(i,:) = ...
	spm_sample_vol(xSDM.Vbeta(Ic(i)),rcp(1,:),rcp(2,:),rcp(3,:),0);
end

% remove confounds and filter
%---------------------------------------------------------------------------
y       = spm_filter('apply',xX.K, y) - xX.xKXs.X(:,Ic)*beta;

% compute regional response in terms of first eigenvariate
%---------------------------------------------------------------------------
[u s v] = svd(y,0);
d       = sign(mean(u(:,1)));
u       = u*d;
v       = v*d;
Y       = u(:,1);
s       = diag(s).^2;


% display MIP of VOI and timecourse
%---------------------------------------------------------------------------
spm_results_ui('Clear',Fgraph,2);
figure(Fgraph);
axes('Position',[0.0500    0.1100    0.5500    0.3388])
spm_mip(v(:,1),XYZ(:,q),VOL.M,VOL.DIM)
title('VOI weighting                                          ')

axes('Position',[0.376    0.13    0.5200    0.294])
plot(Y)
title(['1st eigenvariate: ' Rname],'FontSize',16)
str      = {	'scan',...
		sprintf('Voxels: %0.0f',length(q)),...
		sprintf('Radius {mm}: %0.0f',R),...
		[sprintf('Variance: %0.1f',s(1)*100/sum(s)) '%']};
xlabel(str)

% create structure
%---------------------------------------------------------------------------
xY      = struct('name',	Rname,...
		 'Y',		Y,...
		 'y',		y,...
		 'XYZ',		L,...
		 'radius',	R);		
% save
%---------------------------------------------------------------------------
save([SPM.swd '/VOI_' Rname],'Y','xY')

%-Reset title
%-----------------------------------------------------------------------
spm('FigName',['SPM{',SPM.STAT,'}: Results']);
