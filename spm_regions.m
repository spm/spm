function [Y,xY] = spm_regions(SPM,VOL,xX,xCon,xSDM,hReg,xY)
% VOI time-series extraction of adjusted data (local eigenimage analysis)
% FORMAT [Y xY] = spm_regions(SPM,VOL,xX,xCon,xSDM,hReg,xY);
%
% SPM    - structure containing SPM, distribution & filtering detals
% VOL    - structure containing details of volume analysed
% xX     - Design Matrix structure
% xSDM   - structure containing contents of SPM.mat file
% xCon   - Contrast definitions structure (see spm_FcUtil.m for details)
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% Y      - first scaled eigenvariate of VOI {i.e. weighted mean}
% xY     - VOI structure
%       xY.xyz          - centre of VOI {mm}
%       xY.name         - name of VOI
%       xY.Ic           - contrast used to adjust data (0 - no adjustment)
%       xY.filter       - filtering (yes|no)
%       xY.Sess         - session indices
%       xY.def          - VOI definition
%       xY.spec         - VOI definition parameters
%       xY.XYZmm        - Co-ordinates of VOI voxels {mm}
%       xY.y            - voxel-wise data
%       xY.u            - first eigenvariate {scaled - c.f. mean response}
%       xY.v            - first eigenimage
%       xY.s            - eigenvalues
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
header = get(Finter,'Name');
set(Finter,'Name','VOI time-series extraction')
if ~exist('xY')
	xY = {};
end

%-Find nearest voxel [Euclidean distance] in point list in Y.mad
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

if ~isfield(xY,'xyz')
	[xyz,i] = spm_XYZreg('NearestXYZ',...
		      spm_XYZreg('GetCoords',hReg),SPM.XYZmm(:,find(SPM.QQ)));
	xY.xyz  = xyz;
else
	xyz     = xY.xyz;
end

% and update GUI location
%-----------------------------------------------------------------------
spm_XYZreg('SetCoords',xyz,hReg);


%-Get adjustment options and VOI name
%-----------------------------------------------------------------------
spm_input(sprintf('at [%3.0f %3.0f %3.0f]',xY.xyz),1,'d',...
	'VOI time-series extraction')

if ~isfield(xY,'name')
	xY.name = spm_input('name of region','!+1','s','VOI');
end

if ~isfield(xY,'Ic')
	q     = 0;
	Con   = {'<don''t adjust>'};
	for i = 1:length(xCon)
		if strcmp(xCon(i).STAT,'F')
			  q(end + 1) = i;
			Con{end + 1} = xCon(i).name;
		end
	end
	i     = spm_input('adjust data for (select contrast)','!+1','m',Con);
	xY.Ic = q(i);
end

%-if fMRI data get sessions and filtering options
%-----------------------------------------------------------------------
if isfield(xSDM,'Sess')

	if ~isfield(xY,'filter')
		xY.filter = spm_input('filter','!+1','apply|none|high');
	end

	if ~isfield(xY,'Sess')
		s         = length(xSDM.Sess);
		s         = spm_input('which session[s]','+1','n',1,[1 Inf],s);
		xY.Sess   = s;
	end
end


%-Specify VOI
%-----------------------------------------------------------------------
if ~isfield(xY,'def')
	xY.def    = spm_input('VOI definition...','!+1','b',...
			{'sphere','box','cluster'});
end
Q       = ones(1,size(SPM.XYZmm,2));


switch xY.def

	case 'sphere'
	%---------------------------------------------------------------
	if ~isfield(xY,'spec')
		xY.spec = spm_input('VOI radius (mm)','!+1','r',0,1,[0,Inf]);
	end
	d     = [SPM.XYZmm(1,:)-xyz(1);
		 SPM.XYZmm(2,:)-xyz(2);
		 SPM.XYZmm(3,:)-xyz(3)];
	Q     = find(sum(d.^2) <= xY.spec^2);

	case 'box'
	%---------------------------------------------------------------
	if ~isfield(xY,'spec')
		xY.spec = spm_input('box dimensions [x y z] {mm}',...
			'!+1','r','0 0 0',3);
	end
	Q     = find(all(abs(SPM.XYZmm - xyz*Q) <= xY.spec(:)*Q/2));

	case 'cluster'
	%---------------------------------------------------------------
	[x i] = spm_XYZreg('NearestXYZ',xyz,SPM.XYZmm);
	A     = spm_clusters(SPM.XYZ);
	Q     = find(A == A(i));
end

% voxels defined
%-----------------------------------------------------------------------
spm('Pointer','Watch')
q      = find(SPM.QQ(Q));
if any(SPM.QQ(Q)==0)
    spm('alert"',{...
        sprintf('Incomplete data for all %d suprathreshold',length(Q)),...
	sprintf('voxels.  Proceeding using the %d voxels saved.',length(q)),...
	},mfilename,sqrt(-1));
end



%-Extract required data from results files
%=======================================================================

%-Get (approximate) raw data y from Y.mad file
%-NB: Data in Y.mad file is compressed, and therefore not fully accurate
%-----------------------------------------------------------------------
y        = spm_extract(fullfile(SPM.swd,'Y.mad'),SPM.QQ(Q(q)));
rcp      = VOL.iM(1:3,:)*[SPM.XYZmm(:,Q(q));ones(size(q))];
xY.XYZmm = SPM.XYZmm(:,Q(q));



%-Computation
%=======================================================================

% remove null space of contrast (prior to filering)
%-----------------------------------------------------------------------
if xY.Ic

	%-Parameter estimates: beta = xX.pKX*xX.K*y
	%---------------------------------------------------------------
	nBeta   = length(xSDM.Vbeta);
	beta    = ones(nBeta,size(q,2));
	for   i = 1:nBeta
		beta(i,:) = ...
		spm_sample_vol(xSDM.Vbeta(i),rcp(1,:),rcp(2,:),rcp(3,:),0);
	end

	Fc      = spm_FcUtil('Set','Fc','F','iX0',xCon(xY.Ic).iX0,xX.X);
	y       = y - spm_FcUtil('Y0',Fc,xX.X,beta);
end

% fMRI-specific operations
%-----------------------------------------------------------------------
if isfield(xY,'Sess')

	% filter
	%---------------------------------------------------------------
	y       = spm_filter(xY.filter,xX.K,y);

	% extract sessions
	%---------------------------------------------------------------
	j       = [];
	for   i = 1:length(xY.Sess)
		j = [j xSDM.Sess{xY.Sess(i)}.row];
	end
	y       =  y(j,:);
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
Y       = u*sqrt(s(1)/n);

% set in structure
%-----------------------------------------------------------------------
xY.y    = y;
xY.u    = Y;
xY.v    = v;
xY.s    = s;

%-Display VOI weighting and eigenvariate
%-----------------------------------------------------------------------
if nargin < 7


	% show position
	%---------------------------------------------------------------
	spm_results_ui('Clear',Fgraph);
	figure(Fgraph);
	subplot(2,2,3)
	spm_dcm_display(xY,[],[],[],[[1 0 0];[0 1 0]]',64)


	% show dynamics
	%---------------------------------------------------------------
	subplot(2,2,4)
	if isfield(xX,'RT')
		plot(xX.RT*[1:length(xY.u)],Y)
		str = 'time (seconds}';
	else
		plot(Y)
		str = 'scan';
	end
	title(['1st eigenvariate: ' xY.name],'FontSize',10)
	str = {	str;' ';sprintf(...
		'%d voxels in VOI at [%3.0f %3.0f %3.0f]',...
		length(q),xyz);sprintf('Variance: %0.2f%%',s(1)*100/sum(s))};
	xlabel(str)
	axis tight square
end


% save
%-----------------------------------------------------------------------
str     = ['VOI_' xY.name];
if isfield(xY,'Sess')
	if length(xY.Sess) == 1
		str = sprintf('VOI_%s_%i',xY.name,xY.Sess);
	end
end
save(fullfile(SPM.swd,str),'Y','xY')

%-Reset title
%-----------------------------------------------------------------------
spm('FigName',header);
spm('Pointer','Arrow')
