function [Y,xY] = spm_regions(xSPM,SPM,hReg,xY)
% VOI time-series extraction of adjusted data (& local eigenimage analysis)
% FORMAT [Y xY] = spm_regions(xSPM,SPM,hReg,[xY]);
%
% xSPM   - structure containing specific SPM, distribution & filtering details
% SPM    - structure containing generic analysis details
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% Y      - first scaled eigenvariate of VOI {i.e. weighted mean}
% xY     - VOI structure
%       xY.xyz          - centre of VOI {mm}
%       xY.name         - name of VOI
%       xY.Ic           - contrast used to adjust data (0 - no adjustment)
%       xY.Sess         - session index
%       xY.def          - VOI definition
%       xY.spec         - VOI definition parameters
%       xY.str          - VOI description as a string
%       xY.XYZmm        - Co-ordinates of VOI voxels {mm}
%       xY.y            - [whitened and filtered] voxel-wise data
%       xY.u            - first eigenvariate {scaled - c.f. mean response}
%       xY.v            - first eigenimage
%       xY.s            - eigenvalues
%       xY.X0           - [whitened] confounds (including drift terms)
%
% Y and xY are also saved in VOI_*.mat in the SPM working directory
%
% (See spm_getSPM for details on the SPM & xSPM structures.)
%
%__________________________________________________________________________
%
% spm_regions extracts a representative time course from voxel data in
% terms of the first eigenvariate of the filtered and adjusted response in
% all suprathreshold voxels within a specified VOI centered on the current
% MIP cursor location.
%
% If temporal filtering has been specified, then the data will be filtered.
% Similarly for whitening. Adjustment is with respect to the null space of
% a selected contrast, or can be omitted.
%
% For a VOI of radius 0, the [adjusted] voxel time-series is returned, and
% scaled to have a 2-norm of 1. The actual [adjusted] voxel time series can
% be extracted from xY.y, and will be the same as the [adjusted] data 
% returned by the plotting routine (spm_graph.m) for the same contrast.
%__________________________________________________________________________
% Copyright (C) 1999-2011 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_regions.m 4185 2011-02-01 18:46:18Z guillaume $

if nargin < 4, xY = []; end

%-Get figure handles
%--------------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
if isempty(Finter), noGraph = 1; else noGraph = 0; end
header = get(Finter,'Name');
set(Finter,'Name','VOI time-series extraction');
if ~noGraph, Fgraph = spm_figure('GetWin','Graphics'); end

%-Find nearest voxel [Euclidean distance] in point list
%--------------------------------------------------------------------------
if isempty(xSPM.XYZmm)
    spm('alert!','No suprathreshold voxels!',mfilename,0);
    Y = []; xY = [];
    return
end
try
    xyz    = xY.xyz;
catch
    xyz    = spm_XYZreg('NearestXYZ',...
             spm_XYZreg('GetCoords',hReg),xSPM.XYZmm);
    xY.xyz = xyz;
end

% and update GUI location
%--------------------------------------------------------------------------
spm_XYZreg('SetCoords',xyz,hReg);


%-Get adjustment options and VOI name
%--------------------------------------------------------------------------
if ~noGraph
    if ~isempty(xY.xyz)
        posstr = sprintf('at [%3.0f %3.0f %3.0f]',xY.xyz);
    else
        posstr = '';
    end
    spm_input(posstr,1,'d','VOI time-series extraction');
end

if ~isfield(xY,'name')
    xY.name    = spm_input('name of region','!+1','s','VOI');
end

if ~isfield(xY,'Ic')
    q     = 0;
    Con   = {'<don''t adjust>'};
    for i = 1:length(SPM.xCon)
        if strcmp(SPM.xCon(i).STAT,'F')
            q(end + 1) = i;
            Con{end + 1} = SPM.xCon(i).name;
        end
    end
    i     = spm_input('adjust data for (select contrast)','!+1','m',Con);
    xY.Ic = q(i);
end

%-If fMRI data then ask user to select session
%--------------------------------------------------------------------------
if isfield(SPM,'Sess') && ~isfield(xY,'Sess')
    s       = length(SPM.Sess);
    if s > 1
        s   = spm_input('which session','!+1','n1',s,s);
    end
    xY.Sess = s;
end

%-Specify VOI
%--------------------------------------------------------------------------
xY.M = xSPM.M;
[xY, xY.XYZmm, Q] = spm_ROI(xY, xSPM.XYZmm);
try, xY = rmfield(xY,'M'); end
try, xY = rmfield(xY,'rej'); end

if isempty(xY.XYZmm)
    warning('Empty region.');
    Y = [];
    return;
end


%-Extract required data from results files
%==========================================================================
spm('Pointer','Watch')

%-Get raw data, whiten and filter 
%--------------------------------------------------------------------------
y        = spm_get_data(SPM.xY.VY,xSPM.XYZ(:,Q));
y        = spm_filter(SPM.xX.K,SPM.xX.W*y);


%-Computation
%==========================================================================

%-Remove null space of contrast
%--------------------------------------------------------------------------
if xY.Ic

    %-Parameter estimates: beta = xX.pKX*xX.K*y
    %----------------------------------------------------------------------
    beta  = spm_get_data(SPM.Vbeta,xSPM.XYZ(:,Q));

    %-subtract Y0 = XO*beta,  Y = Yc + Y0 + e
    %----------------------------------------------------------------------
    y     = y - spm_FcUtil('Y0',SPM.xCon(xY.Ic),SPM.xX.xKXs,beta);

end

%-Confounds
%--------------------------------------------------------------------------
xY.X0     = SPM.xX.xKXs.X(:,[SPM.xX.iB SPM.xX.iG]);

%-Extract session-specific rows from data and confounds
%--------------------------------------------------------------------------
try
    i     = SPM.Sess(xY.Sess).row;
    y     = y(i,:);
    xY.X0 = xY.X0(i,:);
end

% and add session-specific filter confounds
%--------------------------------------------------------------------------
try
    xY.X0 = [xY.X0 SPM.xX.K(xY.Sess).X0];
end
try
    xY.X0 = [xY.X0 SPM.xX.K(xY.Sess).KH]; % Compatibility check
end

%-Remove null space of X0
%--------------------------------------------------------------------------
xY.X0     = xY.X0(:,any(xY.X0));


%-Compute regional response in terms of first eigenvariate
%--------------------------------------------------------------------------
[m n]   = size(y);
if m > n
    [v s v] = svd(y'*y);
    s       = diag(s);
    v       = v(:,1);
    u       = y*v/sqrt(s(1));
else
    [u s u] = svd(y*y');
    s       = diag(s);
    u       = u(:,1);
    v       = y'*u/sqrt(s(1));
end
d       = sign(sum(v));
u       = u*d;
v       = v*d;
Y       = u*sqrt(s(1)/n);

%-Set in structure
%--------------------------------------------------------------------------
xY.y    = y;
xY.u    = Y;
xY.v    = v;
xY.s    = s;

%-Display VOI weighting and eigenvariate
%==========================================================================
if ~noGraph
    
    % show position
    %----------------------------------------------------------------------
    spm_results_ui('Clear',Fgraph);
    figure(Fgraph);
    subplot(2,2,3)
    spm_dcm_display(xY)

    % show dynamics
    %----------------------------------------------------------------------
    subplot(2,2,4)
    try
        plot(SPM.xY.RT*[1:length(xY.u)],Y)
        str = 'time (seconds}';
    catch
        plot(Y)
        str = 'scan';
    end
    title(['1st eigenvariate: ' xY.name],'FontSize',10)
    if strcmpi(xY.def,'mask')
        [p,n,e] = fileparts(xY.spec.fname);
        posstr  = sprintf('from mask %s', [n e]);
    else
        posstr  = sprintf('at [%3.0f %3.0f %3.0f]',xY.xyz);
    end
    str = { str;' ';...
        sprintf('%d voxels in VOI %s',length(Q),posstr);...
        sprintf('Variance: %0.2f%%',s(1)*100/sum(s))};
    xlabel(str)
    axis tight square
end

%-Save
%==========================================================================
str = ['VOI_' xY.name '.mat'];
if isfield(xY,'Sess') && isfield(SPM,'Sess')
    str = sprintf('VOI_%s_%i.mat',xY.name,xY.Sess);
end
if spm_check_version('matlab','7') >= 0
    save(fullfile(SPM.swd,str),'-V6','Y','xY')
else
    save(fullfile(SPM.swd,str),'Y','xY')
end

fprintf('   VOI saved as %s\n',spm_str_manip(fullfile(SPM.swd,str),'k55'));

%-Reset title
%--------------------------------------------------------------------------
set(Finter,'Name',header);
spm('Pointer','Arrow')
