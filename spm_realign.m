function P = spm_realign(P,flags)
% Estimation of within modality rigid body movement parameters
% FORMAT P = spm_realign(P,flags)
%
% P     - char array of filenames
%         All operations are performed relative to the first image.
%         ie. Coregistration is to the first image, and resampling
%         of images is into the space of the first image.
%         For multiple sessions, P should be a cell array, where each
%         cell should be a matrix of filenames.
%
% flags - a structure containing various options.  The fields are:
%         quality  - Quality versus speed trade-off.  Highest quality (1)
%                    gives most precise results, whereas lower qualities
%                    gives faster realignment.
%                    The idea is that some voxels contribute little to
%                    the estimation of the realignment parameters.
%                    This parameter is involved in selecting the number
%                    of voxels that are used.
%
%         fwhm     - The FWHM of the Gaussian smoothing kernel (mm) applied
%                    to the images before estimating the realignment
%                    parameters.
%
%         sep      - the default separation (mm) to sample the images.
%
%         rtm      - Register to mean.  If field exists then a two pass
%                    procedure is to be used in order to register the
%                    images to the mean of the images after the first
%                    realignment.
%
%         wrap     - Directions in the volume whose values should wrap
%                    around in. For example, in MRI scans, the images wrap
%                    around in the phase encode direction, so (e.g.) the
%                    subject's nose may poke into the back of the subject's
%                    head.
%
%         PW       -  a filename of a weighting image (reciprocal of
%                    standard deviation).  If field does not exist, then
%                    no weighting is done.
%
%         interp   - B-spline degree used for interpolation
%
%         graphics - display coregistration outputs
%                    default: ~spm('CmdLine')
%
%__________________________________________________________________________
%
% If no output argument, then an updated voxel to world matrix is written
% to the headers of the images (a .mat file is created for 4D images).
% The details of the transformation are displayed in the results window as
% plots of translation and rotation.
% A set of realignment parameters are saved for each session, named:
% rp_*.txt.
%__________________________________________________________________________
%
% Voxel to world mapping:
%
% These are simply 4x4 affine transformation matrices represented in the
% NIFTI headers (see http://nifti.nimh.nih.gov/nifti-1 ).
% These are normally modified by the `realignment' and `coregistration'
% modules.  What these matrices represent is a mapping from the voxel
% coordinates (x0,y0,z0) (where the first voxel is at coordinate (1,1,1)),
% to coordinates in millimeters (x1,y1,z1).
%  
% x1 = M(1,1)*x0 + M(1,2)*y0 + M(1,3)*z0 + M(1,4)
% y1 = M(2,1)*x0 + M(2,2)*y0 + M(2,3)*z0 + M(2,4)
% z1 = M(3,1)*x0 + M(3,2)*y0 + M(3,3)*z0 + M(3,4)
%
% Assuming that image1 has a transformation matrix M1, and image2 has a
% transformation matrix M2, the mapping from image1 to image2 is: M2\M1
% (ie. from the coordinate system of image1 into millimeters, followed
% by a mapping from millimeters into the space of image2).
%
% These matrices allow several realignment or coregistration steps to be
% combined into a single operation (without the necessity of resampling the
% images several times).
%__________________________________________________________________________
%
% Reference:
%
% Friston KJ, Ashburner J, Frith CD, Poline J-B, Heather JD & Frackowiak
% RSJ (1995) Spatial registration and normalization of images Hum. Brain
% Map. 2:165-189
%__________________________________________________________________________
% Copyright (C) 1994-2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_realign.m 7141 2017-07-26 09:05:05Z guillaume $


SVNid = '$Rev: 7141 $';
 
%-Say hello
%--------------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,SVNid);

%-Parameters & Arguments
%==========================================================================

%-Flags
%--------------------------------------------------------------------------
def_flags          = spm_get_defaults('realign.estimate');
def_flags.PW       = '';
def_flags.graphics = ~spm('CmdLine');
def_flags.lkp      = 1:6;
if nargin < 2
    flags = def_flags;
else
    fnms = fieldnames(def_flags);
    for i=1:length(fnms)
        if ~isfield(flags,fnms{i})
            flags.(fnms{i}) = def_flags.(fnms{i});
        end
    end
end

%-Images
%--------------------------------------------------------------------------
if ~iscell(P), P = {P}; end
for i=1:numel(P), if ischar(P{i}), P{i} = spm_vol(P{i}); end; end
P(cellfun(@isempty,P)) = [];
if ~isempty(flags.PW) && ischar(flags.PW), flags.PW = spm_vol(flags.PW); end

if isempty(P), warning('Nothing to do.'); end

%-Perform realignment
%==========================================================================
if numel(P)==1
    P{1} = realign_series(P{1}, flags);
elseif numel(P) > 1
    Ptmp = P{1}(1);
    for s=2:numel(P)
        Ptmp = [Ptmp ; P{s}(1)];
    end
    Ptmp = realign_series(Ptmp,flags);
    for s=1:numel(P)
        M  = Ptmp(s).mat / P{s}(1).mat;
        for i=1:numel(P{s})
            P{s}(i).mat = M * P{s}(i).mat;
        end
    end

    for s=1:numel(P)
        P{s} = realign_series(P{s},flags);
    end
end

%-Save results
%==========================================================================

if ~nargout
    for s=1:numel(P)
        %-Save parameters as rp_*.txt files
        %------------------------------------------------------------------
        save_parameters(P{s});
        
        %-Update voxel to world mapping in images header
        %------------------------------------------------------------------
        for i=1:numel(P{s})
            spm_get_space([P{s}(i).fname ',' num2str(P{s}(i).n)], P{s}(i).mat);
        end
    end
end

%-Display translation and rotation parameters
%--------------------------------------------------------------------------
if flags.graphics, plot_parameters(P); end


if numel(P)==1, P = P{1}; end

fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#

%==========================================================================
% function P = realign_series(P,flags)
%==========================================================================
function P = realign_series(P,flags)
% Realign a time series of 3D images to the first of the series
% FORMAT P = realign_series(P,flags)
% P  - a vector of volumes (see spm_vol)
%--------------------------------------------------------------------------
% P(i).mat is modified to reflect the modified position of the image i.
% The scaling (and offset) parameters are also set to contain the
% optimum scaling required to match the images.
%__________________________________________________________________________

if numel(P)<2, return; end

skip = sqrt(sum(P(1).mat(1:3,1:3).^2)).^(-1)*flags.sep;
d    = P(1).dim(1:3);                                                                                                                        
lkp  = flags.lkp;
st   = rand('state'); % st = rng;
rand('state',0); % rng(0,'v5uniform'); % rng('defaults');
if d(3) < 3
    lkp  = [1 2 6];
    [x1,x2,x3] = ndgrid(1:skip(1):d(1)-.5, 1:skip(2):d(2)-.5, 1:skip(3):d(3));
    x1   = x1 + rand(size(x1))*0.5;
    x2   = x2 + rand(size(x2))*0.5;
else
    [x1,x2,x3] = ndgrid(1:skip(1):d(1)-.5, 1:skip(2):d(2)-.5, 1:skip(3):d(3)-.5);
    x1   = x1 + rand(size(x1))*0.5;
    x2   = x2 + rand(size(x2))*0.5;
    x3   = x3 + rand(size(x3))*0.5; 
end
rand('state',st); % rng(st);

x1 = x1(:);
x2 = x2(:);
x3 = x3(:);

%-Possibly mask an area of the sample volume.
%--------------------------------------------------------------------------
if ~isempty(flags.PW)
    [y1,y2,y3] = coords([0 0 0  0 0 0],P(1).mat,flags.PW.mat,x1,x2,x3);
    wt  = spm_sample_vol(flags.PW,y1,y2,y3,1);
    msk = find(wt>0.01);
    x1  = x1(msk);
    x2  = x2(msk);
    x3  = x3(msk);
    wt  = wt(msk);
else
    wt  = [];
end

%-Compute rate of change of chi2 w.r.t changes in parameters (matrix A)
%--------------------------------------------------------------------------
V   = smooth_vol(P(1),flags.interp,flags.wrap,flags.fwhm);
deg = [flags.interp*[1 1 1]' flags.wrap(:)];

[G,dG1,dG2,dG3] = spm_bsplins(V,x1,x2,x3,deg);
clear V
A0  = make_A(P(1).mat,x1,x2,x3,dG1,dG2,dG3,wt,lkp);

b   = G;
if ~isempty(wt), b = b.*wt; end

%-Remove voxels that contribute very little to the final estimate
%--------------------------------------------------------------------------
if numel(P) > 2
    % Simulated annealing or something similar could be used to
    % eliminate a better choice of voxels - but this way will do for
    % now. It basically involves removing the voxels that contribute
    % least to the determinant of the inverse covariance matrix.

    spm_plot_convergence('Init','Eliminating Unimportant Voxels',...
              'Relative quality','Iteration');
    Alpha = [A0 b];
    Alpha = Alpha'*Alpha;
    det0  = det(Alpha);
    det1  = det0;
    spm_plot_convergence('Set',det1/det0);
    while det1/det0 > flags.quality
        dets = zeros(size(A0,1),1);
        for i=1:size(A0,1)
            tmp     = [A0(i,:) b(i)];
            dets(i) = det(Alpha - tmp'*tmp);
        end
        clear tmp
        [junk,msk] = sort(det1-dets);
        msk        = msk(1:round(length(dets)/10));
         A0(msk,:) = [];   b(msk,:) = [];   G(msk,:) = [];
         x1(msk,:) = [];  x2(msk,:) = [];  x3(msk,:) = [];
        dG1(msk,:) = []; dG2(msk,:) = []; dG3(msk,:) = [];
        if ~isempty(wt),  wt(msk,:) = []; end
        Alpha = [A0 b];
        Alpha = Alpha'*Alpha;
        det1  = det(Alpha);
        spm_plot_convergence('Set',single(det1/det0));
    end
    spm_plot_convergence('Clear');
end

if flags.rtm
    count = ones(size(b));
    ave   = G;
    grad1 = dG1;
    grad2 = dG2;
    grad3 = dG3;
end

%-Loop over images
%--------------------------------------------------------------------------

spm_progress_bar('Init',numel(P)-1,'Registering Images');
for i=2:numel(P)
    V  = smooth_vol(P(i),flags.interp,flags.wrap,flags.fwhm);
    d  = [size(V) 1 1];
    d  = d(1:3);
    ss = Inf;
    countdown = -1;
    for iter=1:64
        [y1,y2,y3] = coords([0 0 0  0 0 0],P(1).mat,P(i).mat,x1,x2,x3);
        msk        = find((y1>=1 & y1<=d(1) & y2>=1 & y2<=d(2) & y3>=1 & y3<=d(3)));
        if length(msk)<32, error_message(P(i)); end

        F          = spm_bsplins(V, y1(msk),y2(msk),y3(msk),deg);
        if ~isempty(wt), F = F.*wt(msk); end

        A          = A0(msk,:);
        b1         = b(msk);
        sc         = sum(b1)/sum(F);
        b1         = b1-F*sc;
        soln       = (A'*A)\(A'*b1);

        p          = [0 0 0  0 0 0  1 1 1  0 0 0];
        p(lkp)     = p(lkp) + soln';
        P(i).mat   = spm_matrix(p) \ P(i).mat;

        pss        = ss;
        ss         = sum(b1.^2)/length(b1);
        if (pss-ss)/pss < 1e-8 && countdown == -1 % Stopped converging.
            countdown = 2;
        end
        if countdown ~= -1
            if countdown==0, break; end
            countdown = countdown -1;
        end
    end
    if flags.rtm
        % Generate mean and derivatives of mean
        tiny = 5e-2; % From spm_vol_utils.c
        msk        = find((y1>=(1-tiny) & y1<=(d(1)+tiny) &...
                           y2>=(1-tiny) & y2<=(d(2)+tiny) &...
                           y3>=(1-tiny) & y3<=(d(3)+tiny)));
        count(msk) = count(msk) + 1;
        [G,dG1,dG2,dG3] = spm_bsplins(V,y1(msk),y2(msk),y3(msk),deg);
        ave(msk)   = ave(msk)   +   G*sc;
        grad1(msk) = grad1(msk) + dG1*sc;
        grad2(msk) = grad2(msk) + dG2*sc;
        grad3(msk) = grad3(msk) + dG3*sc;
    end
    spm_progress_bar('Set',i-1);
end
spm_progress_bar('Clear');

if ~flags.rtm, return; end


%-Register to mean
%--------------------------------------------------------------------------
M  = P(1).mat;
A0 = make_A(M,x1,x2,x3,grad1./count,grad2./count,grad3./count,wt,lkp);
if ~isempty(wt), b = (ave./count).*wt;
else b = (ave./count); end

clear ave grad1 grad2 grad3

% Loop over images
%--------------------------------------------------------------------------
spm_progress_bar('Init',length(P),'Registering Images to Mean');
for i=1:numel(P)
    V  = smooth_vol(P(i),flags.interp,flags.wrap,flags.fwhm);
    d  = [size(V) 1 1 1];
    ss = Inf;
    countdown = -1;
    for iter=1:64
        [y1,y2,y3] = coords([0 0 0  0 0 0],M,P(i).mat,x1,x2,x3);
        msk        = find((y1>=1 & y1<=d(1) & y2>=1 & y2<=d(2) & y3>=1 & y3<=d(3)));
        if length(msk)<32, error_message(P(i)); end

        F          = spm_bsplins(V, y1(msk),y2(msk),y3(msk),deg);
        if ~isempty(wt), F = F.*wt(msk); end

        A          = A0(msk,:);
        b1         = b(msk);
        sc         = sum(b1)/sum(F);
        b1         = b1-F*sc;
        soln       = (A'*A)\(A'*b1);

        p          = [0 0 0  0 0 0  1 1 1  0 0 0];
        p(lkp)     = p(lkp) + soln';
        P(i).mat   = spm_matrix(p) \ P(i).mat;

        pss        = ss;
        ss         = sum(b1.^2)/length(b1);
        if (pss-ss)/pss < 1e-8 && countdown == -1 % Stopped converging.
            % Do three final iterations to finish off with
            countdown = 2;
        end
        if countdown ~= -1
            if countdown==0, break; end
            countdown = countdown -1;
        end
    end
    spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');

% Aligning everything to the first image
%--------------------------------------------------------------------------
M = M / P(1).mat;
for i=1:numel(P)
    P(i).mat = M * P(i).mat;
end


%==========================================================================
% function [y1,y2,y3]=coords(p,M1,M2,x1,x2,x3)
%==========================================================================
function [y1,y2,y3]=coords(p,M1,M2,x1,x2,x3)
% Rigid body transformation of a set of coordinates
M  = inv(M2) * inv(spm_matrix(p)) * M1;
y1 = M(1,1)*x1 + M(1,2)*x2 + M(1,3)*x3 + M(1,4);
y2 = M(2,1)*x1 + M(2,2)*x2 + M(2,3)*x3 + M(2,4);
y3 = M(3,1)*x1 + M(3,2)*x2 + M(3,3)*x3 + M(3,4);


%==========================================================================
% function V = smooth_vol(P,hld,wrp,fwhm)
%==========================================================================
function V = smooth_vol(P,hld,wrp,fwhm)
% Convolve the volume in memory
vx = sqrt(sum(P.mat(1:3,1:3).^2));
s  = vx.^(-1)*(fwhm/sqrt(8*log(2)));
x  = round(6*s(1)); x = -x:x;
y  = round(6*s(2)); y = -y:y;
z  = round(6*s(3)); z = -z:z;

x = spm_smoothkern(fwhm/vx(1),x,1);
y = spm_smoothkern(fwhm/vx(2),y,1);
z = spm_smoothkern(fwhm/vx(3),z,1);

%x  = exp(-(x).^2/(2*(s(1)).^2));
%y  = exp(-(y).^2/(2*(s(2)).^2));
%z  = exp(-(z).^2/(2*(s(3)).^2));
%x  = x/sum(x);
%y  = y/sum(y);
%z  = z/sum(z);

i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;
d  = [hld*[1 1 1]' wrp(:)];
V  = spm_bsplinc(P,d);
spm_conv_vol(V,V,x,y,z,-[i j k]);


%==========================================================================
% function A = make_A(M,x1,x2,x3,dG1,dG2,dG3,wt,lkp)
%==========================================================================
function A = make_A(M,x1,x2,x3,dG1,dG2,dG3,wt,lkp)
% Matrix of rate of change of weighted difference w.r.t. parameter changes
p0 = [0 0 0  0 0 0  1 1 1  0 0 0];
A  = zeros(numel(x1),length(lkp));
for i=1:length(lkp)
    pt         = p0;
    pt(lkp(i)) = pt(i)+1e-6;
    [y1,y2,y3] = coords(pt,M,M,x1,x2,x3);
    tmp        = sum([y1-x1 y2-x2 y3-x3].*[dG1 dG2 dG3],2)/(-1e-6);
    if ~isempty(wt), A(:,i) = tmp.*wt;
    else A(:,i) = tmp; end
end


%==========================================================================
% function error_message(P)
%==========================================================================
function error_message(P)
str = {'There is not enough overlap in the images to obtain a solution.',...
       ' ',...
       'Offending image:',...
       P.fname,...
       ' ',...
       'Please check that your header information is OK.',...
       'The Check Reg utility will show you the initial',...
       'alignment between the images, which must be',...
       'within about 4cm and about 15 degrees in order',...
       'for SPM to find the optimal solution.'};
spm('alert*',str,mfilename,sqrt(-1));
error('Insufficient image overlap.');


%==========================================================================
% function plot_parameters(P)
%==========================================================================
function plot_parameters(P)
fg = spm_figure('FindWin','Graphics');
if isempty(fg), return; end

P = cat(1,P{:});
if length(P)<2, return; end
Params = zeros(numel(P),12);
for i=1:numel(P)
    Params(i,:) = spm_imatrix(P(i).mat/P(1).mat);
end

%-Display results: translation and rotation over time series
%--------------------------------------------------------------------------
spm_figure('Clear','Graphics');
ax = axes('Position',[0.1 0.65 0.8 0.2],'Parent',fg,'Visible','off');
set(get(ax,'Title'),'String','Image realignment',...
    'FontSize',16,'FontWeight','Bold','Visible','on');
x     =  0.1;
y     =  0.9;
for i = 1:min([numel(P) 12])
    text(x,y,[sprintf('%-4.0f',i) P(i).fname ',' num2str(P(i).n(1))],...
        'FontSize',10,'Interpreter','none','Parent',ax);
    y = y - 0.08;
end
if numel(P) > 12
    text(x,y,'................ etc','FontSize',10,'Parent',ax); end

ax = axes('Position',[0.1 0.35 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on',...
    'NextPlot','replacechildren','ColorOrder',[0 0 1;0 0.5 0;1 0 0]);
plot(Params(:,1:3),'Parent',ax)
s  = {'x translation','y translation','z translation'};
%text([2 2 2], Params(2, 1:3), s, 'Fontsize',10,'Parent',ax)
legend(ax, s, 'Location','Best')
set(get(ax,'Title'),'String','translation','FontSize',16,'FontWeight','Bold');
set(get(ax,'Xlabel'),'String','image');
set(get(ax,'Ylabel'),'String','mm');


ax = axes('Position',[0.1 0.05 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on',...
    'NextPlot','replacechildren','ColorOrder',[0 0 1;0 0.5 0;1 0 0]);
plot(Params(:,4:6)*180/pi,'Parent',ax)
s  = {'pitch','roll','yaw'};
%text([2 2 2], Params(2, 4:6)*180/pi, s, 'Fontsize',10,'Parent',ax)
legend(ax, s, 'Location','Best')
set(get(ax,'Title'),'String','rotation','FontSize',16,'FontWeight','Bold');
set(get(ax,'Xlabel'),'String','image');
set(get(ax,'Ylabel'),'String','degrees');

%-Print realigment parameters
%--------------------------------------------------------------------------
spm_print;


%==========================================================================
% function save_parameters(V)
%==========================================================================
function save_parameters(V)
fname = spm_file(V(1).fname, 'prefix','rp_', 'ext','.txt');
n = length(V);
Q = zeros(n,6);
for j=1:n
    qq     = spm_imatrix(V(j).mat/V(1).mat);
    Q(j,:) = qq(1:6);
end
save(fname,'Q','-ascii');
