function x = spm_coreg(varargin)
% Between modality coregistration using information theory
% FORMAT x = spm_coreg(VG,VF,params)
% VG - handle for first image (see spm_vol).
% VF - handle for second image.
% x - the parameters describing the rigid body rotation.
%     such that a mapping from voxels in G to voxels in F
%     is attained by:  VF.mat\spm_matrix(x(:)')*VG.mat
% flags - a structure containing the following elements:
%          sep      - optimisation sampling steps (mm)
%                     default: [4 2]
%          params   - starting estimates (6 elements)
%                     default: [0 0 0  0 0 0]
%          cost_fun - cost function string:
%                      'mi'  - Mutual Information
%                      'nmi' - Normalised Mutual Information
%                      'ecc' - Entropy Correlation Coefficient
%                      'ncc' - Normalised Cross Correlation
%                      default: 'nmi'
%          tol      - tolerences for accuracy of each param
%                     default: [0.02 0.02 0.02 0.001 0.001 0.001]
%          fwhm     - smoothing to apply to 256x256 joint histogram
%                     default: [7 7]
%
% The registration method used here is based on the work described in:
% A Collignon, F Maes, D Delaere, D Vandermeulen, P Suetens & G Marchal
% (1995) "Automated Multi-modality Image Registration Based On
% Information Theory". In the proceedings of Information Processing in
% Medical Imaging (1995).  Y. Bizais et al. (eds.).  Kluwer Academic
% Publishers.
%
% The original interpolation method described in this paper has been
% changed in order to give a smoother cost function.  The images are
% also smoothed slightly, as is the histogram.  This is all in order to
% make the cost function as smooth as possible, to give faster
% convergence and less chance of local minima.
%
% References
% ==========
% Mutual Information
% ------------------
% Collignon, Maes, Delaere, Vandermeulen, Suetens & Marchal (1995).
% "Automated multi-modality image registration based on information theory".
% In Bizais, Barillot & Di Paola, editors, Proc. Information Processing
% in Medical Imaging, pages 263--274, Dordrecht, The Netherlands, 1995.
% Kluwer Academic Publishers.
%
% Wells III, Viola, Atsumi, Nakajima & Kikinis (1996).
% "Multi-modal volume registration by maximisation of mutual information".
% Medical Image Analysis, 1(1):35-51, 1996. 
%
% Entropy Correlation Coefficient
% -------------------------------
% F Maes, A Collignon, D Vandermeulen, G Marchal & P Suetens (1997).
% "Multimodality image registration by maximisation of mutual
% information". IEEE Transactions on Medical Imaging 16(2):187-198
%
% Normalised Mutual Information
% -----------------------------
% Studholme,  Hill & Hawkes (1998).
% "A normalized entropy measure of 3-D medical image alignment".
% in Proc. Medical Imaging 1998, vol. 3338, San Diego, CA, pp. 132-143.             
%
% Optimisation
% ------------
% Press, Teukolsky, Vetterling & Flannery (1992).
% "Numerical Recipes in C (Second Edition)".
% Published by Cambridge.
%
% At the end, the voxel-to-voxel affine transformation matrix is
% displayed, along with the histograms for the images in the original
% orientations, and the final orientations.  The registered images are
% displayed at the bottom.
%_______________________________________________________________________
% %W% John Ashburner %E%

if nargin>=4,
	x = optfun(varargin{:});
	return;
end;

def_flags = struct('sep',[4 2],'params',[0 0 0  0 0 0], 'cost_fun','nmi','fwhm',[7 7],...
	'tol',[0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001],'graphics',1);
if nargin < 3,
	flags = def_flags;
else,
	flags = varargin{3};
	fnms  = fieldnames(def_flags);
	for i=1:length(fnms),
		if ~isfield(flags,fnms{i}), flags = setfield(flags,fnms{i},getfield(def_flags,fnms{i})); end;
	end;
end;
%disp(flags)

if nargin < 1,
	VG = spm_vol(spm_get(1,'IMAGE','Select reference image'));
else,
	VG = varargin{1};
	if ischar(VG), VG = spm_vol(VG); end;
end;
if nargin < 2,
	VF = spm_vol(spm_get(1,'IMAGE','Select moved image'));
else,
	VF = varargin{2};
	if ischar(VF), VF = spm_vol(VF); end;
end;


% Voxel sizes (mm)
vxg   = sqrt(sum(VG.mat(1:3,1:3).^2));
vxf   = sqrt(sum(VF.mat(1:3,1:3).^2));

% Smoothnesses
fwhmg = sqrt(max([1 1 1]*flags.sep(end)^2 - vxg.^2, [0 0 0]))./vxg;
fwhmf = sqrt(max([1 1 1]*flags.sep(end)^2 - vxf.^2, [0 0 0]))./vxf;

if ~isfield(VG, 'uint8'),
	VG.uint8 = loaduint8(VG);
	VG       = smooth_uint8(VG,fwhmg); % Note side effects
end;
if ~isfield(VF, 'uint8'),
	VF.uint8 = loaduint8(VF);
	VF       = smooth_uint8(VF,fwhmf); % Note side effects
end;


sc = flags.tol(:)'; % Required accuracy
sc = sc(1:length(flags.params));
xi = diag(sc*20);
x  = flags.params(:);
for samp=flags.sep(:)',
	[x,fval] = spm_powell(x(:), xi,sc,mfilename,VG,VF,samp,flags.cost_fun,flags.fwhm);
	x        = x(:)';
end;
if flags.graphics,
	display_results(VG,VF,x,flags);
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function o = optfun(x,VG,VF,s,cf,fwhm)
% The function that is minimised.
if nargin<6, fwhm = [7 7];   end;
if nargin<5, cf   = 'mi';    end;
if nargin<4, s    = [1 1 1]; end;

% Voxel sizes
vxg   = sqrt(sum(VG.mat(1:3,1:3).^2));sg = s./vxg;

% Create the joint histogram
H =   spm_hist2(VG.uint8,VF.uint8, VF.mat\spm_matrix(x(:)')*VG.mat ,sg);

% Smooth the histogram
sm = max(fwhm/sqrt(8*log(2)),0.001); % FWHM -> Gaussian param
t  = max(round(3*sm(1)),0); krn1 = exp(-([-t:t].^2)/sm(1)^2) ; krn1 = krn1/sum(krn1) ; H = conv2(H,krn1);
t  = max(round(3*sm(2)),0); krn2 = exp(-([-t:t].^2)/sm(2)^2)'; krn2 = krn2/sum(krn2)'; H = conv2(H,krn2);
%d = 32;
%H = sum(reshape(H,[256/d d 256]),1);
%H = reshape(sum(reshape(H,[d 256/d d]),2),[d d]);

% Compute cost function from histogram
H  = H+eps;
sh = sum(H(:));
H  = H/sh;
s1 = sum(H,1);
s2 = sum(H,2);

switch lower(cf)
	case 'mi',
		% Mutual Information:
		H   = H.*log2(H./(s2*s1));
		mi  = sum(H(:));
		o   = -mi;
	case 'ecc',
		% Entropy Correlation Coefficient of:
		% Maes, Collignon, Vandermeulen, Marchal & Suetens (1997).
		% "Multimodality image registration by maximisation of mutual
		% information". IEEE Transactions on Medical Imaging 16(2):187-198
		H   = H.*log2(H./(s2*s1));
		mi  = sum(H(:));
		ecc = -2*mi/(sum(s1.*log2(s1))+sum(s2.*log2(s2)));
		o   = -ecc;
	case 'nmi',
		% Normalised Mutual Information of:
		% Studholme,  Hill & Hawkes (1998).
		% "A normalized entropy measure of 3-D medical image alignment".
		% in Proc. Medical Imaging 1998, vol. 3338, San Diego, CA, pp. 132-143.
		nmi = (sum(s1.*log2(s1))+sum(s2.*log2(s2)))/sum(sum(H.*log2(H)));
		o   = -nmi;
	case 'ncc',
		% Normalised Cross Correlation
		i     = 1:size(H,1);
		j     = 1:size(H,2);
		m1    = sum(s2.*i');
		m2    = sum(s1.*j);
		sig1  = sqrt(sum(s2.*(i'-m1).^2));
		sig2  = sqrt(sum(s1.*(j -m2).^2));
		[i,j] = ndgrid(i-m1,j-m2);
		ncc   = sum(sum(H.*i.*j))/(sig1*sig2);
		o     = -ncc;
	otherwise,
		error('Invalid cost function specified');
end;

return;
%_______________________________________________________________________

%_______________________________________________________________________
function udat = loaduint8(V)
% Load data from file indicated by V into an array of unsigned bytes.
if size(V.pinfo,2)==1 & V.pinfo(1) == 2,
	mx = 255*V.pinfo(1) + V.pinfo(2);
	mn = V.pinfo(2);
else,
	spm_progress_bar('Init',V.dim(3),...
		['Computing max/min of ' spm_str_manip(V.fname,'t')],...
		'Planes complete');
	mx = -Inf; mn =  Inf;
	for p=1:V.dim(3),
		img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
		mx  = max([max(img(:))+paccuracy(V,p) mx]);
		mn  = min([min(img(:)) mn]);
		spm_progress_bar('Set',p);
	end;
end;
spm_progress_bar('Init',V.dim(3),...
	['Loading ' spm_str_manip(V.fname,'t')],...
	'Planes loaded');

udat = uint8(0);
udat(V.dim(1),V.dim(2),V.dim(3))=0;
rand('state',100);
for p=1:V.dim(3),
	img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
	acc = paccuracy(V,p);
	if acc==0,
		udat(:,:,p) = uint8(round((img-mn)*(255/(mx-mn))));
	else,
		% Add random numbers before rounding to reduce aliasing artifact
		r = rand(size(img))*acc;
		udat(:,:,p) = uint8(round((img+r-mn)*(255/(mx-mn))));
	end;
	spm_progress_bar('Set',p);
end;
spm_progress_bar('Clear');
return;

function acc = paccuracy(V,p)
if ~spm_type(V.dim(4),'intt'),
	acc = 0;
else,
	if size(V.pinfo,2)==1,
		acc = abs(V.pinfo(1,1));
	else,
		acc = abs(V.pinfo(1,p));
	end;
end;
%_______________________________________________________________________

%_______________________________________________________________________
function V = smooth_uint8(V,fwhm)
% Convolve the volume in memory (fwhm in voxels).
s  = fwhm/sqrt(8*log(2));
x  = round(6*s(1)); x = [-x:x];
y  = round(6*s(2)); y = [-y:y];
z  = round(6*s(3)); z = [-z:z];
x  = exp(-x.^2/(2*s(1).^2+eps));
y  = exp(-y.^2/(2*s(2).^2+eps));
z  = exp(-z.^2/(2*s(3).^2+eps));
x  = x/sum(x);
y  = y/sum(y);
z  = z/sum(z);

i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;
spm_conv_vol(V.uint8,V.uint8,x,y,z,-[i j k]);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function display_results(VG,VF,x,flags)
fig = spm_figure('FindWin','Graphics');
if isempty(fig), return; end;
set(0,'CurrentFigure',fig);
spm_figure('Clear','Graphics');

txt = 'Information Theoretic Coregistration';
switch lower(flags.cost_fun)
	case 'mi',  txt = 'Mutual Information Coregistration';
	case 'ecc', txt = 'Entropy Correlation Coefficient Registration';
	case 'nmi', txt = 'Normalised Mutual Information Coregistration';
	case 'ncc', txt = 'Normalised Cross Correlation';
	otherwise, error('Invalid cost function specified');
end;

% Display text
%-----------------------------------------------------------------------
ax = axes('Position',[0.1 0.8 0.8 0.15],'Visible','off','Parent',fig);
text(0.5,0.7, txt,'FontSize',16,...
	'FontWeight','Bold','HorizontalAlignment','center','Parent',ax);

Q = inv(VF.mat\spm_matrix(x(:)')*VG.mat);
text(0,0.5, sprintf('X1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(1,:)),'Parent',ax);
text(0,0.3, sprintf('Y1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(2,:)),'Parent',ax);
text(0,0.1, sprintf('Z1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(3,:)),'Parent',ax);

% Display joint histograms
%-----------------------------------------------------------------------
ax  = axes('Position',[0.1 0.5 0.35 0.3],'Visible','off','Parent',fig);
H   = spm_hist2(VG.uint8,VF.uint8,VF.mat\VG.mat,[1 1 1]);
tmp = log(H+1);
image(tmp*(64/max(tmp(:))),'Parent',ax');
set(ax,'DataAspectRatio',[1 1 1],...
	'PlotBoxAspectRatioMode','auto','XDir','normal','YDir','normal',...
	'XTick',[],'YTick',[]);
title('Original Joint Histogram','Parent',ax);
xlabel(spm_str_manip(VG.fname,'k22'),'Parent',ax);
ylabel(spm_str_manip(VF.fname,'k22'),'Parent',ax);

H   = spm_hist2(VG.uint8,VF.uint8,VF.mat\spm_matrix(x(:)')*VG.mat,[1 1 1]);
ax  = axes('Position',[0.6 0.5 0.35 0.3],'Visible','off','Parent',fig);
tmp = log(H+1);
image(tmp*(64/max(tmp(:))),'Parent',ax');
set(ax,'DataAspectRatio',[1 1 1],...
	'PlotBoxAspectRatioMode','auto','XDir','normal','YDir','normal',...
	'XTick',[],'YTick',[]);
title('Final Joint Histogram','Parent',ax);
xlabel(spm_str_manip(VG.fname,'k22'),'Parent',ax);
ylabel(spm_str_manip(VF.fname,'k22'),'Parent',ax);

% Display ortho-views
%-----------------------------------------------------------------------
spm_orthviews('Reset');
h1 = spm_orthviews('Image',VG.fname,[0.01 0.01 .48 .49]);
h2 = spm_orthviews('Image',VF.fname,[.51 0.01 .48 .49]);
global st
st.vols{h2}.premul = inv(spm_matrix(x(:)'));
spm_orthviews('Space',h1);

spm_print
return;

