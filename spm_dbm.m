function morphin
% Deformation Based Morphometry
%
% **** Extract Shapes *****
% Extracts the brain shapes from the nonlinear deformations. 
%
% The deformation fields are defined by both non-linear and  affine
% components.   In  order  to proceed, it is necessary to decompose
% the transformation into components relating to position and  size
% (uninteresting components), and shape (the components that we are
% interested in).  In order  to  effect  this  decomposition,  each
% deformation  field  is  reconstructed from it's parameters.  Each
% field provides a mapping from points in the template to points in
% the  image, allowing standard landmark based registration methods
% to be used  to  extract  the  size  and  positional  information.
% Rather  than  basing the registration on a few landmarks, all the
% elements of the deformation field corresponding to voxels  within
% the brain are considered.
%
%
% **** Run Stats ****
% For the  multivariate  analysis,   it is  necessary to reduce the
% number of  deformation  coefficients  relative to  the  number of
% images.   Principal  component analysis is  used to compact  this
% information.
%
% The analysis is performed by Multivariate Analysis of Covariance. 
%
%_______________________________________________________________________
% %W% John Ashburner %E%


[Finter,Fgraph,CmdLine] = spm('FnUIsetup','DBM',1);
spm_help('!ContextHelp',mfilename);
opt = 1;
while opt ~= 0,
	opt = spm_input('Option?',1,'m',['Extract shapes|Run Stats|Show Difference|Return'], [1 2 3 0]);
	if opt == 1,
		xtract_shapes;
	elseif opt==2,
		run_stats;
	elseif opt==3,
		show_difference;
	end;
end;
fg = spm_figure('FindWin','Interactive');
if ~isempty(fg), spm_figure('Clear',fg); end;
return;
%=======================================================================


%=======================================================================
function xtract_shapes
Q = spm_input('Shape filename',1,'s');
q = max([find(Q == filesep) 0]);
Q = [spm_str_manip(Q((q+1):length(Q)),'sd') '_shape.mat'];

matfiles = spm_get(Inf,'_sn3d.mat',['Normalization Parameters'],'');
drawnow;
n = size(matfiles,1);
Transforms = zeros(8*8*8*3,n);
spm_progress_bar('Init',n,'Extracting Shapes')
for j=1:n,
	T = procrustes(matfiles(j,:));
	Transforms(:,j) = T(:);
	spm_progress_bar('Set',j);
end;
save(Q,'Transforms','matfiles');
spm_progress_bar('Clear');
return;
%=======================================================================

%=======================================================================
function run_stats
fnames = spm_get(Inf,'_shape.mat',['Shape Files'],'');
T = [];
D = [];
for j=1:size(fnames,1),
	load(deblank(fnames(j,:)));
	D1 = zeros(size(Transforms,2),size(fnames,1));
	D1(:,j) = 1;
	T = [T Transforms];
	D = [D ; D1];
end;

ok=0;
while ~ok,
	nn = abs(spm_input('Number of PCs','+1','e',num2str(floor(size(T,2)/3))));
	if prod(size(nn)) == 1,
		if nn > 0 & nn <= size(T,2)/3 & nn<= size(T,1), ok=1; end;
	end;
end;

[C,G]=enter_design_matrix(D);

if spm_input('Restrict to asymmetries?','+1','y/n',[1 0],2),
	d2 = size(T,2);
	tmp = reshape(T,[8 8 8 3 d2]);
	tmp(2:2:8,:,:,1,:)=0;
	tmp(1:2:7,:,:,2:3,:)=0;
	T = reshape(tmp,[8*8*8*3 d2]);
end;

fg = spm_figure('FindWin','Graphics');
spm_clf(fg);


fprintf('\nPCA and saving %d coefficients - ',nn);
[v1,s1,u1]=svd(T,0);
X=u1(:,1:nn)*s1(1:nn,1:nn);
fprintf('accounting for %g of the variance\n\n',...
	sum(diag(s1(1:nn,1:nn)).^2)/sum(diag(s1).^2));

%text(0.30,1.02,'Deformation Based Morphometry','Fontsize',16,'Fontweight','Bold','Parent',fg);
%text(0.30,1.02,'Design Matrix','Fontsize',16,'Fontweight','Bold','Parent',fg);
ax=axes('Position',[0.2 0.8 0.6 0.1],'Parent',fg,'Visible','off');
text(0,1,'Deformation Based Morphometry','Fontsize',18,'Fontweight','Bold','Parent',ax);

ax=axes('Position',[0.2 0.6 0.6 0.2],'Parent',fg);
s2 = diag(s1).^2;
plot(cumsum(s2)/sum(s2),'o','Parent',ax);
set(ax,'XLim',[0 length(s2)],'YLim',[0 1],'FontSize',14);
set(get(ax,'XLabel'),'String','Number of principal components','FontSize',14);
set(get(ax,'YLabel'),'String','Fraction of variance','FontSize',14);
set(get(ax,'Title'),'String','Principal Components','FontSize',16,'Fontweight','Bold');
text(0.2,0.1,...
	sprintf('%d coefficients account for %.4g%% of variance',...
	nn,100*sum(diag(s1(1:nn,1:nn)).^2)/sum(diag(s1).^2)),...
	'Parent',ax,'FontSize',14);

ax   = axes('Position',[0.2 0.2 0.6 0.3],'Parent',fg);
dmtx = [C G];
tmp  = min(dmtx).*any(diff(dmtx));
dmtx = dmtx - repmat(tmp,size(dmtx,1),1);
dmtx = 64*(dmtx./repmat(max(dmtx),size(dmtx,1),1));
image(dmtx','Parent',ax);
set(ax,'NextPlot', 'add');
set(get(ax,'XLabel'),'String','Subjects','FontSize',14);
set(get(ax,'YLabel'),'String','Conditions & Covariates','FontSize',14);
set(get(ax,'Title'),'String','Design Matrix','Fontsize',16,'Fontweight','Bold');
plot([0 size(dmtx,1)+0.5]',[1 1]'*(size(C,2)+0.5),'r-','Parent',ax,'LineWidth',2);
set(ax,'NextPlot', 'replace','YTick',1:size(dmtx,2));

[P,Lambda,chi2,df] = mancova(X,C,G);

ax = axes('Position',[0.2 0.0 0.6 0.2],'Parent',fg,'Visible','off');
text(0,0.6,...
	['Wilke''s \lambda=' num2str(Lambda) ',  \chi^2=' num2str(chi2)],...
	'Fontsize',16,'Fontweight','Bold','Parent',ax,'Interpreter','tex');
text(0,0.3,...
	['\nu=' num2str(df) ',  p=' num2str(P) ],...
	'Fontsize',16,'Fontweight','Bold','Parent',ax,'Interpreter','tex');

fprintf('Lambda= %g \t chi2= %g \t df= %g \t P= %g\n', Lambda, chi2, df, P);

if P<0.05,
	disp('This looks like a result!');
end;
return;
%=======================================================================

%=======================================================================
function [C,G]=enter_design_matrix(D)
ok = 0;
while ~ok, % (size(tmp,2) ~= size(D,2)) & ~isempty(tmp),
	tmp=spm_input(['Interesting Group weights (' num2str(size(D,2)) ')'],'+1','e');
	if size(tmp,2) == size(D,2), C  = D*tmp'; ok = 1; end;
	if isempty(tmp), C = []; ok = 1; end;
end;

ok = 0;
while ~ok,
	tmp=spm_input(['Interesting Covariates (' num2str(size(D,1)) ')'],'+1','e');
	if (size(tmp,2) == size(D,1)) | isempty(tmp), ok = 1; end;
end;
C = [C tmp'];
if ~isempty(C), C = C - repmat(mean(C),size(C,1),1); end;

ok = 0;
while ~ok,
	tmp=spm_input(['Uninteresting Group weights (' num2str(size(D,2)) ')'],'+1','e');
	if size(tmp,2) == size(D,2), G  = D*tmp'; ok = 1; end;
	if isempty(tmp), G = []; ok = 1; end;
	
end;

ok = 0;
while ~ok,
	tmp=spm_input(['Uninteresting Covariates (' num2str(size(D,1)) ')'],'+1','e');
	if (size(tmp,2) == size(D,1)) | isempty(tmp), ok = 1; end;
end;
G = [G tmp'];
if ~isempty(G), G = G - repmat(mean(G),size(G,1),1); end;
G = [G ones(size(D,1),1)];
return;
%=======================================================================

%=======================================================================
function [P,Lambda,chi2,df]=mancova(X1,C,G)
[n,b]  = size(C);
b      = rank(C);
p      = size(X1,2);
w      = n-b-rank(G);
if ~isempty(G),
	XA     = X1 - G*(pinv(G)*X1);
	CA     = C  - G*(pinv(G)*C);
else,
	XA = X1;
	CA = C;
end;
T      = CA*(pinv(CA)*XA);
B0     = T'*T;
tmp    = XA - T;
W0     = tmp'*tmp;

%Lambda = det(W0)/det(W0+B0);
Lambda = exp(sum(log(svd(W0)))-sum(log(svd(W0+B0))));

chi2   = -(w-(p-b+1)/2)*log(Lambda);
df     = p*b;
P      = (1-spm_Xcdf(chi2,df));
return;
%=======================================================================

%=======================================================================
function [v,u]=cva(X1,C,G)
[n,b]   = size(C);
b       = rank(C);
p       = size(X1,2);
w       = n-b-rank(G);
if ~isempty(G),
	XA     = X1 - G*(pinv(G)*X1);
	CA     = C  - G*(pinv(G)*C);
else,
	XA = X1;
	CA = C;
end;
T       = CA*(pinv(CA)*XA);
B0      = T'*T;
tmp     = XA - T;
W0      = tmp'*tmp;
B       = B0/(b);
W       = W0/(n-b-size(G,2));
[v,d]   = eig(W\B);
[tmp,q] = sort(-real(abs(diag(d))));
b       = rank(B);
d       = real(d(q(1:b),q(1:b)));
v       = real(v(:,q(1:b)));
v       = v(:,1:b)*d;
u       = XA*v;
if size(C,2)==1,
	if u'*C	< 0,
		v  = -v;
		u  = -u;
	end;
end;
return;
%=======================================================================

%=======================================================================
function Transform = procrustes(matfile)
load(deblank(matfile));

global sptl_Ornt
tmp = spm_matrix(sptl_Ornt);
if det(tmp(1:3,1:3))<0,
	% Need to mirror the co-ordinate system.
	tmp = diag([-1 1 1 1]);
else,
	tmp = diag([ 1 1 1 1]);
end;
Mat = (tmp*MF)*Affine*inv(MG);

global SWD;
V = spm_vol(fullfile(SWD,'apriori','brainmask.img'));

Vox = [2 2 2];
bb  = reshape([-90 91 -126 91  -72 109],2,3);

x = (bb(1,1):Vox(1):bb(2,1))/Dims(3,1) + Dims(4,1);
y = (bb(1,2):Vox(2):bb(2,2))/Dims(3,2) + Dims(4,2);
z = (bb(1,3):Vox(3):bb(2,3))/Dims(3,3) + Dims(4,3);

Dim = [length(x) length(y) length(z)];
origin = round(-bb(1,:)./Vox + 1);

basX = spm_dctmtx(Dims(1,1),Dims(2,1),x-1);
basY = spm_dctmtx(Dims(1,2),Dims(2,2),y-1);
basZ = spm_dctmtx(Dims(1,3),Dims(2,3),z-1);

TX = zeros(Dim);
TY = zeros(Dim);
TZ = zeros(Dim);

X = x'*ones(1,Dim(2));
Y = ones(Dim(1),1)*y;

M = V.mat\Mat*V.mat;

% Reconstruct deformation field
for j=1:Dims(1,3),
	tx = reshape( reshape(Transform(:,1),Dims(2,1)*Dims(2,2),Dims(2,3)) *basZ(j,:)', Dims(2,1), Dims(2,2) );
	ty = reshape( reshape(Transform(:,2),Dims(2,1)*Dims(2,2),Dims(2,3)) *basZ(j,:)', Dims(2,1), Dims(2,2) );
	tz = reshape( reshape(Transform(:,3),Dims(2,1)*Dims(2,2),Dims(2,3)) *basZ(j,:)', Dims(2,1), Dims(2,2) );

	ttx = X + basX*tx*basY';
	tty = Y + basX*ty*basY';
	ttz = z(j) + basX*tz*basY';

	TX(:,:,j) = M(1,1)*ttx + M(1,2)*tty + M(1,3)*ttz + M(1,4);
	TY(:,:,j) = M(2,1)*ttx + M(2,2)*tty + M(2,3)*ttz + M(2,4);
	TZ(:,:,j) = M(3,1)*ttx + M(3,2)*tty + M(3,3)*ttz + M(3,4);
end;

% Compute centres of mass.
sm = 0;
c1 = [0 0 0];
c2 = [0 0 0];
[Y,X]=meshgrid(1:Dim(2),1:Dim(1));
for j=1:Dims(1,3),
	msk = spm_slice_vol(V,spm_matrix([0 0 j]),V.dim(1:2),0);
	X1  = TX(:,:,j);
	Y1  = TY(:,:,j);
	Z1  = TZ(:,:,j);
	sm = sm + sum(msk(:));
	c1  = c1 + sum([ X(:).*msk(:)   Y(:).*msk(:)      j*msk(:)]);
	c2  = c2 + sum([X1(:).*msk(:)  Y1(:).*msk(:) Z1(:).*msk(:)]);
end;
c1 = c1/sm;
c2 = c2/sm;

mom1 = 0;
mom2 = 0;
% Compute centroid size.
for j=1:Dims(1,3),
	msk = spm_slice_vol(V,spm_matrix([0 0 j]),V.dim(1:2),0);
	X1  = TX(:,:,j);
	Y1  = TY(:,:,j);
	Z1  = TZ(:,:,j);
	sm = sm + sum(msk(:));
	mom1 = mom1 + sum([( X(:)-c1(1)).^2.*msk(:)  ( Y(:)-c1(2)).^2.*msk(:) (    j-c1(3)).^2 *msk(:)]);
	mom2 = mom2 + sum([(X1(:)-c2(1)).^2.*msk(:)  (Y1(:)-c2(2)).^2.*msk(:) (Z1(:)-c2(3)).^2.*msk(:)]);
end;
zm = sqrt(mom1/mom2);
%zm=1;
ZM = spm_matrix([0 0 0  0 0 0  zm zm zm  0 0 0]);


% Compute rotations
C = zeros(3);
for j=1:Dims(1,3),
	msk = spm_slice_vol(V,spm_matrix([0 0 j]),V.dim(1:2),0);
	X1  = (TX(:,:,j) -c2(1))*zm;
	Y1  = (TY(:,:,j) -c2(2))*zm;
	Z1  = (TZ(:,:,j) -c2(3))*zm;
	C = C + [(X(:)-c1(1)).*msk(:) (Y(:)-c1(2)).*msk(:) (j-c1(3))*msk(:)]'*[X1(:).*msk(:) Y1(:).*msk(:) Z1(:).*msk(:)];
end;
[u,s,v]=svd(C);
R=eye(4);
R(1:3,1:3) = u*v';

% Rotate and translate deformation field.
T = spm_matrix(c1)*R*ZM*spm_matrix(-c2);
for j=1:Dims(1,3),
	X1 = TX(:,:,j);
	Y1 = TY(:,:,j);
	Z1 = TZ(:,:,j);
	TX(:,:,j) = T(1,1)*X1 + T(1,2)*Y1 + T(1,3)*Z1 + T(1,4);
	TY(:,:,j) = T(2,1)*X1 + T(2,2)*Y1 + T(2,3)*Z1 + T(2,4);
	TZ(:,:,j) = T(3,1)*X1 + T(3,2)*Y1 + T(3,3)*Z1 + T(3,4);
end;

basX = spm_dctmtx(length(x),8);
basY = spm_dctmtx(length(y),8);
basZ = spm_dctmtx(length(z),8);

dx = size(basX);
dy = size(basY);
dz = size(basZ);

Transform = zeros(size(basX,2)*size(basY,2)*size(basZ,2),3);
for j=1:Dims(1,3),
	tx = reshape( reshape(Transform(:,1),dx(2)*dy(2),dz(2)) *basZ(j,:)', dx(2), dy(2));
	ty = reshape( reshape(Transform(:,2),dx(2)*dy(2),dz(2)) *basZ(j,:)', dx(2), dy(2));
	tz = reshape( reshape(Transform(:,3),dx(2)*dy(2),dz(2)) *basZ(j,:)', dx(2), dy(2));

	ttx = basX'*(TX(:,:,j)-X)*basY; ttx = ttx(:)*basZ(j,:);
	tty = basX'*(TY(:,:,j)-Y)*basY; tty = tty(:)*basZ(j,:);
	ttz = basX'*(TZ(:,:,j)-j)*basY; ttz = ttz(:)*basZ(j,:);
	Transform(:,1) = Transform(:,1) + ttx(:);
	Transform(:,2) = Transform(:,2) + tty(:);
	Transform(:,3) = Transform(:,3) + ttz(:);
end;
%=======================================================================

%=======================================================================
function pt=setup_orthviews(Transform)
Transform=-Transform;

Dims = [
   91.0000  109.0000   91.0000
    8.0000    8.0000    8.0000
    2.0000    2.0000    2.0000
   46.0000   64.0000   37.0000];

global SWD;
V = spm_vol(fullfile(SWD,'apriori','gray.img'));

Vox = [2 2 2];
bb = reshape([-90 91 -126 91  -72 109],2,3);

x = (bb(1,1):Vox(1):bb(2,1))/Dims(3,1) + Dims(4,1);
y = (bb(1,2):Vox(2):bb(2,2))/Dims(3,2) + Dims(4,2);
z = (bb(1,3):Vox(3):bb(2,3))/Dims(3,3) + Dims(4,3);

Dim = [length(x) length(y) length(z)];
origin = round(-bb(1,:)./Vox + 1);

basX = spm_dctmtx(Dims(1,1),Dims(2,1),x-1);
basY = spm_dctmtx(Dims(1,2),Dims(2,2),y-1);
basZ = spm_dctmtx(Dims(1,3),Dims(2,3),z-1);

TX = zeros(Dims(1,:));
TY = zeros(Dims(1,:));
TZ = zeros(Dims(1,:));

X = x'*ones(1,Dim(2));
Y = ones(Dim(1),1)*y;


vr = 0;
Transform = reshape(Transform,prod(Dims(2,:)),3);
spm_progress_bar('Init',Dims(1,3),'Setting up orthviews');
for j=1:Dims(1,3),
	tx = reshape(reshape(Transform(:,1),Dims(2,1)*Dims(2,2),Dims(2,3))*basZ(j,:)', Dims(2,1), Dims(2,2));
	ty = reshape(reshape(Transform(:,2),Dims(2,1)*Dims(2,2),Dims(2,3))*basZ(j,:)', Dims(2,1), Dims(2,2));
	tz = reshape(reshape(Transform(:,3),Dims(2,1)*Dims(2,2),Dims(2,3))*basZ(j,:)', Dims(2,1), Dims(2,2));
	TX(:,:,j) = basX*tx*basY';
	TY(:,:,j) = basX*ty*basY';
	TZ(:,:,j) = basX*tz*basY';
	spm_progress_bar('Set',j);
end;
displ2 = sum((Vox(1)*TX(:)).^2 + (Vox(2)*TY(:)).^2 + (Vox(3)*TZ(:)).^2,2);
sc = 8/sqrt(sum(displ2)/prod(size(TX)));
for j=1:Dims(1,3),
	TX(:,:,j) =    X + TX(:,:,j)*sc;
	TY(:,:,j) =    Y + TY(:,:,j)*sc;
	TZ(:,:,j) = z(j) + TZ(:,:,j)*sc;
	spm_progress_bar('Set',j);
end;
spm_progress_bar('clear');

wt = zeros(V.dim(1:3));
for j=1:Dims(1,3),
	wt(:,:,j)=spm_slice_vol(V,spm_matrix([0 0 j]),V.dim(1:2),0);
end;

[tmp,pt] = max(displ2.*wt(:));
pt=pt(1);
pt = [rem(pt-1,size(TX,1))+1  rem(floor((pt-1)/size(TX,1)),size(TX,2))+1  floor((pt-1)/(size(TX,1)*size(TX,2)))+1]';
pt = V.mat(1:3,1:3)*pt + V.mat(1:3,4);

global dbm_orth;
dbm_orth = struct('TX',TX,'TY',TY,'TZ',TZ,'V',V,'sc',sc);

return;
%=======================================================================

%=======================================================================
function sample_orthviews(tpt)
global dbm_orth;
M = inv(dbm_orth.V.mat);

pt = M(1:3,1:3)*reshape(tpt,3,1) + M(1:3,4);

Dims = [
   91.0000  109.0000   91.0000
    8.0000    8.0000    8.0000
    2.0000    2.0000    2.0000
   46.0000   64.0000   37.0000];

Vox = [2 2 2];
bb = reshape([-90 91 -126 91  -72 109],2,3);

x = (bb(1,1):Vox(1):bb(2,1))/Dims(3,1) + Dims(4,1);
y = (bb(1,2):Vox(2):bb(2,2))/Dims(3,2) + Dims(4,2);
z = (bb(1,3):Vox(3):bb(2,3))/Dims(3,3) + Dims(4,3);

Dim = [length(x) length(y) length(z)];
origin = round(-bb(1,:)./Vox + 1);

Dim = size(dbm_orth.TX);
X = x'*ones(1,Dim(2));
Y = ones(Dim(1),1)*y;
smp = 8;

[X,Y]= meshgrid(1:Dim(1),1:Dim(2));
Z    = ones(size(X))*pt(3);
X1   = spm_sample_vol(dbm_orth.TX,X,Y,Z,0);
Y1   = spm_sample_vol(dbm_orth.TY,X,Y,Z,0);
Z1   = spm_sample_vol(dbm_orth.TZ,X,Y,Z,0);
imgt = spm_sample_vol(dbm_orth.V,X1,Y1,Z1,1);
mskt = spm_sample_vol(dbm_orth.V,X,Y,Z,1);
msk  = find(mskt & rem(X,smp)==0 & rem(Y,smp)==0);
Xt   = X(msk);
Yt   = Y(msk);
dXt  = X(msk)-X1(msk);
dYt  = Y(msk)-Y1(msk);

[Y,Z]= meshgrid(1:Dim(2),1:Dim(2));
X    = ones(size(Y))*pt(1);
X1   = spm_sample_vol(dbm_orth.TX,X,Y,Z,0);
Y1   = spm_sample_vol(dbm_orth.TY,X,Y,Z,0);
Z1   = spm_sample_vol(dbm_orth.TZ,X,Y,Z,0);
imgs = spm_sample_vol(dbm_orth.V,X1,Y1,Z1,1);
msks = spm_sample_vol(dbm_orth.V,X,Y,Z,1);
msk  = find(msks & rem(Y,smp)==0 & rem(Z,smp)==0);
Xs   = Y(msk);
Ys   = Z(msk);
dXs  = Y(msk)-Y1(msk);
dYs  = Z(msk)-Z1(msk);

[X,Z]= meshgrid(1:Dim(1),1:Dim(2));
Y    = ones(size(X))*pt(2);
X1   = spm_sample_vol(dbm_orth.TX,X,Y,Z,0);
Y1   = spm_sample_vol(dbm_orth.TY,X,Y,Z,0);
Z1   = spm_sample_vol(dbm_orth.TZ,X,Y,Z,0);
imgc = spm_sample_vol(dbm_orth.V,X1,Y1,Z1,1);
mskc = spm_sample_vol(dbm_orth.V,X,Y,Z,1);
msk  = find(mskc & rem(X,smp)==0 & rem(Z,smp)==0);
Xc   = X(msk);
Yc   = Z(msk);
dXc  = X(msk)-X1(msk);
dYc  = Z(msk)-Z1(msk);

fg = spm_figure('FindWin','Graphics');
set(0,'CurrentFigure',fg);
ax=subplot(2,1,2);
image([imgt imgs imgc]*64*0.8,'Parent',ax);
set(ax,'NextPlot','add');
set(0,'CurrentFigure',fg);
set(fg,'CurrentAxes',ax);
contour([mskt msks mskc],1,'y');
q=quiver(	[Xt ; Xs+Dim(1) ; Xc+Dim(1)+Dim(2)],...
		[Yt ; Ys ; Yc], ...
		[dXt ; dXs ; dXc],...
		[dYt ; dYs ; dYc],0);
set(q,'Color',[1 1 0],'Parent',ax);
tmp = plot(pt([1 2 1])+[0 Dim(1) Dim(1)+Dim(2)]', pt([2 3 3]),'+','Parent',ax);
set(tmp,'MarkerSize',20,'Color',[1 1 0]);

hold off;
axis image off xy
set(get(ax,'Xlabel'),'String',['scaled x ' num2str(dbm_orth.sc) ',  coordinate=[' num2str(reshape(tpt,1,3)) ']'],...
	'FontSize',14,'Visible','on');
drawnow;
return;
%=======================================================================

%=======================================================================
function show_difference
fnames = spm_get(Inf,'_shape.mat',['Shape Files'],'');
T = [];
D = [];
for j=1:size(fnames,1),
	load(deblank(fnames(j,:)));
	D1 = zeros(size(Transforms,2),size(fnames,1));
	D1(:,j) = 1;
	T = [T Transforms];
	D = [D ; D1];
end;
[C,G]=enter_design_matrix(D);

if spm_input('Restrict to asymmetries?','+1','y/n',[1 0],2),
	d2 = size(T,2);
	tmp = reshape(T,[8 8 8 3 d2]);
	tmp(2:2:8,:,:,1,:)=0;
	tmp(1:2:7,:,:,2:3,:)=0;
	T = reshape(tmp,[8*8*8*3 d2]);
end;

opt = spm_input('Type of representation?','+1','m',['Difference Representation|Canonical Variate Representation'], [1 2]);

if opt==1,
	transf = pinv([C G])*T';
	if ~isempty(G),  TA = T' - G*(pinv([G])*T');  else,  TA = T;  end;

	if size(C,2)==1,
		transf = transf(1,:)';
	else,
		disp('Using first column');
		transf = transf(1,:)';
	end;

	plx = C(:,1);
	ply = TA*transf;
	ttl = 'Difference Representation';
else
	ok=0;
	while ~ok,
		nn = abs(spm_input('Number of PCs','+1','e',num2str(floor(size(T,2)/3))));
		if prod(size(nn)) == 1,
			if nn > 0 & nn <= size(T,2)/3 & nn<= size(T,1), ok=1; end;
		end;
	end;

	if size(C,2) ~=1,
		str = '1';
		for i=2:size(C,2),
			str = [str '|' num2str(i)];
		end;
		nnn = spm_input('Plot which vector','+1','b',str,1:size(C,2),1);
		tmp = spm_input('Against which covariate of interest','+1','b',str,1:size(C,2),1);
		plx = [C(:,tmp)];
	else,
		nnn = 1;
		plx = [C(:,1)];
	end;
	fprintf('\nPCA and saving %d coefficients - ',nn);
	[v1,s1,u1]=svd(T,0);
	X=u1(:,1:nn)*s1(1:nn,1:nn);
	fprintf('accounting for %g of the variance\n\n',...
		sum(diag(s1(1:nn,1:nn)).^2)/sum(diag(s1).^2));
	[tmp1,tmp2] = cva(X,C,G);
	transf = (v1(:,1:nn)*tmp1);
	transf = transf(:,nnn);
	ply = [tmp2];
	ttl = 'Canonical Variate Representation';
end;

pt = setup_orthviews(transf);

fg = spm_figure('FindWin','Graphics');
spm_clf(fg);
set(0,'CurrentFigure',fg);
ax=axes('Parent',fg,'Position',[.3 .5 .5 .3]);
pl=plot(plx,ply,'o','Parent',ax);
dd=diff(get(gca,'Xlim'))/10;
set(pl,'MarkerSize',6,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0], 'Marker','o');
set(ax,'Xlim',get(gca,'Xlim')+[-dd dd]);
set(get(ax,'Title'),'String',ttl);
set(get(ax,'Xlabel'),'String','Variate of interest');


sample_orthviews(pt);
opt = 1;
while opt ~= 0,
	opt = spm_input('Option?',1,'m',['Change co-ordinate|Return'], [1 0]);
	if opt == 1,
		ctr = [];
		while prod(size(ctr))~=3,
			ctr = spm_input('Co-ordinate',2,'e');
		end;
		sample_orthviews(ctr);
	elseif opt == 0,
		clear global dbm_orth
	end;
end;

return;
