function spm_fmri_mancova
% MAN[C]OVA and CVA for fMRI time-series
% FORMAT spm_fmri_mancova
%___________________________________________________________________________
%
%
% spm_fmri_mancova effects a multivariate analysis of fMRI time-series that
% characterizes spatio-temporal response patterns in terms of (i) a p-value
% (ii) a series of canonical images and (iii) a corresponding set of
% time-dependent variates reflecting the expression of the canonical
% images (these are expressed as a function of time over the entire
% series and a hemodynamic transient associated with each condition).
%
% It is assumed that the data represent [approximately] independent
% observations.
%
% The data that are subject to MANCOVA are the adjusted data matrices
% that obtain after removing any confounds using ANCOVA on a voxel by voxel
% basis (e.g. time and global activity). The data are then reduced using 
% the most prevalent spatial modes.
%
% The effects of interest are modeled as evoked hemodynamic transients
% that can be differemt in form, from voxel to voxel, and from condition
% to condition.  The form of these transients is estimated in the context
% of the general linear model using Fourier basis functions.
%
% spm_fmri_mancova saves the canonical images, values and variates in
% CVA.mat. These spatial modes can be displayed using spm_fmri_cva.
%
% Ref: Friston et al (1995).  Characterizing brain dynamics with fMRI:
% A multivariate approach. NeuroImage 0:00-00
%
%__________________________________________________________________________
% %W% %E%
 

% get filenames and create design matrix
%---------------------------------------------------------------------------
clear; global CWD; cd(CWD)

figure(2); clf; set(2,'Name','MAN[C]OVA - Canonical analysis');

Q     = [];				% matrix of filename strings
X     = [];				% data matrix {q x pixels}
H     = [];				% condition effects
G     = [];				% confounds

% cycle over studies, subjects and conditions
%----------------------------------------------------------------------------
n      = spm_input(['# of subjects'],1);
k      = spm_input(['# of scans/subject'],2);
q      = n*k;
for j  = 1:n
	d = ['select scans 1 - ' num2str(k) ' {subject ' num2str(j) '}'];
 	Q = str2mat(Q,spm_get(k,'.img',d));
end
Q(1,:) = [];
P      = Q;

% design matrix subpartitions - effects of interest
%===========================================================================
h     = spm_input('number of conditions', 3);
l     = spm_input('# scans per condition',4);

% set basis functions for transients
%---------------------------------------------------------------------------
W     = [];
for i = 1:4
	d     = sin(i*pi*[0:(l + 1)]/(l + 1));
	W     = [W d(:)];
end
for i = 1:h
	d       = spm_input(sprintf('onsets of condition %0.0f {scans}',i),5);
	c       = zeros(q,1);
	c(d(:)) = ones(size(d(:)));
	for j = 1:size(W,2)
		d = conv(c,W(:,j));
		d = d(1:q);
		H = [H d(:)];
	end
end


% design matrix subpartitions - confounds
%===========================================================================

% subject effects
%---------------------------------------------------------------------------
for i = 1:k; G = [G; eye(n)]; end
G     = reshape(G,n,q)';
g     = size(G,2);
if (g == q) | (g == 1); G = []; end

% condition subject subpartition
%---------------------------------------------------------------------------
if ~isempty(G)
	d      = G(:,1);
	G(:,1) = [];
	G      = G - d*sum(G)/sum(d);
end

% high pass filter
%----------------------------------------------------------------------------
if spm_input('high pass filter',5,'yes|no',[1 0])
	for v = 0.5:3
		d = sin(2*pi*[1:q]*v/q);
		G = [G d(:)];
		d = cos(2*pi*[1:q]*v/q);
		G = [G d(:)];
	end
end

% column of ones and mean correct H (to ensure orthogonality of G and H)
%----------------------------------------------------------------------------
G      = [G ones(q,1)]; end
H      = H - ones(q,1)*mean(H);


% get Repeat time
%----------------------------------------------------------------------------
RT     = spm_input('Repeat time {secs}',6);



% interactive components now complete
%============================================================================
set(2,'Name','AnCova','Pointer','Watch')


% get file identifiers and global values
%----------------------------------------------------------------------------
V     = zeros(12,q);
for i = 1:q; V(:,i) = spm_map(P(i,:));  end


% check for consistency of image size and voxel size
%----------------------------------------------------------------------------
if ~(all(all(~diff(V([1:3],:)'))))
	error('images are not compatible'); end


% get ORIGIN
%----------------------------------------------------------------------------
[DIM VOX SCALE TYPE OFFSET ORIGIN] = spm_hread(P(1,:));
if DIM(3) == 1; ORIGIN = [1 1 1]; end

% adjust scaling coefficients so that Grand mean - <m> = 100
%----------------------------------------------------------------------------
m      = zeros(q,1);
for i  = 1:q
	m(i) = spm_global(V(:,i)); end

V(7,:) = V(7,:)*100/mean(m);
m      = m*100/mean(m);
G      = [G m];


% degrees of freedom (only G is used in ANCOVA and then H in MANCOVA)
%---------------------------------------------------------------------------
h     = rank(H);					% condition
g     = rank(G);					% confounds
r     = q - rank(G);					% residuals [AnCova]

% vextors of location
%---------------------------------------------------------------------------
[y x] = meshgrid([1:V(2,1)],[1:V(1,1)]');
z     = ones(size(x));
x     = (x(:) - ORIGIN(1))*V(4,1);
y     = (y(:) - ORIGIN(2))*V(5,1);
z     = z(:);


% delete previous results
%---------------------------------------------------------------------------
cd(CWD)
delete XA.mat
delete XYZ.mat
delete SPMF.mat

% set thresholds
%----------------------------------------------------------------------------
UF    = spm_invFcdf(1 - 0.05,[g,r]);				% for SPM{F}
GX    = m*ones(1,V(1,1)*V(2,1));				% global 

% cycle over planes to avoid working memory problems
%---------------------------------------------------------------------------
for i   = 1:V(3,1)

  X     = zeros(q,V(1,1)*V(2,1));				% data matrix
  for j = 1:q
      d = spm_slice_vol(V(:,j),spm_matrix([0 0 i]),[V(1,1) V(2,1)],0);
      X(j,:) = d(:)';
  end

  % eliminate background voxels (less than 0.8 of the global activity GX)
  %-------------------------------------------------------------------------
  Q     = find(all(X > 0.8*GX) & ~all(~diff(X)));

  if length(Q)
    X     = X(:,Q);
    XYZ   = [x(Q) y(Q) z(Q)*(i - ORIGIN(3))*V(6,1)]';		% locatinns

    % estimate parameters and sum of squares due to error
    %-----------------------------------------------------------------------
    BETA  = pinv(G)*X;
    R     = sum((X - G*BETA).^2);

    % test for effects of interest SPM{F} with he F ratio of variances
    %-----------------------------------------------------------------------
    N     = sum(X.^2);
    F     = (N - R)/g./(R/r);

    % adjustment (remove effects of no interest)
    %-----------------------------------------------------------------------
    XA    = X - G*BETA;

    % remove non-significant voxels from further analysis
    %-----------------------------------------------------------------------
    Q     = find(F > UF);

    % cumulate significant voxels
    %-----------------------------------------------------------------------
    if length(Q)
    	spm_append('XA'  ,XA(:,Q));
    	spm_append('SPMF',F(Q) );
    	spm_append('XYZ' ,XYZ(:,Q));
     end				% end conditional on P(F > x) < 0.05
  end					% end conditional on non-zero voxels
end					% end cycle over planes


% unmap volumes and retain image dimensions
%---------------------------------------------------------------------------
for i  = 1:q; spm_unmap(V(:,i)); end
V      = [V(1:6,1); ORIGIN(:)];


% display and print SPM{F}, Design matrix and textual information
%---------------------------------------------------------------------------
figure(3); spm_clf

load XYZ
load SPMF
axes('Position',[0.1 0.5 0.8 0.4]);
spm_mip(sqrt(SPMF),XYZ,V(1:6))
title('SPM{F} from AnCova (p < 0.05)','FontSize',16,'Fontweight','Bold')

axes('Position',[0.1 0.1 0.36 0.36]);
imagesc(-G*inv(diag(max(abs(G)))))
axis image
title ('Design matrix','Fontweight','Bold')
xlabel 'Confounds'; ylabel 'scan'


axes('Position',[0.46 0.1 0.42 0.3]); axis off
text(0,1,'Results directory:')
text(0,.9,CWD,'FontSize',16,'Fontweight','Bold')
text(0,.7,sprintf('Image size  %-4.0f %-4.0f %-4.0f voxels',V(1:3)))
text(0,.6,sprintf('Voxel size  %-4.1f %-4.1f %-4.1f mm',V(4:6)))
text(0,.4,'Degrees of freedom','Fontsize',10)
text(0,.2,sprintf('Confounds  %0.0f ',g),'Fontsize',10)
text(0,.1,sprintf('Residuals  %0.0f ',r),'Fontsize',10)
line([0 1],[1.1 1.1],'LineWidth',3);
line([0 1],[0.8 0.8],'LineWidth',3);
line([0 1],[0.3 0.3],'LineWidth',1);
line([0 1],[0.0 0.0],'LineWidth',1);

spm_print


% MULTIVARIATE ANALYSIS - ManCova and CVA
%===========================================================================
set(2,'Name','ManCova','Pointer','Watch')
clear GX
load XA
load XYZ


% Euclidean normalize
%---------------------------------------------------------------------------
d    = ones(q,1)*sqrt(sum(XA.^2));
XA   = XA./d;

% SVD DIMENSION REDUCTION
%===========================================================================

% compute SSPM
%---------------------------------------------------------------------------
A     = zeros(q);
for i = 1:q
	A(:,i) = XA*(XA(i,:)'); end

% degrees of freedom
%---------------------------------------------------------------------------
h     = rank(H);					% effects   [ManCova]
r     = rank(A) - h;					% residuals [ManCova]

% eigenvector reduction (with temporal smoothing)
%---------------------------------------------------------------------------
[e s] = eig(A);
[i j] = sort(-diag(s));
s     = s(j,j);
e     = e(:,j);

v     = diag(s);
v     = length(v)*v/sum(v);

% reduce dimensionality
%---------------------------------------------------------------------------
Q     = find(v > 1);
X     = e(:,Q)*s(Q,Q)^(0.5);
p     = size(X,2);					% p-variate repsonse

% MANCOVA
%===========================================================================

% estimate parameters
%---------------------------------------------------------------------------
BETA  = pinv(H)*X;					% parameter estimates
T     = H*BETA;						% condition
R     = X - T;						% error


% CANONICAL VARIATES ANALYSIS
%===========================================================================+

% Identify significant canonical images and variates
% Using the generalized eigenvalue solution
%----------------------------------------------------------------------------
[E v] = eig(T'*T,R'*R);
v     = diag(real(v));
[d j] = sort(-v);
v     = v(j,1);
E     = E(:,j);

% test for significance of H with Wilk's Lambda (i.e P{dimensionality > 0})
%----------------------------------------------------------------------------
CHI   = -(r - ((p - h + 1)/2))*log(det(R'*R)/det(X'*X));


% test for the dimensionality of the alternative hypothesis
%----------------------------------------------------------------------------
pV    = [];
for t = 0:(min([h p]) - 1)
	chi  = (r - ((p - h + 1)/2))*log(prod(1 + v((t + 1):p)));
	pV  = [pV (1 - spm_Xcdf(chi, (p - t)*(h - t)) )];
end


% compute canonical images and varies at p < 0.1
%----------------------------------------------------------------------------
E     = real(E(:,find(pV < 0.05)));
CV    = T*E;					% canonical variate
CU    = XA'*(e(:,Q)*s(Q,Q)^(-0.5)*E);		% canonical images

% save variables for subsequnet analysis
%---------------------------------------------------------------------------
save CVA CV CU H G A X T R v pV CHI V XYZ E BETA W RT


% maximum intensity projections and textual information
%===========================================================================
figure(3); spm_clf


% Design matrix
%---------------------------------------------------------------------------
axes('Position',[0.1 0.56 0.3 0.3]);
imagesc(H)
title ('Design matrix','Fontweight','Bold')
xlabel 'condition effects'; ylabel 'scan'


% textual information
%---------------------------------------------------------------------------
axes('Position',[0.48 0.6 0.5 0.16],'Visible','off'); y = 1;
text(0,y,'Canonical analysis','FontSize',16,'Fontweight','Bold'); y = y - 0.2;
text(0,y,CWD); y = y - 0.2;
text(0,y,sprintf('Chi-squared = %0.2f',CHI)); y = y - 0.2;
text(0,y,sprintf('threshold @ p = 0.05 = %0.2f',spm_invXcdf(0.95,h*p)));
y = y - 0.2; text(0,y,sprintf('df{error} = %0.0i, df{Chi} = %0.0i',r,h*p));
y = y - 0.2; text(0,y,sprintf('p value = %0.4f',pV(1)),'FontSize',16);


% canonical value spectrum
%---------------------------------------------------------------------------
subplot(2,2,3)
[x y] = bar(v*r/h);
fill(x,y,[1 1 1]*0.9);
axis square
set(gca,'XLim',[0 (1 + min([h p]))]);
ylabel('Canonical value')
xlabel('Canonical variate')
title('Canonical spectrum')


% canonical value spectrum
%---------------------------------------------------------------------------
subplot(2,2,4)
semilogy(pV,'o');
axis square;
ylabel('P{dimensionality > S}')
xlabel S
title('Dimensionality of the result')
grid on


spm_print



%---------------------------------------------------------------------------
figure(2); clf; set(2,'Name','','Pointer','Arrow'); figure(3)
