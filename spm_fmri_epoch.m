% function spm_fmri_epoch
% MAN[C]OVA and CVA of adjusted data
% FORMAT spm_fmri_epoch
%___________________________________________________________________________
%
% spm_fmri_epoch is a multivariate analysis of fMRI time-series that
% characterizes spatio-temporal response patterns in terms of (i) a p-value
% and (ii) the first canonical image.
% One key aspect of this analysis is that each 'obervation' comprises
% a complete 'epoch' of data (i.e. a series of scans time-locked to some
% event - usually a change in condition).  This multivariate observation
% has both spatial and temporal dimensions.  The canonical images that
% ensue from this analysis are spatio-temporal modes that constitute
% a time-series (that has the same length as an epoch).  This spatio-
% temporal mode can be decomposed into a set of spatial modes in the usual
% way.
%
% The data that are subject to MANCOVA are the adjusted data matrices
% XA (in XA.mat) that obtain after removing any confounds using ANCOVA
% on a scan by scan basis (i.e. using spm_fmri_spm[_ui] with only
% covariates that are considered confounds (e.g. time and global activity)
%
% spm_fmri_mancova saves the eigenimage decomposition of the first canonical
% image in CVA.mat.  The spatial modes can be re-displayed using spm_fmri_cva.
%
%__________________________________________________________________________
% %W% %E%


% get adjusted data
%---------------------------------------------------------------------------
figure(2); clf; set(2,'Name','MAN[C]OVA - Canonical analysis'); clear
tmp   = spm_get(1,'.mat','select adjusted data','XA');
CWD   = strrep(tmp,'/XA.mat','');
cd(CWD);

load XA
load SPM
load XYZ

[m n] = size(XA);
if m*n > 3000000;
	XA  = XA(:,[1:2:n]);
	XYZ = XYZ(:,[1:2:n]);
end


% get design matrix
%---------------------------------------------------------------------------
[m n] = size(XA);
l     = spm_input('number of scans per epoch '  ,1);
q     = spm_input('number of epochs (e.g. 10)'  ,2);
I     = spm_input(sprintf('{%0.0f x 1} start scans',  q),3);
h     = spm_input('# of Covariates of interest',4);
for i = 1:h
	d = spm_input(sprintf('{%0.0f} Covariates {of interest}',q),4);
	d = d(1:q);
	H = [H d(:)];

end
b     = spm_input('# of Covariates of no interest',5);
for i = 1:b
	d = spm_input(sprintf('{%0.0f} Covariates {no interest}',q),5);
	d = d(1:q);
	B = [B d(:)];
end

%---------------------------------------------------------------------------
set(2,'Name','','Pointer','Watch')


% redundant partitions of design matrix
%---------------------------------------------------------------------------
I     = I(1:q);
C     = [];
G     = [];


% Reorganize data for multivariate analysis
%---------------------------------------------------------------------------
XA    = XA - ones(m,1)*mean(XA);
XA    = XA'; XA = XA(:);
X     = zeros(q,n*l);
for j = 1:q
	d        = [1:n*l] + (I(j) - 1)*n;
	d(d > m) = d(d > m) - m;
	X(j,:)   = XA(d)';  end

XA    = X; clear X


% SVD DIMENSION REDUCTION
%===========================================================================

% degrees of freedom (specifed by the sizes of the design matrix partitions)
%---------------------------------------------------------------------------
q     = size([H C B G],1);				% total
h     = rank(H);					% condition
c     = rank(C);					% covariates
b     = rank(B);					% block
g     = rank(G);					% confounds
r     = q - rank([H C B G]);				% residuals [error]

% reduce to the eigenimage space
%---------------------------------------------------------------------------
A     = zeros(q);
for i = 1:q
	A(:,i) = XA*(XA(i,:)'); end

[e s] = eig(A);
[i j] = sort(-diag(s));
s     = s(j,j);
e     = e(:,j);

v     = diag(s);
v     = length(v)*v/sum(v);

Q     = 1:(rank(A) - rank([H C B G]));
X     = e(:,Q)*s(Q,Q)^(0.5);
p     = size(X,2);

% MANCOVA
%===========================================================================


% estimate parameters
%---------------------------------------------------------------------------
BETA  = [H C B G]\X;
T     = [H C]*BETA([1:size([H C],2)],:);
R     = X - [H C B G]*BETA;
N     = X;
if ~isempty([B G])
	N     = X - [B G]*([B G]\X); end

% sum of squares
%----------------------------------------------------------------------------
T     = T'*T;
R     = R'*R;
N     = N'*N;

% test for significance of [H C] with Wilks Lambda
%----------------------------------------------------------------------------
W     = det(R)/det(N);
chi   = -(r - ((p - h + 1)/2))*log(W);   	% as Chi-squared(p*h)
pV    = 1 - spm_Xcdf(chi,h*p);		 	% p value


% CANONICAL VARIATES ANALYSIS
% Identify the first canonical image and variate
%===========================================================================
[E v] = eig(T,R);
[i j] = max(diag(v));
v     = v(j,j)*r/h;
E     = real(E(:,j));

S     = X*E;					% canonical variate
E     = e(:,Q)*s(Q,Q)^(-0.5)*E;
U     = XA'*E;					% canonical mode

% reformat canonical vector into a time-series in voxel space
%---------------------------------------------------------------------------
U     = reshape(U,n,l)';


% maximum intensity projections and textual information
%===========================================================================
clear XA; figure(3); spm_clf


axes('Position',[0.05 0.5 0.5 0.4])
spm_mip(sum(U.^2),XYZ,V(1:6));
title('canonical image {2 norm over time}');

% textual information
%---------------------------------------------------------------------------
axes('Position',[0.58 0.6 0.5 0.16],'Visible','off'); y = 1;
text(0,y,'Canonical analysis','FontSize',16,'Fontweight','Bold'); y = y - 0.2;
text(0,y,CWD); y = y - 0.2;
text(0,y,sprintf('Chi-squared = %0.2f',chi)); y = y - 0.2;
text(0,y,sprintf('threshold @ p = 0.05 = %0.2f',spm_invXcdf(0.95,h*p)));
y = y - 0.2; text(0,y,sprintf('p value = %0.6f',pV),'FontSize',16);

% Design matrix
%---------------------------------------------------------------------------
axes('Position',[0.1 0.16 0.3 0.3]);
imagesc(-[H C B G]*inv(diag(max(abs([H C B G])))))
axis image;
title ('Design matrix','Fontweight','Bold')
xlabel 'effect'; ylabel 'epoch'

% canonical variate
%---------------------------------------------------------------------------
axes('Position',[0.5 0.16 0.4 0.3]);
[x y] = bar(S);
fill(x,y,[1 1 1]*0.9);
xlabel 'epoch'
ylabel 'canonical variate'
set(gca,'XLim',([0 length(S)] + 1/2))




% spatiotemporal characterization of canonical vector using SVD
%===========================================================================


% characterize and display canonical images
%---------------------------------------------------------------------------
[CV s CU] = svd(U',0);
s         = diag(s).^2;
s         = s*length(s)/sum(s);

% save variables for subsequnet analysis
%---------------------------------------------------------------------------
df    = h*p;
save CVA CV CU S H C B G X T R N v s pV chi V df XYZ

% eigenimages
%---------------------------------------------------------------------------
figure(3); spm_clf
if V(3) == 1						% 2-dimensional data
	subplot(2,1,1)
	spm_mip(CV(:,1),XYZ,V(1:6));
	title('first eigenimage');
else
	axes('Position',[0.05 0.5 0.45 0.4])
	spm_mip(CV(:,1),XYZ,V(1:6));
	title('first eigenimage {+ve}');

	axes('Position',[0.5 0.5 0.45 0.4])
	spm_mip(-CV(:,1),XYZ,V(1:6));
	title('first eigenimage {-ve}');
end

axes('Position',[0 0.5 1 1]); axis off
text(0.1,0,'Canonical image analysis:','Fontsize',16,'FontWeight','Bold')
text(0.52,0,CWD)


% eigenvalue spectrum
%---------------------------------------------------------------------------
subplot(2,2,3)
d     = zeros(size(s)); d(1) = s(1);
[x y] = bar(s);
fill(x,y,[1 1 1]*0.9); hold on
[x y] = bar(d);
fill(x,y,[1 0 0]); hold off
axis square;
xlabel eigenimage
ylabel eigenvalue
title(sprintf('%0.1f percent of variance',s(1)/length(s)*100))

% component scores
%---------------------------------------------------------------------------
subplot(2,2,4)
[x y] = bar(CU(:,1));
fill(x,y,[1 1 1]*0.8);
axis square
xlabel 'time'
ylabel 'component score'

%---------------------------------------------------------------------------
figure(2); clf; set(2,'Name','','Pointer','Arrow'); figure(3)
