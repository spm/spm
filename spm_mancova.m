% function spm_mancova
% Statistical analysis with the General linear model
% FORMAT spm_mancova
%____________________________________________________________________________
%
% spm_mancova is the mutlivariate equivalent of spm_spm.  It implements
% a ManCova following dimension reduction using SVD (only eigenimages
% with eigenvalues greater than unity are used).  The effects of interest
% are characterized using CVA or Canonical Variates Analysis.
%
% Canonical imges are displayed post hoc using spm_cva
%
%__________________________________________________________________________
% %W% Andrew Holmes %E%

%---------------------------------------------------------------------------
global CWD


% get filenames and other user specified parameters
%============================================================================

% get filenames and create design matrix
%----------------------------------------------------------------------------
set(2,'Name','MAN[C]OVA'); drawnow

Q     = [];				% matrix of filename strings
X     = [];				% data matrix {q x pixels}
H     = [];				% condition partition of design
C     = [];				% covariates of interest
B     = [];				% block partition of design
G     = [];				% covariates of no interest
T     = [];				% row matrix of contrasts
K     = [];				% string matrix for graphical output


% cycle over studies, subjects and conditions
%----------------------------------------------------------------------------
s     = spm_input('number of studies',1);
for i = 1:s
	n      = spm_input(['study ' num2str(i) ': number of subjects'],2);
	k      = spm_input(['number of conditions {scans}'],3);
	q      = n*k;
	if (k == 1)
		d = ['select subject 1 - ' num2str(n) '  {study '...
		      num2str(i) '}'];
 		Q = str2mat(Q,spm_get(n,'.img',d));
		K = [K; ones(n,1)*i [1:n]' ones(n,1)];
	else
		for j  = 1:n
			d = ['select conditions 1 - ' num2str(k) '  {study '...
		     	      num2str(i) ', subject ' num2str(j) '}'];
 			Q = str2mat(Q,spm_get(k,'.img',d));
			K = [K; ones(k,1)*i ones(k,1)*j [1:k]' ];
		end
	end

	% design matrix subpartitions
	%-------------------------------------------------------------------
	h = []; for j = 1:n; h = [h; eye(k)]; end
	b = []; for j = 1:k; b = [b; eye(n)]; end
	b = reshape(b,n,q)';

	% condition block subpartition
	%-------------------------------------------------------------------
	b([1:k],:) = b([1:k],:) - 1;
	b(:,1)     = [];

	% update H and B
	%-------------------------------------------------------------------
	[u v] = size(H);
	if (length(h))
		H(([1:q] + u),([1:k] + v))  = h;
		B(([1:q] + u),[1:size(B,2)]) = zeros(q,size(B,2));
	end
	if (size(b,2))
		B(([1:q] + u),([1:(n - 1)] + size(B,2))) = b;		
	end
end

Q(1,:) = []; P = Q;

% total observations
%----------------------------------------------------------------------------
q     = size(H,1);

% remove redundant parititions of design matrix
%---------------------------------------------------------------------------
h     = size(H,2);
b     = size(B,2) + s;
if (h == q) | (h == 1); H = []; end
if (b == q) | (b == 1); B = []; end



% covariates of interest
%----------------------------------------------------------------------------
c     = spm_input('number of covariates {of interest}',4);
for i = 1:c
	d = spm_input(sprintf('[%0.0f]- covariate %0.0f',q,i),5);
	while length(d) ~= q
	    d = spm_input([num2str(length(d)) '-variate entry, try again'],5);
	end
    C = [C d(:)];
end


% covariates not of interest
%----------------------------------------------------------------------------
g     = spm_input('number of covariates {no interest}',5);
for i = 1:g
	d = spm_input(sprintf('[%0.0f]- covariate %0.0f',q,i),6);
	while length(d) ~= q
	    d = spm_input([num2str(length(d)) '-variate entry, try again'],6);
	end
    G = [G d(:)];
end

% covariates not of interest - linear de-trend
%----------------------------------------------------------------------------
d     = spm_input('linear de-trend {y/n}',6,'s');
if size(B,2) & (d == 'y')
	d     = find(B(:,1) == -1);
	j     = zeros(q,1);
	j(d)  = d - mean(d);
	G     = [G j];
	for i = 1:size(B,2)
		d     = find(B(:,i) == 1);
		j     = zeros(q,1);
		j(d)  = d - mean(d);
		G     = [G j];
	end
elseif  (d == 'y')
	j     = 1:q;
	G     = [G (j(:) - mean(j))];
end

% global normalization with ANCOVA as the default
%----------------------------------------------------------------------------
GLOBAL = ' ';
while (GLOBAL ~= 'A' & GLOBAL ~= 'S' & GLOBAL ~= 'N')
	GLOBAL = spm_input('M[A]NCOVA [S]caling or [N]one ?',7,'s'); end

if GLOBAL == 'A'; GLOBAL = 'ANCOVA';  end
if GLOBAL == 'S'; GLOBAL = 'Scaling'; end
if GLOBAL == 'N'; GLOBAL = 'None';    end

set(get(2,'CurrentObject'),'String',GLOBAL);


% the interactive parts of spm_spm_ui are now finished
%----------------------------------------------------------------------------
set(2,'Name','thankyou','Pointer','Watch')




% get file identifiers and global values
%============================================================================
	
% 50 is the value of GM to which the data are scaled (c.f Fox and Mintun 1988)
%----------------------------------------------------------------------------
GM    = 50;
V     = zeros(12,q);
for i = 1:q; V(:,i) = spm_map(P(i,:));  end

% check for consistency of image size and voxel size
%----------------------------------------------------------------------------
if ~(all(all(~diff(V([1:6],:)'))))
	error('images are not compatible'); end

% get ORIGIN
%----------------------------------------------------------------------------
[DIM VOX SCALE TYPE OFFSET ORIGIN] = spm_hread(P(1,:));
if DIM(3) == 1; ORIGIN = [1 1 1]; end


% adjust scaling coefficients so that Grand mean - <m> = GM
%----------------------------------------------------------------------------
m      = zeros(q,1);
for i  = 1:q
	m(i) = spm_global(V(:,i)); end


V(7,:) = V(7,:)*GM/mean(m);
m      = m*GM/mean(m);

if strcmp(GLOBAL,'ANCOVA');  G      = [G m];         end
if strcmp(GLOBAL,'Scaling'); V(7,:) = GM*V(7,:)./m'; m = GM*m./m;	end




% display analysis parameters
%============================================================================
figure(3); spm_clf; axis off
text(0.50,0.86,'Filename');
text(0.30,1.02,'Canonical analysis','Fontsize',16,'Fontweight','Bold');
text(0.00,0.86,'Study',     'Rotation',90);
text(0.06,0.86,'Block',     'Rotation',90);
text(0.12,0.86,'Condition', 'Rotation',90);
text(0.18,0.86,'Covariates','Rotation',90);
y     = 0.82;
for i = 1:q
	d   = [];
	if length([C G])
		d = [C G];
		d = sprintf('%-8.2f',d(i,:));
	end
	d   = [sprintf('%-6.0f',K(i,:)) '   ' d '      ' P(i,:)];
	text(0,y,d,'FontSize',10)
	y   = y - 0.02;
	if y < 0;
		spm_print
		spm_clf;axis off
		y   = 0.9;
		text(0.16,1.02,['Statistical analysis (continued)'],...
		    'Fontsize',16,'Fontweight','Bold');
	end
end

spm_print


% GET DATA, REMOVE CONFOUNDS
%===========================================================================

% vextors of location
%---------------------------------------------------------------------------
[y x] = meshgrid([1:V(2,1)],[1:V(1,1)]');
z     = ones(size(x));
x     = (x(:) - ORIGIN(1))*V(4,1);
y     = (y(:) - ORIGIN(2))*V(5,1);
z     = z(:);


% variables saved (at voxels satisfying P{F > x} < 0.05)
%---------------------------------------------------------------------------
eval(['cd ' CWD]);
delete M.mat
delete XYZ.mat
delete CVA.mat

% mean correct the covariates and factors (mean = 0)
%---------------------------------------------------------------------------
if ~isempty(C); C = C - ones(q,1)*mean(C); end
if ~isempty(G); G = G - ones(q,1)*mean(G); end


% get and adjust data at a univariate level
%---------------------------------------------------------------------------
GX    = m(:)*ones(1,V(1,1)*V(2,1));			% global activities

for i = 1:V(3)

  X     = zeros(q,V(1,1)*V(2,1));
  for j = 1:q
      d = spm_slice_vol(V(:,j),spm_matrix([0 0 i]),[V(1,1) V(2,1)],0);
      X(j,:) = d(:)';
  end

  % select voxels for analysis (gray matter > 36 a.u. GM = 50 a.u.)
  %-------------------------------------------------------------------------
  Q     = find(all(X > 0.8*GX) & ~all(~diff(X)));

  if length(Q) > size(X,2)/8

  	XYZ   = [x y z*(i - ORIGIN(3))*V(6,1)]';
 	XYZ   = XYZ(:,Q);
 	X     = X(:,Q);
 	X     = X - ones(q,1)*mean(X);
	X     = X - [B G]*pinv([B G])*X;

  	spm_append('M',X);
 	spm_append('XYZ',XYZ);
  end

end


% SVD DIMENSION REDUCTION
%===========================================================================


% load and normalize data
%---------------------------------------------------------------------------
load M
A     = M*M';

% degrees of freedom (specifed by the sizes of the design matrix partitions)
%---------------------------------------------------------------------------
q     = size(A,1);					% total
h     = rank(H) - 1;					% conditions
c     = rank(C);					% covariates
r     = q - rank([H C B G]);				% residuals [error]


% SVD
%---------------------------------------------------------------------------
[e s] = eig(A);
[i j] = sort(-diag(s));
s     = s(j,j);
e     = e(:,j);

v     = diag(s);
v     = length(v)*v/sum(v);

Q     = 1:min([ (rank(A) - rank([H C B G])) max(find(v > 1))]);
X     = e(:,Q)*s(Q,Q)^(0.5);
p     = size(X,2);


% MANCOVA
%===========================================================================


% estimate parameters and sum of squares due to error
%---------------------------------------------------------------------------
BETA  = pinv([H C B G])*X;
T     = [H C]*BETA(1:size([H C],2),:);
R     = X - [H C B G]*BETA;
N     = X - [B G]*pinv([B G])*X; XA = N;

T     = T'*T;
R     = R'*R;
N     = N'*N;

% test for significance of [H C] with Wilks Lambda
%----------------------------------------------------------------------------
W     = det(R)/det(N);
chi   = -(r - ((p - (h + c) + 1)/2))*log(W);	% distributed Chi-squared(p*h)
pV    = 1 - spm_Xcdf(chi,(h + c)*p);		% p value
UF    = spm_invFcdf((1 - 0.05),[(h + c),r]);


% CANONICAL VARIATES ANALYSIS
%===========================================================================
[E v] = eig(inv(R)*T);
[v j] = sort(-diag(real(v)));
v     = -v*r/(h + c);
E     = E(:,j);
E     = E(:,(v > UF));

d     = e(:,Q)*s(Q,Q)^(-0.5)*E;
CV    = M'*d;
S     = BETA(1:size([H C],2),:)*E;
if size(S,1) == 1
	S    = XA*E; end

% save variables for subsequent analysis
%---------------------------------------------------------------------------
V     = [V(1:6,1); ORIGIN(:)];
df    = r;

save CVA CV S H C B G X T R N v pV chi V df



% textual information
%---------------------------------------------------------------------------
figure(3); spm_clf

axes('Position',[0.3 0.54 0.36 0.16],'Visible','off'); y = 1;
text(0,y,'Multivariate analysis','FontSize',16,'Fontweight','Bold');
y = y - 0.2;
text(0,y,pwd,'Fontweight','Bold'); y = y - 0.2;
text(0,y,sprintf('Chi-squared = %0.2f',chi)); y = y - 0.2;
text(0,y,sprintf('threshold @ p = 0.05 = %0.2f',spm_invXcdf(0.95,(h + c)*p)));
y = y - 0.2; text(0,y,sprintf('p value = %0.6f',pV),'FontSize',16);


% Design matrix
%---------------------------------------------------------------------------
axes('Position',[0.1 0.16 0.26 0.3]);
imagesc(-[H C B G]*inv(diag(max(abs([H C B G])))))
title ('Design matrix','Fontweight','Bold')
xlabel 'effect'; ylabel 'scan'

% canonical values
%---------------------------------------------------------------------------
axes('Position',[0.5 0.16 0.4 0.3]);
[x y] = bar(v.*(v > 1));
fill(x,y,[1 1 1]*0.9); hold on
[x y] = bar(v(v > UF));
fill(x,y,[1 0 0]); hold off
line([1 length(v)],[1 1]*UF,'Color',[1 0 0])
xlabel('Canonical variate')
ylabel('Canonical value')
title ('Canonical spectrum','Fontweight','Bold')


spm_print



% display and print canonical images
%---------------------------------------------------------------------------
load XYZ
d     = flipud(find(v > UF));
for j = 1:length(d)
	i = d(j);
	figure(3); spm_clf
	if V(3,1) == 1					% 2-dimensional data

		subplot(2,1,1)
		spm_mip(CV(:,i),XYZ,V(1:6));
		title(sprintf('canonical image %0.0f',i));
	else
		axes('Position',[0.05 0.5 0.45 0.4])
		spm_mip(CV(:,i),XYZ,V(1:6));
		title(sprintf('canonical image %0.0f {+ve}',i));

		axes('Position',[0.55 0.5 0.45 0.4])
		spm_mip(-CV(:,i),XYZ,V(1:6));
		title(sprintf('canonical image %0.0f {-ve}',i));
	end

	% component scores
	%-------------------------------------------------------------------
	subplot(2,1,2);
	[x y] = bar(S(:,i));
	fill(x,y,[1 1 1]*0.9);
	xlabel 'condition'
	ylabel 'canonical variate'
	title(sprintf('canonical value = %0.0f ',v(i)),'FontSize',16)
	axis square
	set(gca,'Xlim',[0 (size(S,1) + 1)])
	spm_print
end

% clear figure 2
%---------------------------------------------------------------------------
figure(2); spm_clf
set(2,'Name',' ','Pointer','Arrow');
