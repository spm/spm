function spm_segment(PF,PG,opts)
% Segment an MR image into Gray, White, CSF & other.
%
% --- The Prompts Explained ---
%
% 'select MRI(s) for subject '
% If more than one volume is specified (eg T1 & T2), then they must be
% in register (same position, size, voxel dims etc..).
%
% 'select Template(s) '
% If the images have been spatially normalised, then there is no need to
% select any images. Otherwise, select one or more template images which
% will be used for affine normalisation of the images. Note that the
% affine transform is only determined from the first image specified
% for segmentation. 
%
% 'Attempt to correct intensity inhomogenaeties?'
% This uses a Bayesian framework (again) to model intensity
% inhomogeneities in the image(s).  The variance associated with each
% tissue class is assumed to be multiplicative (with the
% inhomogeneities).  The low frequency intensity variability is
% modelled by a linear combination of three dimensional DCT basis
% functions (again), using a fast algorithm (again) to generate the
% curvature matrix.  The regularization is based upon `bending energy'
% (like Bookstein's thin plate splines).  A small amount of
% regularization is used when correcting for `Lots of inhomogenaety',
% whereas more regularization is used for `A little inhomogenaety'.
%
%_______________________________________________________________________
%
%                      The algorithm is three step:
%
% 1) Determine the affine transform which best matches the image with a
%    template image. If the name of more than one image is passed, then
%    the first image is used in this step. This step is not performed if
%    no template images are specified.
%    Note: the templates and probability images are Right-is-Right
%    orientation. The affine normalisation assumes that the input images
%    are the same - although the algorithm works well enough with input
%    images in the Left-is-Right orientation.
%
% 2) Perform Cluster Analysis with a modified Mixture Model and a-priori
%    information about the likelihoods of each voxel being one of a
%    number of different tissue types. If more than one image is passed,
%    then they they are all assumed to be in register, and the voxel
%    values are fitted to multi-normal distributions.
%
% 3) Write the segmented image. The names of these images have
%    "_seg1", "_seg2" & "_seg3" appended to the name of the
%    first image passed.
%
%_______________________________________________________________________
%
% Ref:
% Hartigan, J. A. 1975. Clustering Algorithms, pp. 113-129.
% John Wiley & Sons, Inc., New York.
%
%_______________________________________________________________________
%
% The template image, and a-priori likelihood images are modified
% versions of those kindly supplied by Alan Evans, MNI, Canada
% (ICBM, NIH P-20 project, Principal Investigator John Mazziotta).
%_______________________________________________________________________
% %W% (c) John Ashburner %E%

% Programmers notes
%
% FORMAT spm_segment(PF,PG,opts)
% PF   - name(s) of image(s) to segment (must have same dimensions).
% PG   - name(s) of template image(s) for realignment.
%      - or a 4x4 transformation matrix which maps from the image to
%        the set of templates.
% opts - options string.
%        - 't' - write images called *_seg_tmp* rather than *_seg*
%                (that are smoothed with an 8mm Gaussian).
%        - 'f' - fix number of voxels in each cluster
%        - 'c' - attempt to correct small intensity inhomogeneities
%        - 'C' - attempt to correct large intensity inhomogeneities

debug = 0;

if nargin<3
	% set to ' ' rather than '' to get rid of annoying warnings
	%opts = 'f';
	opts = ' ';
	if nargin<2
		PG = '';
	end
end

global SWD
DIR1   = [SWD '/coreg/'];
DIR2   = [SWD '/apriori/'];

if (nargin==0)
	spm_figure('Clear','Interactive');
	set(spm_figure('FindWin','Interactive'),'Name','Segmentation')
	spm_help('!ContextHelp','spm_segment.m');

	n     = spm_input('number of subjects',1,'e',1);
	if (n < 1)
		spm_figure('Clear','Interactive');
		return;
	end


	for i = 1:n
		PF = spm_get(Inf,'.img',...
			['Select MRI(s) for subject ' num2str(i)]);
		eval(['PF' num2str(i) ' = PF;']);
	end

	PG = '';
	tmp = spm_input('Are they spatially normalised?', 1, 'y/n');

	if (tmp == 'n')
		% Get template
		%-----------------------------------------------------------------------
		templates = str2mat([DIR1 'T1.img'], [DIR1 'T2.img'], [DIR1 'EPI.img']);

		% Get modality of target
		respt = spm_input('Modality of first image?','+1','m',...
			'modality - T1 MRI|modality - T2 MRI|modality - EPI MR|--other--',...
			[1 2 3 0],1);
		if (respt > 0)
			PG = deblank(templates(respt,:));
		else
			ok = 0;
			while (~ok)
				PG = spm_get(Inf,'.img',['Select Template(s) for affine matching'],...
					'', DIR1);
				if (size(PG,1)>0)
					dims = zeros(size(PG,1),9);
					for i=1:size(PG,1)
						[dim vox dummy dummy dummy origin dummy]...
							 = spm_hread(deblank(PG(i,:)));
						dims(i,:) = [dim vox origin];
					end
					if size(dims,1) == 1 | ~any(diff(dims))
						ok = 1;
					end
				else
					ok = 0;
				end
			end
		end
	end

	if (1),
	% This stuff isn't ready for general consumption yet
	tmp = spm_input('Attempt to correct intensity inhomogenaeties?','+1','m',...
		['No inhomogenaety correction|' ...
		 'A little inhomogenaety correction|' ...
		 'Lots of inhomogenaety correction'],...
		[' ' 'c' 'C'],1);
	opts = [tmp opts];
	end

	for i = 1:n
		eval(['PF = PF' num2str(i) ';']);
		set(spm_figure('FindWin','Interactive'),'Name',...
			['Segmenting ' num2str(i) '..'],'Pointer','Watch');
		drawnow;
		if (size(PF,1)~=0) spm_segment(PF,PG,opts); end
	end
	spm_figure('Clear','Interactive');
	return;
end

%_______________________________________________________________________
%_______________________________________________________________________

%- A-Priori likelihood images.
PB    = str2mat([DIR2 'gray.img'],[DIR2 'white.img'],[DIR2 'csf.img']);

niter     = 64;      % Maximum number of iterations of Mixture Model
nc        = [1 1 1 3]; % Number of clusters for each probability image


% Determine matrix MM that will transform the image to Talairach space.
%_______________________________________________________________________
MF=zeros(4,4,size(PF,1));
for i=1:size(PF,1)
	MF(:,:,i) = spm_get_space(PF(i,:));
end
if ~isempty(PG) & isstr(PG)
	% Affine normalisation
	%-----------------------------------------------------------------------
	MG = spm_get_space(PG(1,:));
	spm_smooth(PF(1,:),'./spm_seg_tmp.img',8);

	% perform affine normalisation at different sampling frequencies
	% with increasing numbers of parameters.
	disp('Determining Affine Mapping');
	global sptl_Ornt
	if prod(size(sptl_Ornt)) == 12
		prms = [sptl_Ornt ones(1,size(PG,1))]';
	else
		prms = [0 0 0  0 0 0  1 1 1  0 0 0 ones(1,size(PG,1))]';
	end
	spm_chi2_plot('Init','Affine Registration','Convergence');
	prms = spm_affsub3('affine3', PG, './spm_seg_tmp.img', 1, 8, prms);
	prms = spm_affsub3('affine3', PG, './spm_seg_tmp.img', 1, 6, prms);
	spm_chi2_plot('Clear');
	spm_unlink ./spm_seg_tmp.img ./spm_seg_tmp.hdr ./spm_seg_tmp.mat
	MM = spm_matrix(prms);

elseif all(size(PG) == [4 4])
	% Assume that second argument is a matrix that will do the job
	%-----------------------------------------------------------------------
	MM = PG;
else
	% Assume that image is normalized
	%-----------------------------------------------------------------------
	MM = spm_matrix([0 0 0 0 0 0 1 1 1 0 0 0]');
end


for i=1:size(PF,1)
	VF(:,i) = spm_map(PF(i,:));
end

Mat = zeros(4,4,size(PB,1));
for j=1:size(PB,1)
	VB(:,j)    = spm_map(PB(j,:));
	Mat(:,:,j) = spm_get_space(PB(j,:));
end

m  = size(PF,1);
%-----------------------------------------------------------------------



% Voxels to sample during the cluster analysis
%-----------------------------------------------------------------------

% A bounding box for the brain in Talairach space.
%bb1 = [ [-78 78]' [-112 76]' [-50 85]'];
bb1 = [ [-88 88]' [-122 86]' [-60 95]'];

% A mapping from a unit radius sphere to a hyper-ellipse
% that is just enclosed by the bounding box in Talairach
% space.
M0 = [diag(diff(bb1)/2) mean(bb1)';[0 0 0 1]];

% The mapping from voxels to Talairach space is
% MM*MF(:,:,1), so the ellipse in the space
% of the first image becomes:
M0 = inv(MM*MF(:,:,1))*M0;

% So to work out the bounding box in the space of the
% image that just encloses the hyper-ellipse.
tmp = M0(1:3,1:3);
tmp = diag(tmp*tmp'/diag(sqrt(diag(tmp*tmp'))));
bb  = round([M0(1:3,4)-tmp M0(1:3,4)+tmp])';

% Want to sample about every 4mm
tmp  = inv(MF(:,:,1));
tmp  = tmp(1:3,1:3);
samp = round(max(abs(tmp*[4 4 4]'), [1 1 1]'));
%-----------------------------------------------------------------------



reg = 0;
if any(opts == 'c') | any(opts == 'C')
	reg = 1e6;
	if any(opts == 'C'), reg = 1e5; end

	% Stuff for intensity modulation
	%-----------------------------------------------------------------------
	% Set up basis functions
	nbas = max(round((VF(1:3,1).*VF(4:6,1))/32)',[1 1 1]);
	B1=spm_dctmtx(VF(1,1),nbas(1),bb(1,1):samp(1):bb(2,1));
	B2=spm_dctmtx(VF(2,1),nbas(2),bb(1,2):samp(2):bb(2,2));
	B3=spm_dctmtx(VF(3,1),nbas(3),bb(1,3):samp(3):bb(2,3));

	d = [size(B1,1) size(B2,1) size(B3,1)];

	% Set up a priori covariance matrix (based upon bending energy)
	kx=(pi*((1:nbas(1))'-1)/VF(1,1)).^2; ox=ones(nbas(1),1);
	ky=(pi*((1:nbas(2))'-1)/VF(2,1)).^2; oy=ones(nbas(2),1);
	kz=(pi*((1:nbas(3))'-1)/VF(3,1)).^2; oz=ones(nbas(3),1);
	IC0 = (kron(kz.*kz,kron(oy,ox)) + kron(oz,kron(ky.*ky,ox)) + kron(oz,kron(oy,kx.*kx)) ...
	+ 2*kron(kz,kron(ky,ox)) + 2*kron(kz,kron(oy,kx)) + 2*kron(oz,kron(ky,kx)) )*reg;
	% Assume a tiny variance for the DC coefficient.
	IC0(1) = max(IC0)*1000;
	IC0 = diag(IC0);

	% Mode of the a priori distribution
	X0    = zeros(prod(nbas),1);
	X0(1) = sqrt(prod(VF(1:3,1)));

	% Initial estimate for intensity modulation field
	T =zeros(nbas(1),nbas(2),nbas(3),m);
	T(1,1,1,:)=sqrt(prod(VF(1:3,1)));
	%-----------------------------------------------------------------------
else
	T = sqrt(prod(VF(1:3,1)));
	nbas = [1 1 1];
end


lkp=[];
for i=1:size(nc,2)
	lkp = [lkp ones(1,nc(i))*i];
end

n  = size(lkp,2);
nb = size(VB,2);

sumbp = zeros(1,n);
osumpr = -Inf;

cv = zeros(m,m,n);	% (Co)variances
mn = zeros(m,n);	% Means
mg = zeros(1,n);	% Number of voxels/cluster

z = bb(1,3):samp(3):bb(2,3);
d = [length(bb(1,1):samp(1):bb(2,1)) length(bb(1,2):samp(2):bb(2,2)) length(bb(1,3):samp(2):bb(2,3))];

spm_chi2_plot('Init','Segmenting','Log-likelihood','Iteration #');
disp('Segmenting'); 
for iter = 1:niter

	% Initialize variables that are modified during the loop over each plane
	%-----------------------------------------------------------------------
	sumpr= 0;
	mom0 = zeros(1,n)+eps;
	mom1 = zeros(m,n);
	mom2 = zeros(m,m,n)+eps;

	if reg~=0
		Alpha = zeros(prod(nbas),prod(nbas),m);
		Beta  = zeros(prod(nbas),m);
	end

	%-----------------------------------------------------------------------
	for pp = 1:length(z),
		clear dat bp pr


		B = spm_matrix([-bb(1,1)/samp(1)+1 -bb(1,2)/samp(2)+1 -z(pp)...
			0 0 0 1/samp(1) 1/samp(2) 1]);

		% Ignore voxels of value zero - since we don't know if it is because they
		% are truly zero - or if it is because they are zeros due to some kind of
		% image editing or from outside the original FOV
		msk = zeros(d(1)*d(2),1);


		for i=1:m
			M = inv(B);
			tmp = spm_slice_vol(VF(:,i),M, d(1:2),1);
			msk = msk | tmp(:)==0;
			if reg~=0
				% Non-uniformity correct.
				t = reshape(reshape(T(:,:,:,i),...
					nbas(1)*nbas(2),nbas(3))*B3(pp,:)', nbas(1), nbas(2));
				rawdat(:,:,i) = tmp;
				tmp = tmp.*(B1*t*B2');
				dat(:,i) = tmp(:);
			else
				dat(:,i)=tmp;
			end
		end
		msk = find(~msk);

		% If there are at least some voxels to work with..
		if length(msk)>0,

			% A priori probability data for eg WM GM CSF scalp etc..
			bp = zeros(d(1)*d(2),nb+1);
			for j=1:nb
				M = inv(B*inv(MM*MF(:,:,1))*Mat(:,:,j));
				tmp = spm_slice_vol(VB(:,j),M, d(1:2),1)/nc(j);
				bp(:,j) = tmp(:);
			end
			% Other tissue
			bp(:,nb+1) = abs(ones(d(1)*d(2),1) - bp(:,1:nb)*nc(1:nb)')/nc(nb+1);

			pr = zeros(d(1)*d(2),n);
			if (iter==1)
				% Initial probability estimates based upon
				% a priori knowledge
				%-----------------------------------------------------------------------
				for i=1:n
					pr(msk,i) = bp(msk,lkp(i));
					sumbp(i) = sumbp(i) + sum(bp(msk,lkp(i)));
				end
			else
				% Compute PDFs for each cluster
				%-----------------------------------------------------------------------
				for i=1:n
					amp = 1/sqrt((2*pi)^m * det(cv(:,:,i)));
					dst = (dat(msk,:)-ones(size(msk,1),1)*mn(:,i)')/sqrtm(cv(:,:,i));
					dst = sum(dst.*dst,2);
					pr(msk,i)=amp*exp(-0.5*dst).*(bp(msk,lkp(i))*(mg(1,i)/sumbp(i)));
				end
				%-----------------------------------------------------------------------
			end

			% Compute log likelihood, and normalize likelihoods to sum to unity
			%-----------------------------------------------------------------------
			sp = sum(pr(msk,:),2);
			sumpr = sumpr + sum(log(sp));
			msk2 = find(~sp);
			sp(msk2) = 1;
			for i=1:n
				pr(msk,i) = pr(msk,i)./sp;
			end

			%-----------------------------------------------------------------------


			% Compute new n, mean & var for each cluster - step 1
			%-----------------------------------------------------------------------
			for i=1:n
				mom0(1,i)   = mom0(1,i)   + sum(pr(:,i));
				mom1(:,i)   = mom1(:,i)   + sum((pr(:,i)*ones(1,m)).*dat)';
				mom2(:,:,i) = mom2(:,:,i) + ((pr(:,i)*ones(1,m)).*dat)'*dat;
			end

			if reg~=0 & iter > 1
				% Build up A'*A and A'*b to solve for intensity modulations
				%-----------------------------------------------------------------------
				pr  = reshape(pr ,d(1),d(2),n);
				for j=1:m,
					for i=1:n
						wt = pr(:,:,i)*(cv(j,j,i).^(-0.5));
						if i==1,
							[alpha,beta]=spm_kronutil(wt.*rawdat(:,:,j),wt*mn(j,i),B1,B2);
						else,
							[alph ,bet ]=spm_kronutil(wt.*rawdat(:,:,j),wt*mn(j,i),B1,B2);
							alpha = alpha + alph;
							beta  = beta  + bet;
						end;
					end
					Alpha(:,:,j) = Alpha(:,:,j) + kron(B3(pp,:)'*B3(pp,:),alpha);
					Beta(:,j)    = Beta(:,j)    + kron(B3(pp,:)', beta);
				end
				pr  = reshape(pr ,d(1)*d(2),n);
			end

		end
	end

	% Solve for intensity modulations
	%-----------------------------------------------------------------------
	if reg~=0 & iter>1
		for i=1:m
			x = T(:,:,:,i);x=x(:);
			x = (Alpha(:,:,i) + IC0)\(IC0*X0 + Beta(:,i));
			T(:,:,:,i) = reshape(x,nbas);
		end
	end
	%-----------------------------------------------------------------------

	if iter>2, spm_chi2_plot('Set',sumpr); end


	% Compute new n, mean & var for each cluster - step 2
	%-----------------------------------------------------------------------
	for i=1:n
		mg(1,i) = mom0(1,i);
		if any(opts == 'f') & i<=nb , mg(1,i) = sumbp(i); end;
		mn(:,i) = mom1(:,i)/mom0(1,i);

		tmp = (mom0(1,i).*mn(:,i))*mn(:,i)';
		tmp = tmp-eye(size(tmp))*eps*1000;
		cv(:,:,i) = (mom2(:,:,i) - tmp)/mom0(1,i);
	end
	%-----------------------------------------------------------------------


	if iter==1
		% Split the clusters
		%-----------------------------------------------------------------------
		nn = 0;
		for j=1:length(nc)
			for i=2:nc(j)
				cv(:,:,nn+i) = cv(:,:,nn+i)*0.8^(1-i);
				mn(:,nn+i) = mn(:,nn+i)*0.8^(1-i);
			end
			nn = nn + nc(j);
		end

		% Background Cluster.
		%    Strictly speaking, since voxels contain absolute values,
		%    the distributions should be modified accordingly. However
		%    in practice, the only allowance made for this is for the
		%    distribution of a background cluster. The mean of this
		%    cluster is assumed to be zero, and in order for the model
		%    to fit properly, the number of voxels contained in this 
		%    cluster is doubled.
		%-----------------------------------------------------------------------
		mn(:,n)   = zeros(m,1);
		mg(1,n)   = 2*mg(1,n);
		cv(:,:,n)   = (mom2(:,:,n))/mom0(1,n);
	end

	% Stopping criterion
	%-----------------------------------------------------------------------
	if iter == 2
		sumpr2 = sumpr;
	elseif iter > 3
		if (sumpr-osumpr)/(sumpr-sumpr2) < 0.0003
			break;
		end
	end
	osumpr = sumpr;
end
spm_chi2_plot('Clear');

%save segmentation_results.mat T mg mn cv MM


%-----------------------------------------------------------------------

% Create headers, open files etc...
%-----------------------------------------------------------------------
dm     = VF(1:3,1);
vx     = VF(4:6,1);
orgn   = MF(:,:,1)\[0 0 0 1]';
orgn   = round(orgn(1:3)');
planes = 1:VF(3,1);

k = prod(dm(1:2));

if any(opts == 't')
	app  = '_seg_tmp';
	nimg = 2;
else
	app  = '_seg';
	nimg = nb+1;
end


B1=spm_dctmtx(VF(1,1),nbas(1));
B2=spm_dctmtx(VF(2,1),nbas(2));
B3=spm_dctmtx(VF(3,1),nbas(3));

for j=1:nimg
	iname = [spm_str_manip(PF(1,:),'rd') app num2str(j)];
	fp(j) = fopen([iname '.img'],'w');
	if (fp(j) == -1)
		open_error_message([iname '.img']);
		error(['Failed to open ' iname '.img']);
	end
	spm_hwrite([iname '.img'],dm,vx,1/255,2,0,orgn,...
		'Segmented image');
	spm_get_space(iname,MF(:,:,1));
end


% Write the segmented images.
%-----------------------------------------------------------------------
spm_progress_bar('Init',dm(3),'Writing Segmented','planes completed');
disp('Writing Segmented');
clear pr bp dat

for pp=1:size(planes,2)
	p = planes(pp);
	B  = spm_matrix([0 0 -p]);
	M2 = B*inv(MM*MF(:,:,1));

	for i=1:m
		% The image data
		M1 = inv(B*(MF(:,:,1)\MF(:,:,i)));
		tmp = spm_slice_vol(VF(:,i), M1, dm(1:2), 1);
		if reg~=0
			% Apply non-uniformity correction
			t = reshape(reshape(T(:,:,:,i),...
				nbas(1)*nbas(2),nbas(3))*B3(pp,:)', nbas(1), nbas(2));
			t = B1*t*B2';
			dat(:,i) = tmp(:).*t(:);
		else
			dat(:,i) = tmp(:);
		end
	end

	bp = zeros(size(dat,1),nb+1);
	for j=1:nb
		M = inv(M2*Mat(:,:,j));
		tmp = spm_slice_vol(VB(:,j), M, dm(1:2),1);
		bp(:,j) = tmp(:)/nc(j);
	end

	bp(:,nb+1) = abs(ones(k,1) - bp(:,1:nb)*nc(1:nb)')/nc(nb+1);

	for i=1:n
		amp = 1/sqrt((2*pi)^m * det(cv(:,:,i)));
		dst = (dat-ones(k,1)*mn(:,i)')/sqrtm(cv(:,:,i));
		dst = sum(dst.*dst,2);
		pr(:,i) = amp*exp(-0.5*dst).*(bp(:,lkp(i))*(mg(1,i)/sumbp(i)));
	end

	sp = (sum(pr,2)+eps)/255;

	for j=1:nimg
		tmp = find(lkp(1:(length(lkp)-1))==j);
		if length(tmp) == 1
			dat = pr(:,tmp);
		else
			dat = sum(pr(:,tmp)')';
		end
		dat = round(dat./sp);
		if fwrite(fp(j),dat,'uchar') ~= prod(size(dat))
			write_error_message([spm_str_manip(PF(1,:),'rd') app num2str(j)]);
			error(['Failed to write ' [spm_str_manip(PF(1,:),'rd') app num2str(j)]]);
		end
	end

	spm_progress_bar('Set',pp);
end
spm_progress_bar('Clear');

for j=1:nimg
	fclose(fp(j));
end

for v=[VF VB]
	spm_unmap(v);
end


% Do the graphics
%=======================================================================
spm_figure('Clear','Graphics');
fg=spm_figure('FindWin','Graphics');

% Show some text
%-----------------------------------------------------------------------
ax=axes('Position',[0.05 0.8 0.9 0.2],'Visible','off','Parent',fg);
text(0.5,0.80, 'Segmentation','FontSize',16,'FontWeight','Bold',...
	'HorizontalAlignment','center','Parent',ax);

text(0,0.65, ['Image:  ' spm_str_manip(PF(1,:),'k60d')],...
	'FontSize',14,'FontWeight','Bold','Parent',ax);

text(0,0.40, 'Means:','FontSize',12,'FontWeight','Bold','Parent',ax);
text(0,0.30, 'Std devs:' ,'FontSize',12,'FontWeight','Bold','Parent',ax);
text(0,0.20, 'N vox:','FontSize',12,'FontWeight','Bold','Parent',ax);
for j=1:nb
	text((j+0.5)/(nb+2),0.40, num2str(mn(1,j)),...
		'FontSize',12,'FontWeight','Bold',...
		'HorizontalAlignment','center','Parent',ax);
	text((j+0.5)/(nb+2),0.30, num2str(sqrt(cv(1,j))),...
		'FontSize',12,'FontWeight','Bold',...
		'HorizontalAlignment','center','Parent',ax);
	text((j+0.5)/(nb+2),0.20, num2str(mg(1,j)/sum(mg(1,:))),...
		'FontSize',12,'FontWeight','Bold',...
		'HorizontalAlignment','center','Parent',ax);
end
if m > 1
	text(0,0.10,'Parent',ax, ...
	'Note: only means and variances for the first image are shown','FontSize',12);
end

% and display a few images.
%-----------------------------------------------------------------------
V = spm_map(deblank(PF(1,:)));
for j=1:nimg
	iname = [spm_str_manip(PF(1,:),'rd') app num2str(j) '.img'];
	VS(:,j) = spm_map(iname);
end
M1 = spm_get_space(iname);
M2 = spm_get_space(PF(1,:));
for i=1:5
	M = spm_matrix([0 0 i*V(3)/6]);
	img = spm_slice_vol(V,M,V(1:2),0);
	img(1,1) = eps;
	ax=axes('Position',[0.05 0.75*(1-i/5)+0.05 0.9/(nb+2) 0.75/5],'Visible','off','Parent',fg);
	imagesc(rot90(img), 'Parent', ax);
	set(ax,'Visible','off','DataAspectRatio',[1 1 1]);

	for j=1:nimg
		img = spm_slice_vol(VS(:,j),M2\M1*M,V(1:2),0);
		ax=axes('Position',...
			[0.05+j*0.9/(nb+2) 0.75*(1-i/5)+0.05 0.9/(nb+2) 0.75/5],...
			'Visible','off','Parent',fg);
		image(rot90(img*64), 'Parent', ax);
		set(ax,'Visible','off','DataAspectRatio',[1 1 1]);
	end
end
spm_unmap(V);
for i=1:size(VS,2)
	spm_unmap(VS(:,i));
end
clear VS;

spm_print;
drawnow;


if any(opts == 't')
	fprintf('Smoothing.\n');
	for i=1:nimg
		iname1 = [spm_str_manip(PF(1,:),'rd') app num2str(i)];
		iname2 = [spm_str_manip(PF(1,:),'rd') '_sseg_tmp' num2str(i)];
		spm_smooth([iname1 '.img'],[iname2 '.img'],8);
		spm_unlink([iname1 '.img'], [iname1 '.hdr'], [iname1 '.mat']);
	end
end
return;




function open_error_message(q)
f=spm_figure('findwin','Graphics'); 
if ~isempty(f), 
	figure(f); 
	spm_figure('Clear','Graphics'); 
	spm_figure('Clear','Interactive'); 
	ax=axes('Visible','off','Parent',f); 
	text(0,0.60,'Error opening:', 'FontSize', 25, 'Interpreter', 'none'); 
	text(0,0.55,spm_str_manip(q,'k40d'), 'FontSize', 25, 'Interpreter', 'none'); 
	text(0,0.40,'  Please check that you have write permission.', 'FontSize', 16, 'Interpreter', 'none'); 
end
return;

function write_error_message(q)
f=spm_figure('findwin','Graphics'); 
if ~isempty(f), 
	figure(f); 
	spm_figure('Clear','Graphics'); 
	spm_figure('Clear','Interactive'); 
	ax=axes('Visible','off','Parent',f); 
	text(0,0.60,'Error opening:', 'FontSize', 25, 'Interpreter', 'none'); 
	text(0,0.55,spm_str_manip(q,'k40d'), 'FontSize', 25, 'Interpreter', 'none'); 
	text(0,0.40,'  Please check that you have write permission.', 'FontSize', 16, 'Interpreter', 'none'); 
end
return;

