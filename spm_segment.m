% Segment an MR image into Gray, White, CSF & other.
%
% FORMAT spm_segment
%
% 'select MRI(s) for subject '
% If more than one volume is specified (eg T1 & T2), then they must be
% in register (same position, size, voxel dims etc..).
%
% 'select Template(s) '
% If the images have been spatially normalised, then there is no need to
% select any images. Otherwise, select one or more template images which
% will be used for affine normalisation of the images.
%
% 'Generate:'
%	'Image Segments'
%	Produces segments of gray matter, white matter, csf, and 'other'
%	which have the same image dimensions as the original volume(s).
%
%	'Pseudo Pet'
%	Produces a 'PETlike' image for registration purposes.
%
%	'Reduced Image Segments'
%	Produces reduced size image segments (2mm voxels) which can
%	(possibly) be used for registration purposes. These images can
%	be mapped to the original via their '.mat' files.
%
%	'Reduced Pseudo Pet'
%	Rapidly produces a reduced size 'PETlike' image for quick
%	registration purposes. These images can be mapped to the
%	original via their '.mat' files.
%_______________________________________________________________________
%
% FORMAT spm_segment(PF,PG,opts)
% PF   - name(s) of image(s) to segment (must have same dimensions).
% PG   - name(s) of template image(s) for realignment.
% opts - options string.
%        - 'p' - produce a 'pseudo-pet' image.
%        - 'n' - dont produce the _seg images.
%        - 'r' - write reduced size image (2mm voxels).
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
%    "_seg1", "_seg2", "_seg3" & "_seg4" appended to the name of the
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
%

% %W% (c) John Ashburner %E%

function spm_segment(PF,PG,opts)

debug = 0;

if nargin<3
	opts = '';
end

global SWD
DIR   = [SWD '/mni/'];

if (nargin==0)
	n     = spm_input('number of subjects',1);
	if (n < 1) return; end

	for i = 1:n
		PF = spm_get(Inf,'.img',...
			['select MRI(s) for subject ' num2str(i)]);
		eval(['PF' num2str(i) ' = PF;']);
	end

	% Get template(s)
	ok = 0;
	PG = '';
	while (~ok)
		PG = spm_get(Inf,'.img',['select Template(s) '],...
			'', DIR);
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
			ok = 1; % assume already normalised.
		end
	end
	options = ['   '; 'np '; 'r  '; 'npr'];
	opts = options(spm_input('Generate:',2,'m',...
		'Image Segments|Pseudo Pet|Reduced Image Segments|Reduced Pseudo Pet',...
		[1 2 3 4]),:);

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

%- A-Priori likelihood images
PB    = str2mat([DIR 'pgray.img'],[DIR 'pwhite.img'],[DIR 'pcsf.img']);

niter     = 24;
nc        = [1 1 1 3]; % Number of clusters for each probability image
petscales = [1 0.3 0.1 0.1 0.1 0]';

% Starting parameters for affine normalisation.
% - assume input image uses neurological convention
prms = [0 0 0 0 0 0 1 1 1 0 0 0];

MF = spm_get_space(PF(1,:));
if (~isempty(PG))
	% Affine normalisation
	MG = spm_get_space(PG(1,:));
	spm_smooth(PF(1,:),'./spm_seg_tmp.img',8);
	VF = spm_map('./spm_seg_tmp.img');
	VG = [];
	for i=1:size(PG,1)
		VG = [VG spm_map(PG(i,:))];
	end

	% perform affine normalisation at different sampling frequencies
	% with increasing numbers of parameters.
	[prms,scales] = spm_affine(VG,VF,MG,MF,prms,ones(1,9),1,1.5,8);
	[prms,scales] = spm_affine(VG,VF,MG,MF,prms,ones(1,12),1,0.5,8);
	[prms,scales] = spm_affine(VG,VF,MG,MF,prms,ones(1,12),1,0.5,6);
	[prms,scales] = spm_affine(VG,VF,MG,MF,prms,ones(1,12),1,0.25,4);

	for i=1:size(VG,2)
		spm_unmap(VG(:,i));
	end
	spm_unmap(VF);
	spm_unlink ./spm_seg_tmp.img ./spm_seg_tmp.hdr ./spm_seg_tmp.mat
end
MM = spm_matrix(prms');

for i=1:size(PF,1)
	VF(:,i) = spm_map(PF(i,:));
end

for j=1:size(PB,1)
	VB(:,j)  = spm_map(PB(j,:));
	tmp      = spm_get_space(PB(j,:));
	Mat(:,j) = tmp(:);
end

% Cluster Analysis
%-----------------------------------------------------------------------
% Voxels to sample during the cluster analysis
sxy    = 1/4;
sz     = 6;

%bb1   = [ [-90 91]' [-126 91]' [-72 109]'];
bb1    = [ [-78 78]' [-112 76]' [-50 85]'];

faces = MM*MF*[
1	VF(2,1)/2	VF(3,1)/2	1
VF(1,1)	VF(2,1)/2	VF(3,1)/2	1
VF(1,1)/2	1	VF(3,1)/2	1
VF(1,1)/2	VF(2,1)	VF(3,1)/2	1
VF(1,1)/2	VF(2,1)/2	1	1
VF(1,1)/2	VF(2,1)/2	VF(3,1)	1]';

faces = faces(1:3,:)';
bb    = [max([min(faces); bb1(1,:)]); min([max(faces) ; bb1(2,:)])];
sb    = diff(bb);
d     = round(sb(1:2)*sxy);


lkp=[];
for i=1:size(nc,2)
	lkp = [lkp ones(1,nc(i))*i];
end

n  = size(lkp,2);
m  = size(PF,1);

mn = zeros(m,n);	% Means
cv = zeros(m*m,n);	% (Co)variances
mg = zeros(1,n);	% Number of voxels/cluster

sumbp = zeros(1,n);

spm_progress_bar('Init',niter,'Segmenting','iterations completed');
for iter = 1:niter
	mom0 = zeros(1,n)+eps;
	mom1 = zeros(m,n);
	mom2 = zeros(m*m,n)+eps;

	for p = bb(1,3):sz:bb(2,3)

		B     = spm_matrix([-bb(1,1)*sxy -bb(1,2)*sxy -p 0 0 0 sxy sxy 1]);

		for i=1:m
			% Sample the image data using
			% nearest neighbour to reduce smoothing effects.
			M = inv(B*MM*MF);
			dt = spm_slice_vol(VF(:,i),M, d(1:2),0);
			dat(:,i) = dt(:);
		end
		k = size(dat,1);
		% A priori probability data for eg WM GM CSF scalp etc..
		bp = zeros(k,size(VB,2)+1);
		for j=1:size(VB,2)
			M = inv(B*reshape(Mat(:,j),4,4));
			tmp = spm_slice_vol(VB(:,j),M, d(1:2),1)/nc(j);
			bp(:,j) = tmp(:);
		end
		% Other tissue
		bp(:,size(VB,2)+1) = ...
			abs(ones(k,1) -...
			sum(bp(1:k,1:size(VB,2))')')/nc(size(VB,2)+1);

		if (iter==1)
			% Initial probability estimates based upon
			% a priori knowledge
			for i=1:n
				pr(:,i) = bp(:,lkp(i));
				sumbp(i) = sumbp(i) + sum(bp(:,lkp(i)));
			end
		else
			% Compute PDFs for each cluster
			for i=1:n
				c = reshape(cv(:,i),m,m);
				amp = 1/sqrt((2*pi)^m * det(c)) *mg(i);
				dst = (dat-ones(k,1)*mn(:,i)')/sqrtm(c);
				dst = dst.*dst;
				if (size(dst,2)>1)
					dst = sum(dst')';
				end
				pr(:,i)=amp*exp(-0.5*dst).*bp(:,lkp(i));
			end
		end

		% PDFs must sum to 1
		sp = sum(pr')';
		msk = find(~sp);
		sp(msk) = ones(size(msk));

		for i=1:n
			pr(:,i) = pr(:,i)./sp;

			if (debug ~= 0)
				figure(3);
				subplot(3,3,i);
				image(rot90(reshape(pr(:,i),d(1),d(2))*64));
				drawnow;
			end
		end

		% Compute new n, mean & var for each cluster - step 1
		for i=1:n
			mom0(:,i) = mom0(:,i) + sum(pr(:,i));
			mom1(:,i) = mom1(:,i) + ...
				sum((pr(:,i)*ones(1,m)).*dat)';
			mom2(:,i) = mom2(:,i) + reshape(...
				((pr(:,i)*ones(1,m)).*dat)'*dat,m*m,1);
		end
	end

	% Compute new n, mean & var for each cluster - step 2
	for i=1:n
		mg(:,i) = mom0(:,i)/sumbp(i);
		mn(:,i) = mom1(:,i)/mom0(:,i);
		tmp = (mom0(:,i).*mn(:,i))*mn(:,i)';

		%stabilise the covariance matrix
		tmp = tmp - eye(size(tmp))*eps*10;

		cv(:,i) = (0.000001+mom2(:,i) - tmp(:))/mom0(:,i);
	end

	if iter==1
		nn = 0;
		for j=1:length(nc)
			% Split the clusters
			for i=2:nc(j)
				cv(:,nn+i) = cv(:,nn+i)*0.8^(1-i);
				mn(:,nn+i) = mn(:,nn+i)*0.8^(1-i);
			end
			nn = nn + nc(j);
		end
	end


	% Background Cluster.
	%    Strictly speaking, since voxels contain absolute values,
	%    the distributions should be modified accordingly. However
	%    in practice, the only allowance made for this is for the
	%    distribution of a background cluster. The mean of this
	%    cluster is assumed to be zero, and in order for the model
	%    to fit properly, the number of voxels contained in this 
	%    cluster is doubled.
	mn(:,n)   = zeros(m,1);
	mg(:,n)   = 2*mg(:,n);
	cv(:,n)   = (0.000001+mom2(:,n))/mom0(:,n);

	spm_progress_bar('Set',iter);
end
spm_progress_bar('Clear');



%-----------------------------------------------------------------------

% Create headers, open files etc...
%-----------------------------------------------------------------------
if any(opts == 'r')
	dm     = round(diff(bb)/2);
	vx     = [2 2 2];
	B      = spm_matrix([-bb(1,1)*.5 -bb(1,2)*.5 -bb(1,3)*.5 0 0 0 .5 .5 .5]);
	M      = inv(B*MM);
	orgn   = round(-bb(1,:)/2 + 1);
	planes = bb(1,3):2:bb(2,3);
else
	dm     = VF(1:3);
	vx     = VF(4:6);
	M      = MF;
	orgn   = M\[0 0 0 1]';
	orgn   = round(orgn(1:3)');
	planes = 1:VF(3);
end

k = prod(dm(1:2));

if ~any(opts == 'n')
	for j=1:(size(VB,2)+1)
		iname = [spm_str_manip(PF(1,:),'rd') '_seg' num2str(j)];
		fp(j) = fopen([iname '.img'],'w');
		spm_hwrite([iname '.img'],dm,vx,1/255,2,0,orgn,...
			'Segmented image');
		spm_get_space(iname,M);
	end
end

if any(opts == 'p')
	inamep = [spm_str_manip(PF(1,:),'rd') '_ppet_nsmo'];
	fpp = fopen([inamep '.img'],'w');
	spm_hwrite([inamep '.img'],dm,vx,1/255,2,0,orgn,...
			'Pseudo Pet');
	spm_get_space(inamep,M);
end



% Write the segmented images.
%-----------------------------------------------------------------------
spm_progress_bar('Init',dm(3),'Writing Segmented','planes completed');
clear pr bp dat

for pp=1:size(planes,2)
	p = planes(pp);
	if any(opts == 'r')
		B  = spm_matrix([-bb(1,1)*.5 -bb(1,2)*.5 -p 0 0 0 0.5 0.5 1]);
		M1 = inv(B*MM*MF);
		M2 = B;
	else
		B  = spm_matrix([0 0 -p]);
		M1 = inv(B);
		M2 = B*inv(MM*MF);
	end

	for i=1:m
		% The image data
		tmp = spm_slice_vol(VF(:,i), M1, dm(1:2), 1);
		dat(:,i) = tmp(:);
	end

	bp = zeros(size(dat,1),size(VB,2)+1);
	for j=1:size(VB,2)
		M = inv(M2*reshape(Mat(:,j),4,4));
		tmp = spm_slice_vol(VB(:,j), M, dm(1:2),1);
		bp(:,j) = tmp(:)/nc(j);
	end

	bp(:,size(VB,2)+1) = abs(ones(k,1) -...
		sum(bp(:,1:size(VB,2))')')/nc(size(VB,2)+1);

	for i=1:n
		c = reshape(cv(:,i),m,m);
		amp = 1/sqrt((2*pi)^m * det(c))*mg(i);
		dst = (dat-ones(k,1)*mn(:,i)')/sqrtm(c);
		dst = dst.*dst;
		if (size(dst,2)>1)
			dst = sum(dst')';
		end
		pr(:,i) = amp*exp(-0.5*dst) .* bp(:,lkp(i));
	end

	sp = (sum(pr')'+eps)/255;

	if ~any(opts == 'n')
		nn = 1;
		for j=1:(length(nc))
			dat = pr(:,nn);
			ncj = nc(j);
			if (j==length(nc))
				ncj = ncj - 1;
			end;
			for i=2:ncj
				dat = dat + pr(:,(nn+i-1));
			end
			nn = nn + ncj;
			dat = round(dat./sp);
			fwrite(fp(j),dat,'uchar');
		end
	end

	if any(opts == 'p')
		nn = 1;
		dat = round((pr*petscales)./sp);
		fwrite(fpp,dat,'uchar');
	end

	spm_progress_bar('Set',pp);
end
spm_progress_bar('Clear');

spm_figure('Clear','Graphics');
if ~any(opts == 'n')
	for j=1:(size(VB,2)+1)
		fclose(fp(j));
	end
end
if any(opts == 'p')
	fclose(fpp);

	iname = [spm_str_manip(PF(1,:),'rd') '_ppet'];
	spm_smooth([inamep '.img'],[iname '.img'], 8);
	spm_unlink([inamep '.img'],[inamep '.hdr'],[inamep '.mat']);
end

for v=[VF VB]
	spm_unmap(v);
end
