% Segment a T1 weighted MR image into Gray White & CSF.
% FORMAT spm_segment(PF)
% PF - name(s) of image(s) to segment
% ___________________________________________________________________________
%
%                      The algorithm is three step:
%
% 1) Determine the affine transform which best matches the T1 image with a
%    template image. If the name of more than one image is passed, then the
%    first image is used in this step.
%
% 2) Perform Cluster Analysis with a modified Mixture Model and a-priori
%    information about the likelihoods of each voxel being one of a number
%    of different data types. If more than one image is passed, then they
%    they are all assumed to be in register, and the voxel values are fitted
%    to multi-normal distributions.
%
% 3) Write the segmented image. The names of these images have
%    "_seg1", "_seg2" & "_seg3" appended to the name of the first image passed.
%
% ___________________________________________________________________________
%
% The template image, and a-priori likelihood images are modified versions of
% those kindly supplied by the MNI group in Montreall.

% %W% (c) John Ashburner %E%

function spm_segment(PF)

DIR   = [SWD '/mni/'];

%- A-Priori likelihood images
PB    = str2mat([DIR 'spgray.img'],[DIR 'spwhite.img'],[DIR 'spcsf.img']);

%- Template image
PG    = [DIR 'scolin.img'];
niter = 24;
nc    = [1 1 1 3]; % Number of clusters for each probability image

MF = spm_get_space(PF(1,:));
if (~isempty(PG))
	% Affine normalisation
	MG = spm_get_space(PG(1,:));
	spm_smooth(PF(1,:),'spm_seg_tmp.img',8);
	VF = spm_map('spm_seg_tmp.img');
	VG = [];
	for i=1:size(PG,1)
		VG = [VG spm_map(PG(i,:))];
	end
	Affine = spm_affine(VG,VF,MF\MG);
	for i=1:size(VG,2)
		spm_unmap(VG(:,i));
	end
	spm_unmap(VF);
	delete spm_seg_tmp.img spm_seg_tmp.hdr
else
	% Image already normalised
	MG = eye(4);
	Affine = MF\MG;
end

for i=1:size(PF,1)
	VF(:,i) = spm_map(PF(i,:));
end

for j=1:size(PB,1)
	VB(:,j)  = spm_map(PB(j,:));
	tmp = spm_get_space(PB(j,:));
	Mat(:,j) = tmp(:);
end

% Sampling intervals
sxy = 1/4;
sz  = 6;
bb     = [ [-70 70]' [-115 80]' [-28 72]' ];
sb    = diff(bb);
d = round(sb(1:2)*sxy);

if exist('pr')==1 clear pr dat pb ; end

lkp=[];
for i=1:size(nc,2)
	lkp = [lkp ones(1,nc(i))*i];
end
n  = size(lkp,2);
m  = size(PF,1);
mn = zeros(m,n);
cv = zeros(m*m,n);
mg = zeros(1,n);

sumbp = zeros(1,n);

for iter = 1:niter
	fprintf('iteration %d', iter);
	mom0 = zeros(1,n)+eps;
	mom1 = zeros(m,n);
	mom2 = zeros(m*m,n)+eps;

	for p = bb(1,3):sz:bb(2,3)

		B     = spm_matrix([-bb(1,1)*sxy -bb(1,2)*sxy -p 0 0 0 sxy sxy 1]);

		for i=1:m
			% The image data
			dt = spm_slice_vol(VF(:,i),inv(B*MF), d(1:2),0);
			dat(:,i) = dt(:);
		end
		k = size(dat,1);
		% A priori probability data for eg WM GM CSF scalp etc..
		bp = zeros(k,size(VB,2)+1);
		for j=1:size(VB,2)
			M = inv(B*MF*Affine*inv(MG)*reshape(Mat(:,j),4,4));
			tmp = spm_slice_vol(VB(:,j),M, d(1:2),1)/nc(j);
			bp(:,j) = tmp(:);
		end
		% Other tissue
		bp(:,size(VB,2)+1) = abs(ones(k,1) - sum(bp(1:k,1:size(VB,2))')')/nc(size(VB,2)+1);

		if (iter==1)
			% Initial probability estimates based upon a priori knowledge
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
				pr(:,i) = amp*exp(-0.5*dst) .* bp(:,lkp(i));
			end
		end

		% PDFs must sum to 1
		sp = sum(pr')';
		msk = find(~sp);
		sp(msk) = ones(size(msk));

		for i=1:n
			pr(:,i) = pr(:,i)./sp;
			%subplot(3,3,i);image(rot90(reshape(pr(:,i),d(1),d(2))*64));drawnow;
		end

		% Compute new n, mean & var for each cluster - step 1
		for i=1:n
			mom0(:,i) = mom0(:,i) + sum(      pr(:,i));
			mom1(:,i) = mom1(:,i) + sum(     (pr(:,i)*ones(1,m)).*dat)';
			mom2(:,i) = mom2(:,i) + reshape(((pr(:,i)*ones(1,m)).*dat)'*dat,m*m,1);
		end
		fprintf('.');
	end

	% Compute new n, mean & var for each cluster - step 2
	for i=1:n
		mg(:,i) = mom0(:,i)/sumbp(i);
		mn(:,i) = mom1(:,i)/mom0(:,i);
		tmp = (mom0(:,i).*mn(:,i))*mn(:,i)';
		tmp = tmp - eye(size(tmp))*eps*10; %stabilise the covariance matrix
		cv(:,i) = (0.000001+mom2(:,i) - tmp(:))/(mom0(:,i)-1);
	end

	if iter==1
		% Background
		mn(:,n)=zeros(m,1);

		nn = 0;
		for j=1:length(nc)
			% Split the clusters
			for i=2:nc(j)
				cv(:,nn+i) = cv(:,nn+i)*0.3^(i-1);
				mn(:,nn+i) = mn(:,nn+i)*0.3^(i-1);
			end
			nn = nn + nc(j);
		end
	end		
	fprintf('\n');
end

% Create headers, open files etc...
[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(PF(1,:));
M = MF;
for j=1:size(VB,2)
	iname = [spm_str_manip(PF(1,:),'rd') '_seg' num2str(j)];
	fp(j) = fopen([iname '.img'],'w');
	spm_hwrite([iname '.img'],DIM,VOX,1/255,2,0,ORIGIN,'Segmented image');
	save(iname,'M');
end

% Write the segmented images.
clear pr bp dat
k = prod(VF(1:2));
for p=1:VF(3)
	for i=1:m
		% The image data
		tmp = spm_slice_vol(VF(:,i),spm_matrix([0 0 p]), VF(1:2),0);
		dat(:,i) = tmp(:);
	end

	bp = zeros(size(dat,1),size(VB,2)+1);
	for j=1:size(VB,2)
		M = inv(spm_matrix([0 0 -p])*Affine*inv(MG)*reshape(Mat(:,j),4,4));
		tmp = spm_slice_vol(VB(:,j),M, VF(1:2),1)/nc(j);
		bp(:,j) = tmp(:);
	end
	bp(:,size(VB,2)+1) = abs(ones(k,1) - sum(bp(1:k,1:size(VB,2))')')/nc(size(VB,2)+1);

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

	sp = sum(pr')'/255;

	nn = 1;
	for j=1:(length(nc)-1)
		dat = pr(:,nn);
		for i=2:nc(j)
			dat = dat + pr(:,(nn+i-1));
		end
		nn = nn + nc(j);
		dat = round(dat./sp);
		fwrite(fp(j),dat,'uchar');
	end

	fprintf('.');
end
fprintf('\n');

for j=1:size(VB,2)
	fclose(fp(j));
end

for v=[VF VB]
	spm_unmap(v);
end
