function spm_realign(P,Flags)
% within mode image realignment routine
%
% FORMAT spm_realign
% With no arguments, spm_realign acts as a user interface for setting
% up the realignment parameters.
%
%____________________________________________________________________________
% The prompts are:
% 'number of subjects'
% Enter the number of subjects you wish to realign.
%
% 'select scans for subject ..'
% Select the scans you wish to realign. All operations are relative
% to the first image selected.
%
% 'Which option?'
% 	'Coregister only'
% 	Only determine the parameters required to transform each of the
% 	images 2..n to the same space as image 1.
% 	The determined parameters for image XXXX.img are saved in the
%	file XXXX.mat. This is a Matlab file, containing the matrix 'M'.
% 	The location of an image voxel (in mm) can be determined by
% 	computing M(1:3,:)*[xcoord ycoord zcoord 1].
%	Note that if the coregistration is performed more than once on
%	the unresliced data, the starting estimates are obtained from
% 	parameters stored in the '.mat' files.
%
% 	'Reslice Only'
% 	Reslice the specified images according to the contents of the
% 	previously determined parameters. The images are resliced to be
% 	in the same space as the first one selected.
%
% 	'Coregister & Reslice'
% 	Combine the above two steps together.
%
% 'Interpolation Method?'
% 	'Bilinear Interpolation'
% 	Use bilinear interpolation to sample the images when
% 	determining parameters/reslicing.
%
% 	'Sinc Interpolation'
% 	Use a modified sinc interpolation to sample the images.
% 	This is slower than bilinear interpolation, but produces better
% 	results. It is especially recommended for fMRI time series.
%
% Options for coregistering:
%
% 'Coregister 2 with 1 only?'
% This is for if you wish to register a group of scans, to a scan of the
% same subject which was performed at a different time.
% This option allows the user to determine the parameters which would
% transform image 2 to the space of image 1, and apply same transformation
% parameters to images 3..n.
%
% 'More iterations?'
% By default, the coregistration uses 5 iterations. This is usually fine for
% simple subject movement when in the scanner.
% However, it is sometimes desirable to realign two data sets aquired at
% different times, with large differences in subject positioning.
% This option allows for 25 iterations of the algorithm to be performed,
% which should be enough to correct large displacements.
%
%
% Options for reslicing:
%
% 'Create what?'
% 	'All Images (1..n)'
% 	This reslices all the images - including the first image selected
% 	- which will remain in it's original position.
%
%	'Images 2..n'
% 	Reslices images 2..n only. Useful for if you wish to reslice
% 	(for example) a PET image to fit a structural MRI, without
% 	creating a second identical MRI volume.
%
%	'All Images + Mean Image'
% 	In addition to reslicing the images, it also creates a mean of the
% 	resliced image.
%
%	'Mean Image Only'
% 	Creates the mean image only.
%
% 'Mask the images?'
% To avoid artifactual movement-related variance.
% Because of subject motion, different images are likely to have different
% patterns of zeros from where it was not possible to sample data.
% With masking enabled, the program searches through the whole time series
% looking for voxels which need to be sampled from outside the original
% images. Where this occurs, that voxel is set to zero for the whole set
% of images.
%
% 'Adjust for motion in Z?' (fMRI only)
% Adjust the data (fMRI) to remove movement-related components.
% The adjustment procedure is based on a autoregression-moving average-like
% model of the effect of position on signal and explicitly includes a spin
% excitation history effect. Do not use this option if you are reslicing
% any less than about 20 images.
%
%____________________________________________________________________________
% TO OBTAIN SIMILAR RESULTS TO SPM95, select:
% 'Which option?'        	'Coregister & Reslice'
% 'Interpolation Method?'	'Bilinear Interpolation' (for PET)
% 'Interpolation Method?'	'Sinc Interpolation'     (for fMRI)
%
% 'Coregister 2 with 1 only?'	'NO'
% 'More iterations?'         	'NO'
%
% 'Create what?'        	'All Images + Mean Image'
% 'Mask the images?'       	'YES'
% 'Adjust for motion in Z?'	'YES' (for fMRI)
%____________________________________________________________________________
%
% FORMAT spm_realign(P, Flags)
% P     - matrix of filenames {one string per row}
%         All operations are performed relative to the first image.
%         ie. Coregistration is to the first image, and resampling
%         of images is into the space of the first image.
%
% Flags - options flags
%         e - compute realignment parameters
%             Realigns a time-series of scans whose filenames are
%             in the input matrix (P).
%             The realignment uses least squares solutions and partial
%             derivatives of the image process with respect to the spatial 
%             transformations constituting a rigid-body transformation.
%             These movement parameter estimates are based on a first order
%             Taylor expansion in terms of the partial derivatives measured 
%             using smoothed scans.
%             Subsequent scans are realigned with the first.
%
%         E - same as "e" but with more iterations.
%             For correcting more than just patient movement.
%
%         a - only register image 2 with image 1, and apply the 
%             transformations to the subsequent images (which should be in
%             the same space as image 1).
%
%         m - write transformation matrixes
%             Matrixes are written which define the space of the images.
%             From these matrixes it is possible to reslice any of the
%             images, to the same space as any other image.
%
%         r - reslice images
%             The spatially realigned and adjusted images are written to
%             the orginal subdirectory with the same filename but prefixed
%             with a 'r'.
%             They are all aligned with the first.
%
%         c - adjust the data (fMRI) to remove movement-related components
%             The adjustment procedure is based on a autoregression-moving 
%             average-like model of the effect of position on signal and 
%             explicitly includes a spin excitation history effect.
%
%         k - mask output images
%             To avoid artifactual movement-related variance the realigned
%             set of images can be internally masked, within the set (i.e.
%             if any image has a zero value at a voxel than all images have
%             zero values at that voxel).  Zero values occur when regions
%             'outside' the image are moved 'inside' the image during
%             realignment.
%
%         i - write mean image
%             The average of all the realigned scans is written to
%             mean*.img.
%
%         s - use sinc interpolation for parameter estimates (7x7x7).
%
%         S - use sinc interpolation for reslicing (11x11x11).
%
%         n - don't reslice the first image
%             The first image is not actually moved, so it may not be
%             necessary to resample it.
%
%         N - don't reslice any of the images - except possibly create a
%             mean image.
%____________________________________________________________________________
%
%
% Refs:
%
% Friston KJ, Ashburner J, Frith CD, Poline J-B, Heather JD & Frackowiak
% RSJ (1995) Spatial registration and normalization of images Hum. Brain
% Map. 2:165-189
%
% Friston KJ, Williams SR, Howard R Frackowiak RSJ and Turner R (1995)
% Movement-related effect in fMRI time-series.  Mag. Res. Med. 35:346-355
%
%__________________________________________________________________________
% %W% Karl Friston - modified by John Ashburner %E%


global MODALITY

if (nargin == 0)
	spm_figure('Clear','Interactive');
	set(spm_figure('FindWin','Interactive'),'Name','Realignment')

	n     = spm_input('number of subjects',1);
	for i = 1:n
		P = spm_get(Inf,'.img',...
			['select scans for subject ' num2str(i)]);
		eval(['P' num2str(i) ' = P;']);
	end

	Flags = 'm';
	p = spm_input('Which option?',1,'m',...
		'Coregister only|Reslice Only|Coregister & Reslice',...
		[1 2 3]);
	if (p == 1 | p == 3) Flags = [Flags 'e']; end
	if (p == 2 | p == 3) Flags = [Flags 'r']; end
	p = spm_input('Interpolation Method?',2,'m',...
		'Bilinear Interpolation|Sinc Interpolation',[1 2]);
	if (p == 2) Flags = [Flags 'sS']; end

	if (any(Flags == 'e'))
		if (spm_input('Coregister 2 with 1 only?',3,'y/n') == 'y')
			Flags = [Flags 'a']; end
		if (spm_input('More iterations?'         ,4,'y/n') == 'y')
			Flags = [Flags 'E']; end
	end

	if (any(Flags == 'r'))
		p = spm_input('Create what?',5,'m',...
			[' All Images (1..n)| Images 2..n|'...
			 ' All Images + Mean Image| Mean Image Only'],...
			[1 2 3 4]);
		if (p == 2) Flags = [Flags 'n']; end
		if (p == 3) Flags = [Flags 'i']; end
		if (p == 4) Flags = [Flags 'Ni']; end
		if (~any(Flags == 'N'))
			if (spm_input('Mask the images?'       ,6,'y/n')...
				== 'y') Flags = [Flags 'k']; end
			if (strcmp(MODALITY,'FMRI'))
				if (spm_input(...
					'Adjust for motion in Z?'...
					,7,'y/n') == 'y')
					Flags = [Flags 'c'];
				end
			end
		end
	end

	set(spm_figure('FindWin','Interactive'),'Name','executing',...
		'Pointer','Watch'); drawnow
	for i = 1:n
		eval(['P = P' num2str(i) ';'])
		spm_realign(P,Flags);
	end
	spm_figure('Clear','Interactive');
	set(spm_figure('FindWin','Interactive'),'Name','',...
		'Pointer','Arrow');
	drawnow
	return;
end


%---------------------------------------------------------------------------
global PRINTSTR

% set up file identifiers and check image and voxel sizes are consistent
%---------------------------------------------------------------------------
Q     = zeros(size(P,1),6);		% initial estimate/default values

nreg = size(P,1);
if any(Flags == 'a') nreg = 2; end

% Computation of Realignment Parameters - uses P & Q
%---------------------------------------------------------------------------
if (any(Flags == 'e') | any(Flags == 'E'))

	Hold = 1;
	if (any(Flags == 's')) Hold = 3; end

	spm_smooth(spm_str_manip(P(1,:),'d'),'spm_ref.img',8);
	V1 = spm_map('spm_ref.img');

	% center, bounding box and margins
	%-------------------------------------------------------------------
	bb    = spm_bb(P(1,:));		% bb = bounding box
	sb    = diff(bb);		% size of bb


	% compute matrices
	% dQ  = 6 ortholinear transformation parameters
	%-------------------------------------------------------------------
	a     = pi/180;			% dQ (rotation) radians
	b     = 1;			% dQ (translation) mm
	dQ    = [b 0 0 0 0 0;		% unit transformations
 		 0 b 0 0 0 0;
		 0 0 b 0 0 0;
	         0 0 0 a 0 0;
	         0 0 0 0 a 0;
	         0 0 0 0 0 a];


	% define height of transverse slices (S) used in
	% subsampling the volume
	%-------------------------------------------------------------------
	% transverse
	S     = linspace((bb(1,3) + 6/V1(6)),(bb(2,3) - 8/V1(6)),8);	
	h     = 5;				% number of recursions
	if (any(Flags == 'E')) h = 25; end;
	M     = sb(1);				% rows per slice
	N     = sb(2);				% columns per slice
	n     = M*N;				% voxels per slice

	U = [1 2 3 4 5 6];
	if V1(3,1) == 1; U = [1 2 6]; end	% 3 params for slices

	if sb(3) == 0; S = 1; end		% 1 section for slices

	C1 = spm_get_space(spm_str_manip(P(1,:), 'd'));

	% compute X (the reference image) and dX/dQ
	% (the effects of moving X).
	%-------------------------------------------------------------------
	X     = zeros(n*length(S),1);
	Y     = zeros(n*length(S),1);
	dXdQ  = zeros(n*length(S),size(U,2));
	for i = 1:length(S)
		B     = spm_matrix([-bb(1,1) -bb(1,2) -S(i) 0 0 0 1 1 1]);
		x     = spm_slice_vol(V1,inv(B),[M N], Hold);
	  	X([1:n] + (i - 1)*n) = x(:);
	        for j = 1:size(U,2)
			A      = spm_matrix(dQ(U(j),:));
			d      = spm_slice_vol(V1,inv(B*inv(C1)*A*C1),...
				[M N],Hold);
			dX     = d - x;
			dXdQ(([1:n] + (i - 1)*n),j) = dX(:);
		end
	end


	% least squares solution for Q the movements where:
	% Y = X + dX/dQ.Q  => Q = [-dX/dQ Y]\X
	% The estimate is repeated iteratively h times and the estimates of
	% Q are cumulated arithmetically (an aproximation but good enough)
	%-------------------------------------------------------------------
	spm_progress_bar('Init',nreg-1,'Coregistering','volumes completed');

	for k = 2:nreg

	    % decided against basing starting estimates on
	    % solution for previous run
	    q = zeros(1,size(dQ,1));

	    C2 = spm_get_space(spm_str_manip(P(k,:), 'd'));

	    spm_smooth(spm_str_manip(P(k,:),'d'),'spm_mov.img',8);
	    V2 = spm_map('spm_mov.img');

	    for i = 1:h
	  	for j = 1:length(S)
			B      = spm_matrix(...
				[-bb(1,1) -bb(1,2) -S(j) 0 0 0 1 1 1]);
			y      = spm_slice_vol(V2,...
				inv(B*inv(C1)*spm_matrix(q)*C2),[M N],Hold);
			Y(([1:n] + (j - 1)*n)) = y(:);
	  	end
		q0 = [-dXdQ Y]\X;
		q0 = q0(1:size(U,2));
	 	q      = q - q0'*dQ(U,:);
	  	spm_progress_bar('Set',k-2+i/h);
	    end
	    Q(k,:) = q;

	    spm_unmap(V2); spm_unlink spm_mov.img spm_mov.hdr spm_mov.mat
	end
	Q((nreg+1):size(P,1),:) = ones(size(P,1)-nreg,1)*Q(nreg,:);
	spm_unmap(V1); spm_unlink spm_ref.img spm_ref.hdr spm_ref.mat
	spm_figure('Clear','Interactive');

	% display results
	% translation and rotation over time series
	%-------------------------------------------------------------------
	figure(spm_figure('FindWin','Graphics'));
	spm_figure('Clear','Graphics');
	subplot(3,1,1);
	axis off
	title('Image realignment','FontSize',16,'FontWeight','Bold')
	x     = -0.1;
	y     =  0.9;
	for i = 1:min([size(P,1) 12])
		text(x,y,[sprintf('%-4.0f',i) P(i,:)],'FontSize',10);
		y = y - 0.08;
	end
	if size(P,1) > 12
		text(x,y,'................ etc','FontSize',10); end


	subplot(3,1,2)
	plot(Q(:,1:3))
	s = ['x translation';'y translation';'z translation'];
	text([2 2 2], Q(2, 1:3), s, 'Fontsize',10)
	title('translation','FontSize',16,'FontWeight','Bold')
	xlabel('image')
	ylabel('mm')
	grid 


	subplot(3,1,3)
	plot(Q(:,4:6)*180/pi)
	s = ['pitch';'roll ';'yaw  '];
	text([2 2 2], Q(2, 4:6)*180/pi, s, 'Fontsize',10)
	title('rotation','FontSize',16,'FontWeight','Bold')
	xlabel('image')
	ylabel('degrees')
	grid

	% print realigment parameters and figure 2 attributes
	spm_print
end



% Application of Realignment Parameters - uses P & Q
%---------------------------------------------------------------------------
if any(Flags == 'r')

	Hold = 1;
	if (any(Flags == 'S')) Hold = 5; end

	% Get properties of image to realign to
	%-------------------------------------------------------------------
	p = spm_str_manip(P(1,:), 'd');
	M = spm_get_space(p);
	[DIM VOX SCALE TYPE OFFSET ORIGIN] = spm_hread(p);

	% Get properties of all the images
	%-------------------------------------------------------------------
	Headers = zeros(size(P,1),(3+1+1));
	Matrixes = zeros(16,size(P,1));
	V     = zeros(12,size(P,1));
	for i = 1:size(P,1)
		p  = spm_str_manip(P(i,:), 'd');

		A  = spm_matrix(Q(i,:));
		C2 = spm_get_space(p);
		M1  = inv(M)*A*C2;
		Matrixes(:,i) = M1(:);

		[dim dummy scale typ dummy dummy] = spm_hread(p);
		Headers(i,:) = [dim scale typ];

		V(:,i) = spm_map(spm_str_manip(P(i,:),'d'));
	end


	spm_progress_bar('Init',DIM(3),'Reslicing','planes completed');

	Xm = kron(ones(DIM(2),1),[1:DIM(1)]');
	Ym = kron([1:DIM(2)]',ones(DIM(1),1));
	open_mode = 'w';
	if (any(Flags == 'i') | any(Flags == 'c'))
		Integral = zeros(prod(DIM(1:2)),DIM(3));
	end

	start_vol = 1;
	if (any(Flags == 'n') & ~any(Flags == 'i')) start_vol = 2; end

	for j = 1:DIM(3)
		B     = spm_matrix([0 0 -j 0 0 0 1 1 1]);
		if (any(Flags == 'c'))
			X = zeros(DIM(1)*DIM(2),size(P,1));
		end

		% get masks for this plane
		%-----------------------------------------------------------
		if (any(Flags == 'i') | any(Flags == 'k') |...
			any(Flags == 'c'))
			Count = zeros(prod(DIM(1:2)),1);
			for i = 1:size(P,1)
				M0  = reshape(Matrixes(:,i),4,4);
				M1 = inv(B*M0);

				tiny = 5e-2; % From spm_vol_utils.c
				% Find range of slice
				tmp = M1(1,1)*Xm + M1(1,2)*Ym +...
					(M1(1,3)*0 + M1(1,4));
				Mask = (tmp >= (1-tiny) &...
					 tmp <= (Headers(i,1)+tiny));

				tmp  = M1(2,1)*Xm + M1(2,2)*Ym +...
					(M1(2,3)*0 + M1(2,4));
				Mask = Mask & (tmp >= (1-tiny) &...
					tmp <= (Headers(i,2)+tiny));

				tmp = M1(3,1)*Xm + M1(3,2)*Ym +...
					(M1(3,3)*0 + M1(3,4));
				Mask = Mask & (tmp >= (1-tiny) &...
					tmp <= (Headers(i,3)+tiny));

				Count = Count + Mask;
			end
			Mask = (Count == size(P,1));
		end

		% get values for this plane
		%-----------------------------------------------------------
		for i = start_vol:size(P,1)
			M0  = reshape(Matrixes(:,i),4,4);
			M1 = inv(B*M0);

			d  = spm_slice_vol(V(:,i),M1,DIM(1:2),Hold);

			if (any(Flags == 'i') | any(Flags == 'c'))
				Integral(:,j) = Integral(:,j) + d(:);
			end

			if (~any(Flags == 'c'))
				% don't need to load the whole time series.
				if (i > 1 | ~any(Flags == 'n'))
					d = d/Headers(i,4);
					if any(Flags == 'k')
						d = d(:).*Mask;
					end

	% Deal with different data types
	if (Headers(i,5) == 2)
		d = round(d);
		tmp = find(d > 255);
		d(tmp) = zeros(size(tmp))+255;
		tmp = find(d < 0);
		d(tmp) = zeros(size(tmp));
	elseif (Headers(i,5) == 4)
		d = round(d);
		tmp = find(d > 32767);
		d(tmp) = zeros(size(tmp))+32767;
		tmp = find(d < -32768);
		d(tmp) = zeros(size(tmp))-32768;
	elseif (Headers(i,5) == 8)
		d = round(d);
		tmp = find(d > 2^31-1);
		d(tmp) = zeros(size(tmp))+(2^31-1);
		tmp = find(d < -2^31);
		d(tmp) = zeros(size(tmp))-2^31;
	end

				end

				if (~any(Flags == 'N'))
					p  = spm_str_manip(P(i,:), 'd');
					q  = max([find(p == '/') 0]);
					q  = [p(1:q) 'r' p((q + 1):length(p))];
					fp = fopen(q, open_mode);
					fwrite(fp,d,spm_type(Headers(i,5)));
					fclose(fp);
				end
			else
				X(:,i) = d(:);
			end
		end

		if (any(Flags == 'i') | any(Flags == 'c'))
			Integral(:,j) = Integral(:,j)./(Count + eps);
		end


		if (any(Flags == 'c'))
			% Adjust for motion in Z
			%---------------------------------------------------
			Mask = (Count == size(P,1));

			X = X'; % Transpose to reduce paging
			lzm1 = size(P,1)-1;
			dZmotion = Matrixes(3,:)' - Matrixes(3,1);
			for y=1:DIM(2)
				Zmotiony  = Matrixes(7,:)*y +...
					Matrixes(11,:)*j + Matrixes(15,:);
				Zmotiony  = Zmotiony' - Zmotiony(1);
				for x = 1:DIM(1)
					xy = x+(y-1)*DIM(1);
					if (Mask(xy))

	Tac = X(:,xy) - Integral(xy,j);
	Zmotion = Zmotiony + dZmotion*x;
	A = [Zmotion Zmotion.^2 [0 ; Zmotion(1:lzm1)] ...
		[0 ; Zmotion(1:lzm1)].^2];
	X(:,xy) = X(:,xy) - A*((A'*A)\(A'*Tac));

					end
				end
			end
			X = X';
			start_vol = 1;
			if (any(Flags == 'n')) start_vol = 2; end

			for i = start_vol:size(P,1)
				p = deblank(P(i,:));

				d = X(:,i)/Headers(i,4);
				if any(Flags == 'k') d = d.*Mask; end

				% Deal with different data types
				if (Headers(i,5) == 2)
					d = round(d);
					tmp = find(d > 255);
					d(tmp) = zeros(size(tmp))+255;
					tmp = find(d < 0);
					d(tmp) = zeros(size(tmp));
				elseif (Headers(i,5) == 4)
					d = round(d);
					tmp = find(d > 32767);
					d(tmp) = zeros(size(tmp))+32767;
					tmp = find(d < -32768);
					d(tmp) = zeros(size(tmp))-32768;
				elseif (Headers(i,5) == 8)
					d = round(d);
					tmp = find(d > 2^31-1);
					d(tmp) = zeros(size(tmp))+(2^31-1);
					tmp = find(d < -2^31);
					d(tmp) = zeros(size(tmp))-2^31;
				end

				q      = max([find(p == '/') 0]);
				q      = [p(1:q) 'r' p((q + 1):length(p))];
				fp = fopen(q, open_mode);
				fwrite(fp,d,spm_type(Headers(i,5)));
				fclose(fp);
			end
		end
		open_mode = 'a';
		spm_progress_bar('Set',j);
	end

	% Write headers and matrixes
	%------------------------------------------------------------------
	if (~any(Flags == 'N'))
		for i = start_vol:size(P,1)
			p = spm_str_manip(P(i,:), 'd');
			q = max([find(p == '/') 0]);
			q = [p(1:q) 'r' p((q + 1):length(p))];
			spm_hwrite(q,DIM,VOX,Headers(i,4),Headers(i,5),...
				0,ORIGIN,'spm - realigned');
			if (any(Flags == 'm')) spm_get_space(q,M); end
		end
	end

	if any(Flags == 'i')
		% Write integral image (16 bit signed)
		%-----------------------------------------------------------
		mx = max(max(Integral));
		SCALE  = mx/32767;
		p      = P(1,:);
		p      = p(p ~= ' ');
		q      = max([find(p == '/') 0]);
		q      = [p(1:q) 'mean' p((q + 1):length(p))];
		fp = fopen(q,'w');
		for j = 1:DIM(3)
			d = round(Integral(:,j)/SCALE);
			tmp = find(d > 32767);
			d(tmp) = zeros(size(tmp))+32767;
			tmp = find(d < -32768);
			d(tmp) = zeros(size(tmp))-32768;
			fwrite(fp,d,spm_type(4));
		end
		spm_hwrite(q,DIM,VOX,SCALE,4,0,ORIGIN,'spm - mean image');
		if (any(Flags == 'm')) spm_get_space(q,M); end
	end

	for i = 1:size(P,1); spm_unmap(V(:,i)); end
end

% Saving of Realignment Parameters - uses P & Q
%---------------------------------------------------------------------------
if (any(Flags == 'm') & any(Flags == 'e'))
	for k = 1:size(P,1)
		C2 = spm_get_space(spm_str_manip(P(k,:), 'd'));
		spm_get_space(P(k,:), spm_matrix(Q(k,:))*C2);
	end
end

spm_figure('Clear','Interactive');
