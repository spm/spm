function spm_realign(P,Flags)
% within mode image realignment routine
% FORMAT spm_realign(P, Flags)
% P     - matrix of filenames {one string per row}
%         All operations are performed relative to the first image.
%         ie. Coregistration is to the first image, and resampling
%         of images is into the space of the first image.
%
% Flags - options flags
%         e - compute realignment parameters
%             Realigns a time-series of scans whose filenames are in the input matrix (P).
%             The realignment uses least squares solutions and partial derivatives of the
%             image process with respect to the spatial transformations constituting a
%             rigid-body transformation.  These movement parameter estimates are based
%             These movement parameter estimates are based on a first order Taylor expansion
%             in terms of the partial derivatives measured using smoothed scans.
%             Subsequent scans are realigned with the first.
%
%         E - same as "e" but with more iterations.
%             For correcting more than just patient movement.
%
%         m - write transformation matrixes
%             Matrixes are written which define the space of the images. From these matrixes
%             it is possible to reslice any of the images, to the same space as any other image.
%
%         r - reslice images
%             The spatially realigned and adjusted images are written to
%             the orginal subdirectory with the same filename but prefixed with a 'r'.
%             They are all aligned with the first.
%
%         c - adjust the data (fMRI) to remove movement-related components
%             The adjustment procedure is based on a autoregression-moving average
%             -like model of the effect of position on signal and explicitly includes
%             a spin excitation history effect.
%
%         k - mask output images
%             To avoid artifactual movement-related variance the realigned set of images
%             are internally masked, within the set (i.e. if any image has a zero value at
%             a voxel than all images have zero values at that voxel).  Zero values
%             occur when regions 'outside' the image are moved 'inside' the image during
%             realignment.
%
%         i - write mean image
%             The average of all the realigned scans is written to mean*.img.
%
%         s - use sinc interpolation for parameter estimates (7x7x7).
%
%         S - use sinc interpolation for reslicing (11x11x11).
%
%         n - don't reslice the first image
%             The first image is not actually moved, so it may not be necessary to
%             resample it.
%
%         N - don't reslice any of the images - except possibly create a mean image.
%
%____________________________________________________________________________
%
%
% Refs:
%
% Friston KJ, Ashburner J, Frith CD, Poline J-B, Heather JD & Frackowiak
% RSJ (1995) Spatial registration and normalization of images Hum. Brain
% Map. 0:00-00
%
% Friston KJ, Williams SR, Howard R Frackowiak RSJ and Turner R (1995)
% Movement-related effect in fMRI time-series.  Mag. Res. Med. 0:00-00
%
%__________________________________________________________________________
% %W% Karl Friston - modified by John Ashburner %E%


%----------------------------------------------------------------------------
global PRINTSTR

% set up file identifiers and check image and voxel sizes are consistent
%----------------------------------------------------------------------------

Q     = zeros(size(P,1),6);			% initial estimate/default values


% Computation of Realignment Parameters - uses P & Q
%----------------------------------------------------------------------------
if (any(Flags == 'e') | any(Flags == 'E'))

	Hold = 1;
	if (any(Flags == 's')) Hold = 3; end

	spm_smooth(spm_str_manip(P(1,:),'d'),'spm_ref.img',8);
	V1 = spm_map('spm_ref.img');

	% center, bounding box and margins
	%----------------------------------------------------------------------------
	bb    = spm_bb(P(1,:));					% bb = bounding box
	sb    = diff(bb);					% size of bb


	% compute matrices
	% dQ  = 6 ortholinear transformation parameters
	%----------------------------------------------------------------------------
	a     = pi/180;						% dQ (rotation) radians
	b     = 1;						% dQ (translation) mm
	dQ    = [b 0 0 0 0 0;					% unit transformations
 		 0 b 0 0 0 0;
		 0 0 b 0 0 0;
	         0 0 0 a 0 0;
	         0 0 0 0 a 0;
	         0 0 0 0 0 a];


	% define height of transverse slices (S) used in subsampling the volume
	%----------------------------------------------------------------------------
	S     = linspace((bb(1,3) + 6/V1(6)),(bb(2,3) - 8/V1(6)),8);	% transverse
	h     = 5;						% number of recursions
	if (any(Flags == 'E')) h = 25; end;
	M     = sb(1);						% rows per slice
	N     = sb(2);						% columns per slice
	n     = M*N;						% voxels per slice

	if V1(3,1) == 1; dQ = dQ([1 2 6],:); end			% 3 params for slices
	if sb(3) == 0; S = 1; end				% 1 section for slices

	C1 = spm_get_space(spm_str_manip(P(1,:), 'd'));

	% compute X (the reference image) and dX/dQ (the effects of moving X)
	%----------------------------------------------------------------------------
	X     = zeros(n*length(S),1);
	Y     = zeros(n*length(S),1);
	dXdQ  = zeros(n*length(S),size(dQ,1));
	for i = 1:length(S)
		B     = spm_matrix([-bb(1,1) -bb(1,2) -S(i) 0 0 0 1 1 1]);
		x     = spm_slice_vol(V1,inv(B),[M N], Hold);
	  	X([1:n] + (i - 1)*n) = x(:);
	        for j = 1:size(dQ,1);
			A      = spm_matrix(dQ(j,:));
			d      = spm_slice_vol(V1,inv(B*inv(C1)*A*C1),[M N],Hold);
			dX     = d - x;
			dXdQ(([1:n] + (i - 1)*n),j) = dX(:);
		end
	end


	% least squares solution for Q the movements where:
	% Y = X + dX/dQ.Q  => Q = [-dX/dQ Y]\X
	% The estimate is repeated iteratively h times and the estimates of Q
	% are cumulated arithmetically (an aproximation but good enough)
	%---------------------------------------------------------------------------
	figure(2);
	delete(get(2,'Children'));
	ax = axes('Position', [0.45 0.2 0.1 0.6],...
		'XTick',[],...
		'Xlim', [0 1],...
		'Ylim', [0 size(P,1)]);
	xlabel('Coregistering');
	ylabel('volumes completed');
	drawnow;

	% Base starting estimates on solution for previous run

	q = zeros(1,size(dQ,1));
	for k = 2:(size(P,1))
	    C2 = spm_get_space(spm_str_manip(P(k,:), 'd'));

	    spm_smooth(spm_str_manip(P(k,:),'d'),'spm_mov.img',8);
	    V2 = spm_map('spm_mov.img');

	    for i = 1:h;
	  	for j = 1:length(S)
			B      = spm_matrix([-bb(1,1) -bb(1,2) -S(j) 0 0 0 1 1 1]);
			y      = spm_slice_vol(V2,inv(B*inv(C1)*spm_matrix(q)*C2),[M N],Hold);
			Y(([1:n] + (j - 1)*n)) = y(:);
	  	end
		q0 = [-dXdQ Y]\X;
	 	q      = q - q0(1:size(dQ,1))'*dQ;
	    end
	    Q(k,:) = q;

	    spm_unmap(V2); delete spm_mov.img spm_mov.hdr

	    line(...
		'Xdata',[0.5 0.5],...
		'Ydata',[0 k],...
		'LineWidth',16,...
		'Color', [1 0 0]);
	    drawnow;
	end

	spm_unmap(V1); delete spm_ref.img spm_ref.hdr

	% display results
	% translation and rotation over time series
	%---------------------------------------------------------------------------
	figure(3); spm_clf; subplot(3,1,1); axis off

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
	%----------------------------------------------------------------------
	spm_print
end



% Application of Realignment Parameters - uses P & Q
%----------------------------------------------------------------------------
if any(Flags == 'r')

	Hold = 1;
	if (any(Flags == 'S')) Hold = 5; end

	% Get properties of image to realign to
	%----------------------------------------------------------------------------
	p = spm_str_manip(P(1,:), 'd');
	M = spm_get_space(p);
	[DIM VOX SCALE TYPE OFFSET ORIGIN] = spm_hread(p);

	% Get properties of all the images
	%----------------------------------------------------------------------------
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


	% Start progress plot
	%----------------------------------------------------------------------------
	figure(2);
	delete(get(2,'Children'));
	ax = axes('Position', [0.45 0.2 0.1 0.6],...
		'XTick',[],...
		'Xlim', [0 1],...
		'Ylim', [0 DIM(3)]);
	xlabel('Reslicing');
	ylabel('planes completed');
	drawnow;

	Xm = kron(ones(DIM(2),1),[1:DIM(1)]');
	Ym = kron([1:DIM(2)]',ones(DIM(1),1));
	open_mode = 'w';
	if (any(Flags == 'i') | any(Flags == 'c')) Integral = zeros(prod(DIM(1:2)),DIM(3)); end

	start_vol = 1;
	if (any(Flags == 'n') & ~any(Flags == 'i')) start_vol = 2; end

	for j = 1:DIM(3)
		B     = spm_matrix([0 0 -j 0 0 0 1 1 1]);
		if (any(Flags == 'c')) X = zeros(DIM(1)*DIM(2),size(P,1)); end

		% get masks for this plane
		%------------------------------------------------------------------
		if (any(Flags == 'i') | any(Flags == 'k') | any(Flags == 'c'))
			Count = zeros(prod(DIM(1:2)),1);
			for i = 1:size(P,1)
				M0  = reshape(Matrixes(:,i),4,4);
				M1 = inv(B*M0);

				tiny = 5e-2; % From spm_vol_utils.c
				% Find range of slice
				tmp = M1(1,1)*Xm + M1(1,2)*Ym + (M1(1,3)*0 + M1(1,4));
				Mask =        (tmp >= (1-tiny) & tmp <= (Headers(i,1)+tiny));

				tmp = M1(2,1)*Xm + M1(2,2)*Ym + (M1(2,3)*0 + M1(2,4));
				Mask = Mask & (tmp >= (1-tiny) & tmp <= (Headers(i,2)+tiny));

				tmp = M1(3,1)*Xm + M1(3,2)*Ym + (M1(3,3)*0 + M1(3,4));
				Mask = Mask & (tmp >= (1-tiny) & tmp <= (Headers(i,3)+tiny));

				Count = Count + Mask;
			end
			Mask = (Count == size(P,1));
		end

		% get values for this plane
		%------------------------------------------------------------------
		for i = start_vol:size(P,1)
			M0  = reshape(Matrixes(:,i),4,4);
			M1 = inv(B*M0);

			d  = spm_slice_vol(V(:,i),M1,DIM(1:2),Hold);

			if (any(Flags == 'i') | any(Flags == 'c')) Integral(:,j) = Integral(:,j) + d(:); end

			if (~any(Flags == 'c'))
				% don't need to load the whole time series.
				if (i > 1 | ~any(Flags == 'n'))
					d = d/Headers(i,4);
					if any(Flags == 'k') d = d(:).*Mask; end

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

		if (any(Flags == 'i') | any(Flags == 'c')) Integral(:,j) = Integral(:,j)./(Count + eps); end


		if (any(Flags == 'c'))
			% Adjust for motion in Z
			%------------------------------------------------------------------
			Mask = (Count == size(P,1));

			X = X'; % Transpose to reduce paging
			lzm1 = size(P,1)-1;
			for y=1:DIM(2)
				Zmotion  = Matrixes(7,:)*y + Matrixes(11,:)*j + Matrixes(15,:);
				Zmotion  = Zmotion' - Zmotion(1);
				dZmotion = Matrixes(3,:)' - Matrixes(3,1);
				for x = 1:DIM(1)
					xy = x+(y-1)*DIM(1);
					if (Mask(xy))
						Tac = X(:,xy) - Integral(xy,j);
						Zmotion = Zmotion + dZmotion*x;
						A = [Zmotion Zmotion.^2 [0 ; Zmotion(1:lzm1)] [0 ; Zmotion(1:lzm1)].^2];
						X(:,xy) = Tac + Integral(xy,j) - A*((A'*A)\(A'*Tac));
					end
				end
			end
			X = X';

			start_vol = 1;
			if (any(Flags == 'n')) start_vol = 2; end

			for i = start_vol:size(P,1)
				p = spm_str_manip(P(i,:), 'd');

				d = X(:,i)/Headers(i,4);
				if any(Flags == 'k') d = X(:,i).*Mask; end

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

		line('Xdata',[0.5 0.5], 'Ydata',[0 j],...
			'LineWidth',16, 'Color', [1 0 0]);
		drawnow;
	end

	% Write headers and matrixes
	%------------------------------------------------------------------
	if (~any(Flags == 'N'))
		for i = start_vol:size(P,1)
			p = spm_str_manip(P(i,:), 'd');
			q = max([find(p == '/') 0]);
			q = [p(1:q) 'r' p((q + 1):length(p))];
			spm_hwrite(q,DIM,VOX,Headers(i,4),Headers(i,5),0,ORIGIN,'spm - realigned');
			if (any(Flags == 'm')) spm_get_space(q,M); end
		end
	end

	if any(Flags == 'i')
		% Write integral image (16 bit signed)
		%------------------------------------------------------------------
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
%----------------------------------------------------------------------------
if (any(Flags == 'm') & any(Flags == 'e'))
	for k = 1:size(P,1)
		C2 = spm_get_space(spm_str_manip(P(k,:), 'd'));
		spm_get_space(P(k,:), spm_matrix(Q(k,:))*C2);
	end
end


