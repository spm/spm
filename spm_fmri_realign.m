function spm_fmri_realign(P)
% image realignment routine
% FORMAT spm_fmri_realign(P)
% P - matrix of filenames {one string per row}
%____________________________________________________________________________
%
% spm_fmri_realign realigns a time-series of scans whose filenames
% are in the input matrix (P). The realignment uses explicit and unique
% least squares solutions and partial derivatives of the image process
% with respect to the spatial transformations constituting a rigid-
% body transformation.  These movement parameter estimates are based
% on a first order Taylor expansion in terms of the partial derivatives
% measured using smoothed scans.  These estimates are used to:
%
% i) realign the scans (using sinc interpolation)
% ii) adjust the data to remove movement-related components
%
% The spatially realigned and adjusted images are written to
% the orginal subdirectory with the same filename but prefixed with a 'r'.
% The adjustment procedure is based on a autoregression-moving average
% -like model of the effect of position on signal and explicitly includes
% a spin excitation history effect.
%
% The average of all the realigned scans is written to mean*.img
%
% To avoid artifactual movement-related variance the realigned set of images
% are internally masked, within the set (i.e. if any image has a zero value at
% a voxel than all images have zero values at that voxel).  Zero values
% occur when regions 'outside' the image are moved 'inside' the image during
% realignment.
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
% %W% Karl Friston %E%


%----------------------------------------------------------------------------
global PRINTSTR

% set up file identifiers and check image and voxel sizes are consistent
%----------------------------------------------------------------------------
V     = zeros(12,size(P,1));
for i = 1:size(P,1)
	V(:,i) = spm_map(P(i,:)); end

if ~(all(all(~diff(V([1:6],:)'))))
	error('images are not compatible'); end


% center, bounding box and margins
%----------------------------------------------------------------------------
bb    = spm_bb(P(1,:));					% bb = bounding box
sb    = diff(bb);					% size of bb


% compute matrices
% dQ  = 6 ortholinear transformation parameters
%----------------------------------------------------------------------------
C     = spm_matrix([-mean(bb) 0 0 0 1 1 1]);		% center of bb
D     = spm_matrix([0 0 0 0 0 0 V(4:6,1)']);		% voxel anisotropy 
C     = D*C;						% combine
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
S     = linspace((bb(1,3) + 8/V(6)),(bb(2,3) - 8/V(6)),8);	% transverse
Q     = zeros(size(P,1),size(dQ,1));			% initial estimate
h     = 3;						% number of recursions
M     = sb(1);						% rows per slice
N     = sb(2);						% columns per slice
n     = M*N;						% voxels per slice
FWHM  = [8/V(4) 8/V(5)];
						% in plane smoothing
if V(3,1) == 1; dQ = dQ([1 2 6],:); end			% 3 params for slices
if sb(3) == 0; S = 1; end				% 1 section for slices

% compute X (the reference image) and dX/dQ (the effects of moving X)
%----------------------------------------------------------------------------
X     = zeros(n*length(S),1);
Y     = zeros(n*length(S),1);
dXdQ  = zeros(n*length(S),size(dQ,1));
for i = 1:length(S)
	B     = spm_matrix([-bb(1,1) -bb(1,2) -S(i) 0 0 0 1 1 1]);
	x     = spm_slice_vol(V(:,1),inv(B),[M N],3);
	x     = spm_conv(x,FWHM(1),FWHM(2));
  	X([1:n] + (i - 1)*n) = x(:);
        for j = 1:size(dQ,1);
		A      = spm_matrix(dQ(j,:));
		d      = spm_slice_vol(V(:,1),inv(B*inv(C)*A*C),[M N],3);
		d      = spm_conv(d,FWHM(1),FWHM(2));
		dX     = d - x;
		dXdQ(([1:n] + (i - 1)*n),j) = dX(:);
	end

end


% least squares solution for Q the movements where:
% Y = X + dX/dQ.Q  => Q = [-dX/dQ Y]\X
% The estimate is repeated iteratively h times and the estimates of Q
% are cumulated arithmetically (an aproximation but good enough)
%---------------------------------------------------------------------------
figure(2);clf
axes('Position',[0.45 0.2 0.1 0.6]);
axis([0 1 0 size(P,1)]);
xlabel('estimation of movement parameters')
ylabel('scans completed')
for k = 2:(size(P,1))
    line([0 0],[0 k],'LineWidth',16); drawnow
    for i = 1:h;
  	for j = 1:length(S)
		B      = spm_matrix([-bb(1,1) -bb(1,2) -S(j) 0 0 0 1 1 1]);
		A      = spm_matrix(Q(k,:));
		y      = spm_slice_vol(V(:,k),inv(B*inv(C)*A*C),[M N],3);
		y      = spm_conv(y,FWHM(1),FWHM(2));
		Y(([1:n] + (j - 1)*n)) = y(:);
  	end
 	q      = [-dXdQ Y]\X;
	Q(k,:) = Q(k,:) - q(1:size(dQ,1))'*dQ;
    end
end


% open output file headers and delete existing files (prefixed with a 'r')
%--------------------------------------------------------------------------
U     = P;
for i = 1:size(P,1)
        p      = P(i,:);
        p      = p(p ~= ' ');
        q      = max([find(p == '/') 0]);
        q      = [p(1:q) 'r' p((q + 1):length(p))];
	U(i,[1:length(q)]) = q;
	fid    = fopen(q,'w');
	fclose(fid);
	
	[DIM VOX SCALE TYPE OFFSET ORIGIN] = spm_hread(p);
	spm_hwrite(q,DIM,VOX,SCALE,TYPE,0,ORIGIN,'spm - realigned');
end

% movement-related covariates to be removed
%---------------------------------------------------------------------------
G     = [zeros(1,size(Q,2)); Q([1:(size(Q,1) - 1)],:)];
G     = [Q Q.^2 G G.^2 ones(size(Q,1),1)];


% write the realigned images to output files 
%---------------------------------------------------------------------------
cla
axis([0 1 0 100]);
xlabel('adjusting and writing images')
ylabel('% completed'); drawnow

SSB   = 0;						% sum of squares before
SSA   = 0;						% sum of squares after
for j = 1:DIM(3)
	B     = spm_matrix([0 0 -j 0 0 0 1 1 1]);
	X     = zeros(size(P,1),DIM(1)*DIM(2));

	% get time-series for this plane
	%------------------------------------------------------------------
	for i = 1:size(P,1)
		A  = spm_matrix(Q(i,:));
		d  = spm_slice_vol(V(:,i),inv(B*inv(C)*A*C),DIM(1:2),3);
		X(i,:) = d(:)';
	end

	% remove any components correlated with movement and mask
	%------------------------------------------------------------------
	d     = ~any(~X);
	X     = X(:,d);
	E     = ones(size(X,1),1)*mean(X);
	X     = X - E;
	SSB   = SSB + sum(sum(X.^2));
	X     = X - G*(G\X);
	SSA   = SSA + sum(sum(X.^2));
	X     = X + E;

	% write data
	%------------------------------------------------------------------
	D     = zeros(size(d));
	for i = 1:size(P,1)
        	q   = U(i,:);
        	q   = q(q ~= ' ');
  		fid = fopen(q,'a');
		if ~isempty(X); D(d) = X(i,:)'; end
		fwrite(fid,D/V(7,i),spm_type(TYPE));
		fclose(fid);
	end
        line([0 0],[0 100*j/DIM(3)],'LineWidth',16); drawnow
end

% write mean image
%---------------------------------------------------------------------------
p      = P(1,:);
p      = p(p ~= ' ');
q      = max([find(p == '/') 0]);
CWD    = p(1:(q - 1));
q      = [p(1:q) 'mean' p((q + 1):length(p))];

spm_mean(prod(DIM),TYPE,q,U);
spm_hwrite(q,DIM,VOX,SCALE,TYPE,0,ORIGIN,'spm - mean image');


% unmap and close files
%---------------------------------------------------------------------------
for i = 1:size(P,1); spm_unmap_vol(V(:,i)); end
fclose('all');


% save results
% translation and rotation over time series
%---------------------------------------------------------------------------
eval(['cd ' CWD]);
save REALIGN Q G P SSA SSB

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
figure(2); clf
set(2,'Name','','Pointer','Arrow');
