function spm_mip(X,L,V)
% SPM maximum intensity projection
% FORMAT spm_mip(X,L,V);
% V  -  SPM
% L  -  Talairach coordinates
% V  -  {6 x 1} vector of image & voxel sizes [DIM VOX]
%    or {9 x 1} vector of image & voxel sizes, and origin [DIM VOX ORIGIN]
%       (ORIGIN is required for 2D MIPs)
%_______________________________________________________________________
%
% If the data are 2 dimensional [V(3) = 1] the projection is simply an
% image, otherwise:
%
% spm_mip creates and displays a maximum intensity projection of a point
% list of voxel values (X) and their location (L) in three orthogonal
% views of the brain.  It is assumed voxel locations conform to the space
% defined in the atlas of Talairach and Tournoux (1988), otherwise disable
% the GRID.
%
% This routine loads mip (in MIP.mat).  mip is a 360 x 352 matrix with
% contours and grids defining the space of Talairach & Tournoux (1988)
% A default colormap of 64 levels is assumed.
%
% If global XVIEWER is set to a PGM viewer program, images are also
% displayed with that program.
%
%_______________________________________________________________________
% %W% Karl Friston %E%

global GRID XVIEWER TWD

% default GRID value
%-----------------------------------------------------------------------
if isempty(GRID)
	GRID = 0.6; end


% single slice case
%=======================================================================
% NB: Doesn't work for FLIPped data, since ORIGIN is not mirrored
% within image matrix
if V(3) == 1
	if length(V)==6, V=[V; 0;0;0]; end
	M = [ [diag(V(4:6)), -(V(7:9).*V(4:6))]; [zeros(1,3), 1]];
	rcp = inv(M)*[L; ones(1,size(L,2))];
%	d   = rcp(1,:) > 1 & rcp(1,:) < V(1) & rcp(2) > 1 & rcp(2,:) < V(2);
%	mip = full(sparse(rcp(1,d),rcp(2,d),X(d),V(1),V(2)));
	mip = full(sparse(rcp(1,:),rcp(2,:),X,V(1),V(2)));
	%-NB: rot90(X) & axis xy is equivalent to X' & axis ij
	imagesc([V(4)*(1-V(7)):V(4):V(4)*(V(1)-V(7))],...
		[V(5)*(1-V(8)):V(5):V(5)*(V(2)-V(8))],(1 - mip)')
	axis xy image off
	xlabel('x'), ylabel('y')
	return
end

% 3d case
%=======================================================================

% remove negtive values from point list and scale to a maximium of unity
%-----------------------------------------------------------------------
X    = X(:)';
d    = X > 0;
L    = round(L(:,d));
X    = X(d);
X    = X/max(X);

% load mip and create maximum intensity projection
%-----------------------------------------------------------------------
load MIP
mip  = mip96*GRID;
d    = zeros(size(mip));
spm_project(X,L,d,V(1:6));
mip  = max(d,mip);
image(rot90((1 - mip)*64)); axis image; axis off;

%
%
%
if ~isempty(XVIEWER) & isstr(XVIEWER)
        if isempty(TWD); TWD = '/tmp'; end
        % This should really be a subroutine...
        t = clock;
        TmpNm = sprintf('/%s/mip%02d%02d%02d%02d%02d.pgm',TWD, ...
                                  floor(t(3:6)),floor(100*(t(6)-floor(t(6)))));
        fid = fopen(TmpNm,'w');
        fprintf(fid,'P5\n%d %d\n255\n',size(mip));
        fwrite(fid,fliplr(mip*255),'uchar');
        fclose(fid);
        unix(['(' XVIEWER ' ' TmpNm ';\rm ' TmpNm ' ) &']);
end
