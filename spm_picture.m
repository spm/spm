% rendering of regional effects [SPM{Z}] on drawing of the cortical surface
% FORMAT spm_picture
%____________________________________________________________________________
%
% spm_picture is called by spm_sections_ui and uses variables in working
% memory to load images/pictures of the cortical surface.  Regional
% foci from the selected SPM{Z} are rendered on these drawings.
%
% These are crude renders usually used for didactic purposes only.  The
% voxels are assigned to four parasaggistal blocks (x < -24, 0 > x > -24
% 0 < x x< 24 and x > 24) and are displayed (as a MIP) on the appropriate
% surface.
%__________________________________________________________________________

set(2,'Pointer','Watch');
load Split
colormap default
colormap gray
split = [colormap; ones(64,1)*[1 0 0]];
colormap(split)

%----------------------------------------------------------------------------
load spm_rctx
load spm_lctx

T     = t - min(t);
T     = t*length(colormap)/max(t)/2;
T     = T + length(colormap)/2;
sr    = sr*length(colormap)/300;
sl    = sl*length(colormap)/300;

% right medial
%----------------------------------------------------------------------------
d     = XYZ(1,:) > 0 & XYZ(1,:) < 24;
D     = T(d);
P     = round(XYZ(:,d));
for n = 1:length(D)
      y = P(2,n);
      z = P(3,n);
      q = ones(2,4)*D(n);
      if abs(D(n)) > sr(96 - y,228  +  z)
              sr((96 - y):(97 - y),(228 + z):(231 + z)) = q;
      end
end

% right lateral
%----------------------------------------------------------------------------
d     = XYZ(1,:) >= 24;
D     = T(d);
P     = XYZ(:,d);
for n = 1:length(D)
      y = P(2,n);
      z = P(3,n);
      q = ones(2,4)*D(n);
      if abs(D(n))>sr(134 + y,82 + z)
              sr(134 + y:135 + y,82 + z:85 + z) = q;
      end
end

% left medial
%----------------------------------------------------------------------------
d = XYZ(1,:) <= 0 & XYZ(1,:) > (-24);
D = T(d);
P = XYZ(:,d);

for n = 1:length(D)
      y = P(2,n);
      z = P(3,n);
      q = ones(2,4)*D(n);
      if abs(D(n)) > sl(134 + y,228 + z)
              sl(134 + y:135 + y,228 + z:231 + z) = q;
      end
end

% left lateral
%----------------------------------------------------------------------------
d = XYZ(1,:) <= ( -24);
D = T(d);
P = XYZ(:,d);

for n = 1:length(D)
      y = P(2,n);
      z = P(3,n);
      q = ones(2,4)*D(n);
      if abs(D(n)) > sl(96 - y,82 + z)
              sl(96 - y:97 - y,82 + z:85 + z) = q;
      end
end

% delete previous axis
%----------------------------------------------------------------------------
figure(3)
subplot(2,1,2); delete(gca)
subplot(2,2,4); image(sr');axis xy; axis image; axis off
subplot(2,2,3); image(sl');axis xy; axis image; axis off

% unmap and reset pointer (and x locations is necessary)
%----------------------------------------------------------------------------
set(2,'Pointer','Arrow')

if spm_input('write to disk ?',6,'b','yes|no',[1 0]);
	d = spm_input('filename',6,'s');
	d = [CWD '/' d '.img'];
	A = [sl; sr];
	fid = fopen(d,'w');
	fwrite(fid,A,'uint8');
	fclose(fid);
	spm_hwrite(d,[size(A) 1],[1 1 1],1,spm_type('uint8'),0);
end

