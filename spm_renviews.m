function spm_renviews(P0,thresh)
% Produce images for rendering activations to
%
% FORMAT spm_renviews(P0,thresh)
% P0     - filename of image to render.
% thresh - threshold above which image is defined as solid.
%_______________________________________________________________________
%
% Produces a matrix file "render.mat" which contains everything that
% "spm_render" is likely to need.
% The rendering only works for 8 bit images.
% There is no menu interface for this program. It is intended for
% expert users only.
%_______________________________________________________________________
% %W% John Ashburner %E%

V      = spm_vol(P0);
v      = V.dim(1:3).*sqrt(sum(V.mat(1:3,1:3).^2));
M      = V.mat;
shift0 = inv(spm_matrix(v/2));
shift  = spm_matrix(V.dim(1:3)/2);
zoom   = spm_matrix([0 0 0 0 0 0 sqrt(sum(V.mat(1:3,1:3).^2))]);

fprintf('Transverse 1..');
MT0 = zoom;
[ren, dep] = spm_render_vol(V,MT0,v([1 2]),[thresh 2]);
tmp=find(ren<0);ren(tmp)=0;ren=ren/max(ren(:));
rend{1} = struct('M',MT0/V.mat,'ren',ren,'dep',dep);

fprintf('2.. ');
MT1 = MT0 * shift * spm_matrix([0 0 0 0 pi 0]) * inv(shift);
[ren, dep] = spm_render_vol(V,MT1,v([1 2]),[thresh 2]);
tmp=find(ren<0);ren(tmp)=0;ren=ren/max(ren(:));
rend{2} = struct('M',MT1/V.mat,'ren',ren,'dep',dep);


fprintf('Saggital 1..');
MS0 = spm_matrix(v([3 2 1])/2) * spm_matrix([0 0 0 0 pi/2]) ...
	* shift0 * zoom;
[ren, dep] = spm_render_vol(V,MS0,v([3 2]),[thresh 2]);
tmp=find(ren<0);ren(tmp)=0;ren=ren/max(ren(:));
rend{3} = struct('M',MS0/V.mat,'ren',ren,'dep',dep);

fprintf('2.. ');
MS1 = MS0 * shift * spm_matrix([0 0 0 0 0 pi]) * inv(shift);
[ren, dep] = spm_render_vol(V,MS0,v([3 2]),[thresh 2]);
tmp=find(ren<0);ren(tmp)=0;ren=ren/max(ren(:));
rend{4} = struct('M',MS1/V.mat,'ren',ren,'dep',dep);


fprintf('Coronal 1..');
MC0 = spm_matrix(v([3 1 2])/2) * spm_matrix([0 0 0 0 0 pi/2]) ...
	* spm_matrix([0 0 0 pi/2]) * shift0 * zoom;
[ren, dep] = spm_render_vol(V,MC0,v([3 1]),[thresh 2]);
tmp=find(ren<0);ren(tmp)=0;ren=ren/max(ren(:));
rend{5} = struct('M',MC0/V.mat,'ren',ren,'dep',dep);

fprintf('2.. ');
MC1 = MC0 * shift * spm_matrix([0 0 0 0 0 pi]) * inv(shift);
[ren, dep] = spm_render_vol(V,MC1,v([3 1]),[thresh 2]);
tmp=find(ren<0);ren(tmp)=0;ren=ren/max(ren(:));
rend{6} = struct('M',MC1/V.mat,'ren',ren,'dep',dep);

fprintf('saving.. ');
save('render.mat','rend');
fprintf('done\n');
