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


V      = spm_map(P0);
v      = V(1:3).*V(4:6);
M      = spm_get_space(P0);
shift0 = inv(spm_matrix(v(1:3)'/2));
shift  = spm_matrix(V(1:3)'/2);
zoom   = spm_matrix([0 0 0 0 0 0 V(4:6)']);


MT0 = zoom;

MT1 = MT0 * shift * spm_matrix([0 0 0 0 pi 0]) * inv(shift);

MS0 = spm_matrix(v([3 2 1])'/2) * spm_matrix([0 0 0 0 pi/2]) ...
	* shift0 * zoom;

MS1 = MS0 * shift * spm_matrix([0 0 0 0 0 pi]) * inv(shift);

MC0 = spm_matrix(v([3 1 2])'/2) * spm_matrix([0 0 0 0 0 pi/2]) ...
	* spm_matrix([0 0 0 pi/2]) * shift0 * zoom;

MC1 = MC0 * shift * spm_matrix([0 0 0 0 0 pi]) * inv(shift);

fprintf('Transverse 1..');
[tra0, tz0] = spm_render_vol(V,MT0,v([1 2]),[thresh 2]);
fprintf('2.. ');
[tra1, tz1] = spm_render_vol(V,MT1,v([1 2]),[thresh 2]);
fprintf('Saggital 1..');
[sag0, sz0] = spm_render_vol(V,MS0,v([3 2]),[thresh 2]);
fprintf('2.. ');
[sag1, sz1] = spm_render_vol(V,MS1,v([3 2]),[thresh 2]);
fprintf('Coronal 1..');
[cor0, cz0] = spm_render_vol(V,MC0,v([3 1]),[thresh 2]);
fprintf('2.. ');
[cor1, cz1] = spm_render_vol(V,MC1,v([3 1]),[thresh 2]);
fprintf('done\n');

Matrixes = [ 'MT0'; 'MT1'; 'MS0'; 'MS1'; 'MC0'; 'MC1'];
Rens     = ['tra0';'tra1';'sag0';'sag1';'cor0';'cor1'];
Depths   = [ 'tz0'; 'tz1'; 'sz0'; 'sz1'; 'cz0'; 'cz1'];

for i=1:6
	eval([Matrixes(i,:) ' = ' Matrixes(i,:) ' *inv(M);']);
	eval([Rens(i,:) ' = ' Rens(i,:) ' * (64/(eps+max(max(' Rens(i,:) '))));']);
end

save render.mat Matrixes Rens Depths MT0 MT1 MS0 MS1 MC0 MC1 tra0 tra1 sag0 sag1 cor0 cor1 tz0 tz1 sz0 sz1 cz0 cz1

