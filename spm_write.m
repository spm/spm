function spm_write(P,Q,A) 
% Applies a transformation matrix to an image
% FORMAT spm_write(P,Q,A);
% P    - image filename
% Q    - output filename
% A    - transformation matrix
%___________________________________________________________________________
%
% Applies transformation matrix A to an image P and writes the result to Q
% The matrix A is applied after centering the image and correcting for voxel
% anisotropy
%
%__________________________________________________________________________
% %W% %E%


% remove spaces from filename strings
%---------------------------------------------------------------------------
P     = P(P ~= ' ');
Q     = Q(Q ~= ' ');

% read header, memory map and compute A after centering and anisotropic
% correction
%---------------------------------------------------------------------------
[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(P);
V     = spm_map_vol(P,[DIM VOX 1 TYPE OFFSET]);			% memory map
C     = spm_matrix(-DIM/2);					% image center
I     = spm_matrix([0 0 0 0 0 0 VOX]);				% isotropy
A     = inv(C)*inv(I)*A*I*C;					% transformation

% open output file and write transformed image
%---------------------------------------------------------------------------
fid   = fopen(Q,'w');
for i = 1:V(3)
   B  = spm_matrix([0 0 -i 0 0 0 1 1 1]);
   x  = spm_slice_vol(V,inv(B*A),[V(1) V(2)],1);
   fwrite(fid,x,spm_type(TYPE));
end

% close files, unmap and write header
%---------------------------------------------------------------------------
fclose(fid);
spm_hwrite(Q,DIM,VOX,SCALE,TYPE,0,ORIGIN,DESCRIP);
spm_unmap_vol(V);

