Part of nonlinear spatial normalisation.
[Alpha,Beta,Var] = spm_brainwarp2(VG,VF,Affine,basX,basY,basZ,T,fwhm)
VG	- Mapped templates (must have same dimensions)
VF	- Mapped object image
Affine	- The affine transformation which maps between the object
	  and template.
basX	- Basis vectors in X. # rows must eq. VG(1)
basY	- Basis vectors in Y. # rows must eq. VG(2)
basZ	- Basis vectors in Z. # rows must eq. VG(3)
T	- The current parameter estimates.
fwhm	- The approximate smoothness of the images.

Alpha	- A*A - where A is the design matrix
Beta	- A*b - where f is the object image
Var	- the approximate chi^2 (corrected for number of resels).
-----------------------------------------------------------------------

The voxels of g1, g2.. are sampled according to the smoothness of the
image (fwhm). The corresponding voxels of f are determined according to
the current parameter estimates and the affine transform.
See "spm_write_sn.m" for more details about how this is done.


-----------------------------------------------------------------------

The design matrix A is generated internally from:

[diag(df/dx)*B diag(df/dy)*B diag(df/dz)*B ...
	diag(g1)*[1 x y z] ...
	diag(g2)*[1 x y z] ...]

where	df/dx, df/dy & df/dz are column vectors containing the gradient
	of image f with respect to displacements in x, y & z
	(in the space of g).

	B is generated from kron(basZ,kron(basY,BasX)). Each column of
	B is a basis image.

	g1, g2.. are template images.

	x, y & z are simply the spatial coordinates of the voxels of f.

	s1, s2.. are the current estimates for the required scaling
	factors. These are derived from T(3*prod(VG(1:3))+1),
	T(3*prod(VG(1:3))+2)...


The vector b contains [(f - diag(g1)*s1 - diag(g1)*x*s2 - ...)].
