function def = spm_get_def(Bx,By,Bz,beta)
% Calculating deformation field from coefficient vector.
%
% FORMAT [def] = get_def(Bx,By,Bz,beta);
% Bx, By, Bz  - Separable basis sets such that B=kron(Bz,kron(By,Bx));
% beta        - Coefficient vector for basis set.
% or
% FORMAT [def] = get_def(dim,order,beta);
% dim         - Dimensionality of image in x,y and z directions.
% order       - Order of DCT set in x,y and z directions.
% beta        - Coefficient vector for basis set.
%_______________________________________________________________________
%
% Calculates the equivalent to kron(Bz,kron(By,Bx))*beta
% in a faster and more efficient way. Note that this routine
% can be used to calculate AtY efficiently as well since 
% A'*y = (diag(dfdy)*kron(Bz,kron(By,Bx)))'*y = 
% kron(Bz',kron(By',Bx'))*(diag(dfdy)*y) = get_def(Bx',By',Bz',dfdy.*y)
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson
% $Id: spm_get_def.m 1143 2008-02-07 19:33:33Z spm $


if nargin == 4
   [nx,mx] = size(Bx);
   [ny,my] = size(By);
   [nz,mz] = size(Bz);
   if mx*my*mz ~= size(beta,1)
      warning('get_def: Size mismatch between beta and basis-set');
      return
   end
elseif nargin == 3
   if numel(Bx) ~= 3 || numel(By) ~= 3
      warning('get_def: Wrong dimensionality on input');
   elseif prod(By) ~= size(Bz,1)
      warning('get_def: Size mismatch between beta and basis-set');
      return
   end
   nx = Bx(1); ny = Bx(2); nz = Bx(3);
   mx = By(1); my = By(2); mz = By(3);
   beta = Bz;
   Bx = spm_dctmtx(nx,mx);
   By = spm_dctmtx(ny,my);
   Bz = spm_dctmtx(nz,mz);
end

def = zeros(nx*ny*nz,1);
for bf = 1:mz
   mbeta = reshape(beta((bf-1)*my*mx+1:bf*my*mx),mx,my);
   tmp = reshape(Bx*mbeta*By',nx*ny,1);
   for sl = 1:nz
      def((sl-1)*ny*nx+1:sl*ny*nx) = def((sl-1)*ny*nx+1:sl*ny*nx) + Bz(sl,bf)*tmp;
   end
end

return

