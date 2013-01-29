function [def,jac] = spm_get_image_def(P,ds,defa,ddefa)
% Estimation of deformation field (and optionally Jacobian field) 
% for specified image.
%
% FORMAT [def] = spm_get_image_def(i,ds,(defa));
% or
% FORMAT [def] = spm_get_image_def(P,ds,(defa));
% or
% FORMAT [def,jac] = spm_get_image_def(i,ds,(defa),(ddefa));
% or
% FORMAT [def,jac] = spm_get_image_def(P,ds,(defa),(ddefa));
%
%
% i          - Index into array of file handles given in ds.
% P          - File-name or -handle of file that was aquired
%              in same session as the files in ds.P. Note that
%              P does not have to be one of the files used to
%              estimate the partial derivatives of the
%              deformation fields.
% ds         - Structure returned by spm_FindFields.m
% defa       - Array of partial derivative deformation fields
%              sacled to mm^-1 or deg^-1 (or squares thereof).
%              If not provided it will be calculated, but it can
%              be a good idea to calculate it once and for all
%              if one does repeated calls with the same ds.
% ddefa      - Array of partial derivative in the phase encoding
%              direction of partial derivatives (w.r.t. movement
%              parameters) of deformation fields. Used when local
%              Jacobians are to be estimated along with deformation
%              fields.
%
%
% def        - Deformation field for file given by P, or by
%              ds.P(i). Add to xyz when calling spm_sample_vol.m
% jac        - Field of determinants of local Jacobians, i.e.
%              determinants of array of partial derivatives of
%              for dx'/dx, dy'/dy, dy'/dx etc where x', y'... are
%              transformed coordinates and x, y... are original
%              coordinates.
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson
% $Id: spm_get_image_def.m 5219 2013-01-29 17:07:07Z spm $



if nargin > 2 && ~isempty(defa)
   def = zeros(size(defa,1),1);
else 
   def = zeros(prod(ds.P(1).dim(1:3)),1);
end
if nargout > 1
   if nargin > 3 && ~isempty(ddefa)
      jac = ones(size(ddefa,1),1);
   else
      jac = ones(prod(ds.P(1).dim(1:3)),1);
   end
end

%
% Determine weights for linear combination of partial
% deriviatives of deformation field.
%

if isreal(P)   
   if round(P) ~= P
      warning('Index has to be integer');
      return;
   elseif P < 1 || P > size(ds.P,1)
      warning('Index outside range of ds');
      return;
   end
   T = inv(ds.P(P).mat) * ds.P(1).mat;
   q = ds.q(P,:);
else
   if ~isfield(P,'mat'), P = spm_vol(P); end
   T = inv(P.mat) * ds.P(1).mat;
   p = [1 1 1 180/pi 180/pi 180/pi zeros(1,6)] .* spm_imatrix(T);
   q = zeros(1,prod(size(ds.fot))+size(ds.sot,1));
   for i=1:prod(size(ds.fot))
      q(i) = p(ds.fot(i)) - ds.p0(ds.fot(i));
   end
   for j=1:size(ds.sot,1)  % Continue with i
      i = i + 1;
      q(i) = (p(ds.sot(j,1))-ds.p0(ds.sot(j,1))) *...
             (p(ds.sot(j,2))-ds.p0(ds.sot(j,2)));
   end
end

%
% Get deformation field in phase encoding direction.
%
if nargin < 3 || isempty(defa)
   Bx = spm_dctmtx(ds.P(1).dim(1),ds.order(1));
   By = spm_dctmtx(ds.P(1).dim(2),ds.order(2));
   Bz = spm_dctmtx(ds.P(1).dim(3),ds.order(3));
   for i=1:prod(size(ds.fot))+size(ds.sot,1)
      def = def + q(i) * spm_get_def(Bx,By,Bz,ds.beta(:,i));
   end
else
   def = defa * q';
end

%
% Add static field if one was supplied.
%
if isfield(ds,'sfield')
   if ~isempty(ds.sfield)
      def = def + ds.sfield; % Sign change - JAndersson
   end
end

%
% Get Jacobian field in phase encode direction.
%
if nargout > 1
   if nargin <  4 || isempty(ddefa)
      Bx = spm_dctmtx(ds.P(1).dim(1),ds.order(1));
      dBy = spm_dctmtx(ds.P(1).dim(2),ds.order(2),'diff');
      Bz = spm_dctmtx(ds.P(1).dim(3),ds.order(3));
      for i=1:1:prod(size(ds.fot))+size(ds.sot,1)
         jac = jac + q(i) * spm_get_def(Bx,dBy,Bz,ds.beta(:,i));
      end
   else
      jac = jac + ddefa * q';
   end
   %
   % Add Jacobian for static field if one was supplied.
   %
   if isfield(ds,'sjac')
      if ~isempty(ds.sjac)
         jac = jac + ds.sjac;
      end
   end
end
