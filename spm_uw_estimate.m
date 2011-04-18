function ds = spm_uw_estimate(P,par)
% Estimation of partial derivatives of EPI deformation fields.
%
% FORMAT [ds] = spm_uw_estimate((P),(par))
%
% P            - List of file names or headers.
% par          - Structure containing parameters governing the specifics
%                of how to estimate the fields.
% .M           - When performing multi-session realignment and Unwarp we
%                want to realign everything to the space of the first
%                image in the first time-series. M defines the space of
%                that.
% .order       - Number of basis functions to use for each dimension.
%                If the third dimension is left out, the order for 
%                that dimension is calculated to yield a roughly
%                equal spatial cut-off in all directions.
%                Default: [12 12 *]
% .sfP         - Static field supplied by the user. It should be a 
%                filename or handle to a voxel-displacement map in
%                the same space as the first EPI image of the time-
%                series. If using the FieldMap toolbox, realignment
%                should (if necessary) have been performed as part of
%                the process of creating the VDM. Note also that the
%                VDM mut be in undistorted space, i.e. if it is
%                calculated from an EPI based field-map sequence
%                it should have been inverted before passing it to
%                spm_uw_estimate. Again, the FieldMap toolbox will
%                do this for you.
% .regorder    - Regularisation of derivative fields is based on the
%                regorder'th (spatial) derivative of the field.
%                Default: 1
% .lambda      - Fudge factor used to decide relative weights of
%                data and regularisation.
%                Default: 1e5
% .jm          - Jacobian Modulation. If set, intensity (Jacobian)
%                deformations are included in the model. If zero,
%                intensity deformations are not considered. 
% .fot         - List of indexes for first order terms to model
%                derivatives for. Order of parameters as defined
%                by spm_imatrix. 
%                Default: [4 5]
% .sot         - List of second order terms to model second 
%                derivatives of. Should be an nx2 matrix where
%                e.g. [4 4; 4 5; 5 5] means that second partial
%                derivatives of rotation around x- and y-axis
%                should be modelled.
%                Default: []
% .fwhm        - FWHM (mm) of smoothing filter applied to images prior
%                to estimation of deformation fields.
%                Default: 6
% .rem         - Re-Estimation of Movement parameters. Set to unity means
%                that movement-parameters should be re-estimated at each
%                iteration.
%                Default: 0
% .noi         - Maximum number of Iterations.
%                Default: 5
% .exp_round   - Point in position space to do Taylor expansion around.
%                'First', 'Last' or 'Average'.
%                Default: 'Average'.
% ds           - The returned structure contains the following fields
% .P           - Copy of P on input.
% .sfP         - Copy of sfP on input (if non-empty).
% .order       - Copy of order on input, or default.
% .regorder    - Copy of regorder on input, or default.
% .lambda      - Copy of lambda on input, or default.
% .fot         - Copy of fot on input, or default.
% .sot         - Copy of sot on input, or default.
% .fwhm        - Copy of fwhm on input, or default.
% .rem         - Copy of rem on input, or default.
% .p0          - Average position vector (three translations in mm
%                and three rotations in degrees) of scans in P.
% .q           - Deviations from mean position vector of modelled
%                effects. Corresponds to deviations (and deviations
%                squared) of a Taylor expansion of deformation fields.
% .beta        - Coeffeicents of DCT basis functions for partial
%                derivatives of deformation fields w.r.t. modelled
%                effects. Scaled such that resulting deformation 
%                fields have units mm^-1 or deg^-1 (and squares 
%                thereof).
% .SS          - Sum of squared errors for each iteration.
%
%
% This is a major rewrite which uses some new ideas to speed up
% the estimation of the field. The time consuming part is the
% evaluation of A'*A where A is a matrix with the partial
% derivatives for each scan with respect to the parameters
% describing the warp-fields. If we denote the derivative
% matrix for a single scan by Ai, then the estimation of A'*A
% is A'*A = A1'*A1 + A2'*A2 + ... +An'*An where n is the number
% of scans in the time-series. If we model the partial-derivative
% fields w.r.t. two movement parameters (e.g. pitch and roll), each
% by [8 8 8] basis-functions and the image dimensions are
% 64x64x64 then each Ai is a 262144x1024 matrix and we need to
% evaluate and add n of these 1024x1024 matrices. It takes a while.
%
% The new idea is based on the realisation that each of these
% matrices is the kroneceker-product of the relevant movement
% parameters, our basis set and a linear combination (given by
% one column of the inverse of the rotation matrix) of the image 
% gradients in the x-, y- and z-directions. This means that 
% they really aren't all that unique, and that the amount of
% information in these matrices doesn't really warrant all 
% those calculations. After a lot of head-scratching I arrived
% at the following
%
% First some definitions
%
% n: no. of voxels
% m: no. of scans
% l: no. of effects to model
% order: [xorder yorder zorder] no. of basis functions
% q: mxl matrix of scaled realignment parameters for effects to model.
% T{i}: inv(inv(P(i).mat)*P(1).mat);
% t(i,:) = T{i}(1:3,2)';
% B: kron(Bz,By,Bx);
% Ax = repmat(dx,1,prod(order)).*B;
% Ay = repmat(dy,1,prod(order)).*B;
% Az = repmat(dz,1,prod(order)).*B;
%
% Now, my hypothesis is that A'*A is given by
%
% AtA = kron((q.*kron(ones(1,l),t(:,1)))'*(q.*kron(ones(1,l),t(:,1))),Ax'*Ax) +...
% kron((q.*kron(ones(1,l),t(:,2)))'*(q.*kron(ones(1,l),t(:,2))),Ay'*Ay) +...
% kron((q.*kron(ones(1,l),t(:,3)))'*(q.*kron(ones(1,l),t(:,3))),Az'*Az) +...
% 2*kron((q.*kron(ones(1,l),t(:,1)))'*(q.*kron(ones(1,l),t(:,2))),Ax'*Ay) +...
% 2*kron((q.*kron(ones(1,l),t(:,1)))'*(q.*kron(ones(1,l),t(:,3))),Ax'*Az) +...
% 2*kron((q.*kron(ones(1,l),t(:,2)))'*(q.*kron(ones(1,l),t(:,3))),Ay'*Az);
%
% Which turns out to be true. This means that regardless of how many
% scans we have we will always be able to create AtA as a sum of
% six individual AtAs. It has been tested against the previous
% implementation and yields exactly identical results.
%
% I know this isn't much of a derivation, but there will be a paper soon. Sorry.
%
% Other things that are new for the rewrite is
% 1. We have removed the possibility to try and estimate the "static" field.
%    There simply isn't enough information about that field, and for a
%    given time series the estimation is just too likely to fail.
% 2. New option to pass a static field (e.g. estimated from a dual
%    echo-time measurement) as an input parameter.
% 3. Re-estimation of the movement parameters at each iteration.
%    Let us say we have a nodding motion in the time series, and
%    at each nod the brain appear to shrink in the phase-encode
%    direction. The realignment will find the nods, but might in
%    addition mistake the "shrinks" for translations, leading to
%    biased estimates. We attempt to correct that by, at iteration,
%    updating both parameters pertaining to distortions and to
%    movement. In order to do so in an unbiased fashion we have
%    also switched from sampling of images (and deformation fields)
%    on a grid centered on the voxel-centers, to a grid whith a 
%    different sampling frequency.   
% 4. Inclusion of a regularisation term. The speed-up has facilitated
%    using more basis-functions, which makes it neccessary to impose
%    regularisation (i.e. punisihing some order derivative of the
%    deformation field).
% 5. Change of interpolation model from tri-linear/Sinc to 
%    tri-linear/B-spline.
% 6. Option to include the Jacobian compression/stretching effects
%    in the model for the estimation of the displacement fields.
%    Our tests have indicated that this is NOT a good idea though. 
%
%_______________________________________________________________________
%
% Definition of some of the variables in this routine.
%
% 
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson
% $Id: spm_uw_estimate.m 4310 2011-04-18 16:07:35Z guillaume $

if nargin < 1 || isempty(P), P = spm_select(Inf,'image'); end
if ~isstruct(P), P = spm_vol(P); end

%
% Hardcoded default input parameters.
%
defpar = struct('sfP',             [],...
                'M',               P(1).mat,...
                'fot',             [4 5],...
                'sot',             [],...
                'hold',            [1 1 1 0 1 0]);

%
% Merge hardcoded defaults and spm_defaults. Translate spm_defaults
% settings to internal defaults.
%
ud = spm_get_defaults('unwarp.estimate');
if isfield(ud,'basfcn'),    defpar.order = ud.basfcn; end
if isfield(ud,'regorder'),  defpar.regorder = ud.regorder; end
if isfield(ud,'regwgt'),    defpar.lambda = ud.regwgt; end
if isfield(ud,'jm'),        defpar.jm = ud.jm; end
if isfield(ud,'fwhm'),      defpar.fwhm = ud.fwhm; end
if isfield(ud,'rem'),       defpar.rem = ud.rem; end
if isfield(ud,'noi'),       defpar.noi = ud.noi; end
if isfield(ud,'expround'),  defpar.exp_round = ud.expround; end

defnames = fieldnames(defpar);

%
% Go through input parameters, chosing the default
% for any parameters that are missing, warning the 
% user if there are "unknown" parameters (probably
% reflecting a misspelling).
%

if nargin < 2 || isempty(par)
   par = defpar;
end
ds = [];
for i=1:length(defnames)
   if isfield(par,defnames{i}) && ~isempty(par.(defnames{i}))
      ds.(defnames{i}) = par.(defnames{i});
   else
      ds.(defnames{i}) = defpar.(defnames{i});
   end
end
parnames = fieldnames(par);
for i=1:length(parnames)
   if ~isfield(defpar,parnames{i})
      warning('Unknown par field %s',parnames{i});
   end
end

%
% Resolve ambiguities.
%
if length(ds.order) == 2
   mm          = sqrt(sum(P(1).mat(1:3,1:3).^2)).*P(1).dim(1:3);
   ds.order(3) = round(ds.order(1)*mm(3)/mm(1));
end
if isfield(ds,'sfP') && ~isempty(ds.sfP)
   if ~isstruct(ds.sfP)
      ds.sfP = spm_vol(ds.sfP);
   end
end
nscan = length(P);
nof   = prod(size(ds.fot)) + size(ds.sot,1);
ds.P  = P;

%
% Get matrix of 'expansion point'-corrected movement parameters 
% for which we seek partial derivative fields.
%
[ds.q,ds.ep] = make_q(P,ds.fot,ds.sot,ds.exp_round);

%
% Create matrix for regularisation of warps.
%
H = make_H(P,ds.order,ds.regorder);

%
% Create a temporary smooth time series to use for
% the estimation.
%
old_P = P;
if ds.fwhm ~= 0
   spm_uw_show('SmoothStart',length(P));
   for i=1:length(old_P)
      spm_uw_show('SmoothUpdate',i);
      sfname(i,:) = [tempname '.img,1,1'];
      to_smooth = sprintf('%s,%d,%d',old_P(i).fname,old_P(i).n);
      spm_smooth(to_smooth,sfname(i,:),ds.fwhm);
   end
   P = spm_vol(sfname);
   spm_uw_show('SmoothEnd');
end

% Now that we have littered the disk with smooth
% temporary files we  should use a try-catch
% block for the rest of the function, to ensure files
% get deleted in the event of an error.

try        % Try block starts here

%
% Initialize some stuff.
%
beta0    = zeros(nof*prod(ds.order),10);
beta     = zeros(nof*prod(ds.order),1);
old_beta = zeros(nof*prod(ds.order),1);

%
% Sample images on irregular grid to avoid biasing
% re-estimated movement parameters with smoothing
% effect from voxel-interpolation (See Andersson ,
% EJNM (1998), 25:575-586.).
%

ss         = [1.1 1.1 0.9];
xs         = 1:ss(1):P(1).dim(1); xs = xs + (P(1).dim(1)-xs(end)) / 2;
ys         = 1:ss(2):P(1).dim(2); ys = ys + (P(1).dim(2)-ys(end)) / 2;
zs         = 1:ss(3):P(1).dim(3); zs = zs + (P(1).dim(3)-zs(end)) / 2;
M          = spm_matrix([-xs(1)/ss(1)+1 -ys(1)/ss(2)+1 -zs(1)/ss(3)+1])*inv(diag([ss 1]))*eye(4);
nx         = length(xs);
ny         = length(ys);
nz         = length(zs);
[x,y,z]    = ndgrid(xs,ys,zs);
Bx         = spm_dctmtx(P(1).dim(1),ds.order(1),x(:,1,1));
By         = spm_dctmtx(P(1).dim(2),ds.order(2),y(1,:,1));
Bz         = spm_dctmtx(P(1).dim(3),ds.order(3),z(1,1,:));
dBy        = spm_dctmtx(P(1).dim(2),ds.order(2),y(1,:,1),'diff');
xyz        = [x(:) y(:) z(:) ones(length(x(:)),1)]; clear x y z;
def        = zeros(size(xyz,1),nof);
ddefa      = zeros(size(xyz,1),nof);

%
% Create file struct for use with spm_orthviews to draw
% representations of the field.
%
dispP       = P(1);
dispP       = rmfield(dispP,{'fname','descrip','n','private'});
dispP.dim   = [nx ny nz];
dispP.dt    = [64 spm_platform('bigend')];
dispP.pinfo = [1 0]';
p          = spm_imatrix(dispP.mat); p = p.*[zeros(1,6) ss 0 0 0];
p(1)       = -mean(1:nx)*p(7);
p(2)       = -mean(1:ny)*p(8);
p(3)       = -mean(1:nz)*p(9);
dispP.mat  = spm_matrix(p); clear p;

%
% We will need to resample the static field (if one was supplied)
% on the same grid (given by xs, ys and zs) as we are going to use
% for the time series. We will assume that the fieldmap has been
% realigned to the space of the first EPI image in the time-series.
%
if isfield(ds,'sfP') && ~isempty(ds.sfP)
   T = ds.sfP.mat\ds.M;
   txyz = xyz*T(1:3,:)';
   c = spm_bsplinc(ds.sfP,ds.hold);
   ds.sfield = spm_bsplins(c,txyz(:,1),txyz(:,2),txyz(:,3),ds.hold);
   ds.sfield = ds.sfield(:);
   clear c txyz;
else
   ds.sfield = [];
end

msk = get_mask(P,xyz,ds,[nx ny nz]);
ssq = [];
%
% Here starts iterative search for deformation fields.
%
for iter=1:ds.noi
   spm_uw_show('NewIter',iter);

   [ref,dx,dy,dz,P,ds,sf,dispP] = make_ref(ds,def,ddefa,P,xyz,M,msk,beta,dispP);
   AtA       = build_AtA(ref,dx,dy,dz,Bx,By,Bz,dBy,ds,P);
   [Aty,yty] = build_Aty(ref,dx,dy,dz,Bx,By,Bz,dBy,ds,P,xyz,beta,sf,def,ddefa,msk);

   % Clean up a bit, cause inverting AtA may use a lot of memory
   clear ref dx dy dz

   % Check that residual error still decreases.
   if iter > 1 && yty > ssq(iter-1)
      %
      % This means previous iteration was no good,
      % and we should go back to old_beta.
      %
      beta = old_beta;
      break;
   else
      ssq(iter) = yty;

      spm_uw_show('StartInv',1);

      % Solve for beta
      Aty = Aty + AtA*beta;
      AtA = AtA + ds.lambda * kron(eye(nof),diag(H)) * ssq(iter)/(nscan*sum(msk));

      try % Fastest if it works
         beta0(:,iter) = AtA\Aty;
      catch % Sometimes necessary
         beta0(:,iter) = pinv(AtA)*Aty;
      end

      old_beta = beta;
      beta     = beta0(:,iter);

      for i=1:nof
         def(:,i)   = spm_get_def(Bx, By,Bz,beta((i-1)*prod(ds.order)+1:i*prod(ds.order)));
         ddefa(:,i) = spm_get_def(Bx,dBy,Bz,beta((i-1)*prod(ds.order)+1:i*prod(ds.order)));
      end

      % If we are re-estimating the movement parameters, remove any DC
      % components from the deformation fields, so that they end up in
      % the movement-parameters instead. It shouldn't make any difference
      % to the variance reduction, but might potentially lead to a better
      % explanation of the variance components.
      % Note that since we sub-sample the images (and the DCT basis set)
      % it is NOT sufficient to reset beta(1) (and beta(prod(order)+1) etc,
      % instead we explicitly subtract the DC component. Note that a DC
      % component does NOT affect the Jacobian.
      %
      if ds.rem ~= 0
         def = def - repmat(mean(def),length(def),1);
      end
      spm_uw_show('EndInv');
      tmp = dispP.mat;
      dispP.mat = P(1).mat;
      spm_uw_show('FinIter',ssq,def,ds.fot,ds.sot,dispP,ds.q);
      dispP.mat = tmp;
   end
   clear AtA
end

ds.P = old_P;
for i=1:length(ds.P)
   ds.P(i).mat = P(i).mat;    % Save P with new movement parameters.
end
ds.beta = reshape(refit(P,dispP,ds,def),prod(ds.order),nof);
ds.SS   = ssq;
if isfield(ds,'sfield');
   ds = rmfield(ds,'sfield');
end

cleanup(P,ds)
spm_uw_show('FinTot');

% Document outcome

spm_print

catch % Try block ends here
   cleanup(P,ds)
   spm_uw_show('FinTot');
   fprintf('procedure terminated abnormally:\n%s',lasterr);
end  % Catch block ends here.

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Utility functions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [q,ep] = make_q(P,fot,sot,exp_round)
%
% P         : Array of file-handles
% fot       : nx1 array of indicies for first order effect.
% sot       : nx2 matrix of indicies of second order effects.
% exp_round : Point in position space to perform Taylor-expansion
%             around ('First','Last' or 'Average'). 'Average' should
%             (in principle) give the best variance reduction. If a
%             field-map acquired before the time-series is supplied
%             then expansion around the 'First' MIGHT give a slightly
%             better average geometric fidelity.

if strcmpi(exp_round,'average');
   %
   % Get geometric mean of all transformation matrices.
   % This will be used as the zero-point in the space 
   % of object position vectors (i.e. the point at which
   % the partial derivatives are estimated). This is presumably
   % a little more correct than taking the average of the
   % parameter estimates. 
   %
   mT = zeros(4);
   for i=1:length(P)
      mT = mT+logm(inv(P(i).mat) * P(1).mat);
   end
   mT = real(expm(mT/length(P)));

elseif strcmpi(exp_round,'first');
   mT = eye(4);
elseif strcmpi(exp_round,'last');
   mT = inv(P(end).mat) * P(1).mat;
else
   warning(sprintf('Unknown expansion point %s',exp_round));
end

%
% Rescaling to degrees makes translations and rotations
% roughly equally scaled. Since the scaling influences
% strongly the effects of regularisation, it would surely
% be nice with a more principled approach. Suggestions
% anyone?
%
ep = [1 1 1 180/pi 180/pi 180/pi zeros(1,6)] .* spm_imatrix(mT);

%
% Now, get a nscan-by-nof matrix containing 
% mean (expansion point) corrected values for each effect we
% model.
%
q = zeros(length(P),prod(size(fot)) + size(sot,1));
for i=1:length(P)
   p = [1 1 1 180/pi 180/pi 180/pi zeros(1,6)] .*...
      spm_imatrix(inv(P(i).mat) * P(1).mat);
   for j=1:prod(size(fot))
      q(i,j) = p(fot(j)) - ep(fot(j));
   end
   for j=1:size(sot,1)
      q(i,prod(size(fot))+j) = (p(sot(j,1)) - ep(sot(j,1))) * (p(sot(j,2)) - ep(sot(j,2)));
   end
end
return
%_______________________________________________________________________

%_______________________________________________________________________
function H = make_H(P,order,regorder)
% Utility Function to create Regularisation term of AtA.

mm = sqrt(sum(P(1).mat(1:3,1:3).^2)).*P(1).dim(1:3);
kx=(pi*((1:order(1))'-1)/mm(1)).^2;
ky=(pi*((1:order(2))'-1)/mm(2)).^2;
kz=(pi*((1:order(3))'-1)/mm(3)).^2;

if regorder == 0
   %
   % Cost function based on sum of squares
   %
   H = (1*kron(kz.^0,kron(ky.^0,kx.^0)) +...
        1*kron(kz.^0,kron(ky.^0,kx.^0)) +...
        1*kron(kz.^0,kron(ky.^0,kx.^0)) );
elseif regorder == 1
   %
   % Cost function based on sum of squared 1st derivatives
   %
   H  = (1*kron(kz.^1,kron(ky.^0,kx.^0)) +...
         1*kron(kz.^0,kron(ky.^1,kx.^0)) +...
         1*kron(kz.^0,kron(ky.^0,kx.^1)) );
elseif regorder == 2
   %
   % Cost function based on sum of squared 2nd derivatives
   %
   H = (1*kron(kz.^2,kron(ky.^0,kx.^0)) +...
        1*kron(kz.^0,kron(ky.^2,kx.^0)) +...
        1*kron(kz.^0,kron(ky.^0,kx.^2)) +...
        3*kron(kz.^1,kron(ky.^1,kx.^0)) +...
        3*kron(kz.^1,kron(ky.^0,kx.^1)) +...
        3*kron(kz.^0,kron(ky.^1,kx.^1)) );
elseif regorder == 3
   %
   % Cost function based on sum of squared 3rd derivatives
   %
   H = (1*kron(kz.^3,kron(ky.^0,kx.^0)) +...
        1*kron(kz.^0,kron(ky.^3,kx.^0)) +...
        1*kron(kz.^0,kron(ky.^0,kx.^3)) +...
        3*kron(kz.^2,kron(ky.^1,kx.^0)) +...
        3*kron(kz.^2,kron(ky.^0,kx.^1)) +...
        3*kron(kz.^1,kron(ky.^2,kx.^0)) +...
        3*kron(kz.^0,kron(ky.^2,kx.^1)) +...
        3*kron(kz.^1,kron(ky.^0,kx.^2)) +...
        3*kron(kz.^0,kron(ky.^1,kx.^2)) +...
        6*kron(kz.^1,kron(ky.^1,kx.^1)) );
else
   error('Invalid order of regularisation');
end
return
%_______________________________________________________________________

%_______________________________________________________________________
function AtA = build_AtA(ref,dx,dy,dz,Bx,By,Bz,dBy,ds,P)
% Now lets build up the design matrix A, or rather A'*A  since
% A itself would be forbiddingly large.

if ds.jm, spm_uw_show('StartAtA',10);
else spm_uw_show('StartAtA',6); end

nof   = prod(size(ds.fot)) + size(ds.sot,1);

AxtAx = uwAtA1(dx.*dx,Bx,By,Bz); spm_uw_show('NewAtA',1);
AytAy = uwAtA1(dy.*dy,Bx,By,Bz); spm_uw_show('NewAtA',2);
AztAz = uwAtA1(dz.*dz,Bx,By,Bz); spm_uw_show('NewAtA',3);
AxtAy = uwAtA1(dx.*dy,Bx,By,Bz); spm_uw_show('NewAtA',4);
AxtAz = uwAtA1(dx.*dz,Bx,By,Bz); spm_uw_show('NewAtA',5);
AytAz = uwAtA1(dy.*dz,Bx,By,Bz); spm_uw_show('NewAtA',6);

if ds.jm
   AjtAj = uwAtA1(ref.*ref,Bx,dBy,Bz); spm_uw_show('NewAtA',7);
   AxtAj = uwAtA2( dx.*ref,Bx, By,Bz,Bx,dBy,Bz); spm_uw_show('NewAtA',8);
   AytAj = uwAtA2( dy.*ref,Bx, By,Bz,Bx,dBy,Bz); spm_uw_show('NewAtA',9);
   AztAj = uwAtA2( dz.*ref,Bx, By,Bz,Bx,dBy,Bz); spm_uw_show('NewAtA',10);
end

R = zeros(length(P),3);
for i=1:length(P)
   tmp    = inv(P(i).mat\P(1).mat);
   R(i,:) = tmp(1:3,2)';
end

tmp  = ones(1,nof);
tmp1 = ds.q.*kron(tmp,R(:,1));
tmp2 = ds.q.*kron(tmp,R(:,2));
tmp3 = ds.q.*kron(tmp,R(:,3));
AtA  =    kron(tmp1'*tmp1,AxtAx) + kron(tmp2'*tmp2,AytAy) + kron(tmp3'*tmp3,AztAz) +...
       2*(kron(tmp1'*tmp2,AxtAy) + kron(tmp1'*tmp3,AxtAz) + kron(tmp2'*tmp3,AytAz));

if ds.jm
   tmp  = [ds.q.*kron(tmp,ones(length(P),1))];
   AtA  = AtA + kron(tmp'*tmp,AjtAj) +...
             2*(kron(tmp1'*tmp,AxtAj) + kron(tmp2'*tmp,AytAj) + kron(tmp3'*tmp,AztAj));
end
spm_uw_show('EndAtA');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function AtA = uwAtA1(y,Bx,By,Bz)
% Calculating off-diagonal block of AtA.

[nx,mx] = size(Bx);
[ny,my] = size(By);
[nz,mz] = size(Bz);
AtA     = zeros(mx*my*mz);
for sl =1:nz
   tmp  = reshape(y((sl-1)*nx*ny+1:sl*nx*ny),nx,ny);
   spm_krutil(Bz(sl,:)'*Bz(sl,:),spm_krutil(tmp,Bx,By,1),AtA);
%    AtA  = AtA + kron(Bz(sl,:)'*Bz(sl,:),spm_krutil(tmp,Bx,By,1));
end
return
%_______________________________________________________________________

%_______________________________________________________________________
function AtA = uwAtA2(y,Bx1,By1,Bz1,Bx2,By2,Bz2)
% Calculating cross-term of diagonal block of AtA
% when A is a sum of the type A1+A2 and where both
% A1 and A2 are possible to express as a kronecker
% product of lesser matrices.

[nx,mx1] = size(Bx1,1); [nx,mx2] = size(Bx2,1);
[ny,my1] = size(By1,1); [ny,my2] = size(By2,1);
[nz,mz1] = size(Bz1,1); [nz,mz2] = size(Bz2,1);
AtA      = zeros(mx1*my1*mz1,mx2*my2*mz2);
for sl =1:nz
   tmp  = reshape(y((sl-1)*nx*ny+1:sl*nx*ny),nx,ny);
   spm_krutil(Bz1(sl,:)'*Bz2(sl,:),spm_krutil(tmp,Bx1,By1,Bx2,By2),AtA);
%   AtA  = AtA + kron(Bz1(sl,:)'*Bz2(sl,:),spm_krutil(tmp,Bx1,By1,Bx2,By2));
end
return
%_______________________________________________________________________

%_______________________________________________________________________
function [Aty,yty] = build_Aty(ref,dx,dy,dz,Bx,By,Bz,dBy,ds,P,xyz,beta,sf,def,ddefa,msk)
% Building Aty.

nof = prod(size(ds.fot)) + size(ds.sot,1);
Aty = zeros(nof*prod(ds.order),1);
yty = 0;

spm_uw_show('StartAty',length(P));
for scan = 1:length(P)
   spm_uw_show('NewAty',scan);
   T    = P(scan).mat\ds.M;
   txyz = xyz*T(1:3,:)';
   if ~(all(beta == 0) && isempty(ds.sfield))
      [idef,jac] = spm_get_image_def(scan,ds,def,ddefa);
      txyz(:,2)  = txyz(:,2)+idef;
   end;
   c    = spm_bsplinc(P(scan),ds.hold);
   y    = sf(scan) * spm_bsplins(c,txyz(:,1),txyz(:,2),txyz(:,3),ds.hold);
   if ds.jm && ~(all(beta == 0) && isempty(ds.sfield))
      y = y .* jac;
   end;

   y_diff       = (ref - y).*msk;
   indx         = find(isnan(y));
   y_diff(indx) = 0;
   iTcol        = inv(T(1:3,1:3));
   tmpAty       = spm_get_def(Bx',By',Bz',([dx dy dz]*iTcol(:,2)).*y_diff);
   if ds.jm ~= 0
      tmpAty = tmpAty + spm_get_def(Bx',dBy',Bz',ref.*y_diff);
   end
   for i=1:nof
      rindx      = (i-1)*prod(ds.order)+1:i*prod(ds.order);
      Aty(rindx) = Aty(rindx) + ds.q(scan,i)*tmpAty;
   end
   yty = yty + y_diff'*y_diff;
end
spm_uw_show('EndAty');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [ref,dx,dy,dz,P,ds,sf,dispP] = make_ref(ds,def,ddefa,P,xyz,M,msk,beta,dispP)
% First of all get the mean of all scans given their
% present deformation fields. When using this as the
% "reference" we will explicitly be minimising variance.
%
% First scan in P is still reference in terms of
% y-direction (in the scanner framework). A single set of
% image gradients is estimated from the mean of all scans,
% transformed into the space of P(1), and used for all scans.
%
spm_uw_show('StartRef',length(P));

rem = ds.rem;
if all(beta==0), rem=0; end;

if rem
   [m_ref,D] = get_refD(ds,def,ddefa,P,xyz,msk);
end

ref = zeros(size(xyz,1),1);
for i=1:length(P)
   spm_uw_show('NewRef',i);
   T    = P(i).mat\ds.M;
   txyz = xyz*T(1:3,:)';
   if ~(all(beta == 0) && isempty(ds.sfield))
      [idef,jac] = spm_get_image_def(i,ds,def,ddefa);
      txyz(:,2)  = txyz(:,2)+idef;
   end;
   c     = spm_bsplinc(P(i),ds.hold);
   f     = spm_bsplins(c,txyz(:,1),txyz(:,2),txyz(:,3),ds.hold);
   if ds.jm ~= 0 && exist('jac')==1
      f = f .* jac;
   end
   indx  = find(~isnan(f));
   sf(i) = 1 / (sum(f(indx).*msk(indx)) / sum(msk(indx)));
   ref   = ref + f;

   if rem
      indx   = find(isnan(f)); f(indx) = 0;
      Dty{i} = (((m_ref - sf(i)*f).*msk)'*D)';
   end
end
ref  = ref / length(P);
ref  = reshape(ref,dispP.dim(1:3));
indx = find(~isnan(ref));

% Scale to roughly 100 mean intensity to ensure
% consistent weighting of regularisation.
gl             = spm_global(ref);
sf             = (100 * (sum(ref(indx).*msk(indx)) / sum(msk(indx))) / gl) * sf;
ref            = (100 / gl) * ref;
dispP.dat      = ref;
c              = spm_bsplinc(ref,ds.hold);
txyz           = xyz*M(1:3,:)';
[ref,dx,dy,dz] = spm_bsplins(c,txyz(:,1),txyz(:,2),txyz(:,3),ds.hold);
ref            = ref.*msk;
dx             =  dx.*msk;
dy             =  dy.*msk;
dz             =  dz.*msk;

% Re-estimate (or rather nudge) movement parameters.
if rem ~= 0
   iDtD = inv(D'*D);
   for i=2:length(P)
      P(i).mat = inv(spm_matrix((iDtD*Dty{i})'))*P(i).mat;
   end
   [ds.q,ds.p0] = make_q(P,ds.fot,ds.sot,ds.exp_round);
end
spm_uw_show('EndRef');
return
%_______________________________________________________________________

%_______________________________________________________________________
function [m_ref,D] = get_refD(ds,def,ddefa,P,xyz,msk)
% Get partials w.r.t. movements from first scan.

[idef,jac] = spm_get_image_def(1,ds,def,ddefa);
T    = P(1).mat\ds.M;
txyz = xyz*T';
c    = spm_bsplinc(P(1),ds.hold);
[m_ref,dx,dy,dz] = spm_bsplins(c,txyz(:,1),...
                   txyz(:,2)+idef,txyz(:,3),ds.hold);
indx    = ~isnan(m_ref);
mref_sf = 1 / (sum(m_ref(indx).*msk(indx)) / sum(msk(indx)));
m_ref(indx) = m_ref(indx) .* mref_sf; m_ref(~indx) = 0;
dx(indx)    =    dx(indx) .* mref_sf; dx(~indx)    = 0;
dy(indx)    =    dy(indx) .* mref_sf; dy(~indx)    = 0;
dz(indx)    =    dz(indx) .* mref_sf; dz(~indx)    = 0;
if ds.jm ~= 0
   m_ref = m_ref .* jac;
   dx    = dx    .* jac;
   dy    = dy    .* jac;
   dz    = dz    .* jac;
end
D = make_D(dx.*msk,dy.*msk,dz.*msk,txyz,P(1).mat);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function D = make_D(dx,dy,dz,xyz,M)
% Utility function that creates a matrix of partial
% derivatives in units of /mm and /radian.
%
% dx,dy,dz   - Partial derivatives (/pixel) wrt x-, y- and z-translation.
% xyz        - xyz-matrix (original positions).
% M          - Current voxel->world matrix.

D = zeros(length(xyz),6);
tiny = 0.0001;
for i=1:6
   p      = [0 0 0 0 0 0 1 1 1 0 0 0];
   p(i)   = p(i)+tiny;
   T      = M\spm_matrix(p)*M;
   dxyz   = xyz*T';
   dxyz   = dxyz(:,1:3)-xyz(:,1:3);
   D(:,i) = sum(dxyz.*[dx dy dz],2)/tiny;
end
return
%_______________________________________________________________________

%_______________________________________________________________________
function cleanup(P,ds)
% Delete temporary smooth files.
%
if ~isempty(ds.fwhm) && ds.fwhm > 0
   for i=1:length(P)
      spm_unlink(P(i).fname);
      [fpath,fname,ext] = fileparts(P(i).fname);
      spm_unlink(fullfile(fpath,[fname '.hdr']));
      spm_unlink(fullfile(fpath,[fname '.mat']));
   end
end
return;
%_______________________________________________________________________

%_______________________________________________________________________
function beta = refit(P,dispP,ds,def)
% We have now estimated the fields on a grid that
% does not coincide with the voxel centers. We must
% now calculate the betas that correspond to a
% a grid centered on the voxel centers.
%
spm_uw_show('StartRefit',1);
rsP     = P(1);
p       = spm_imatrix(rsP.mat);
p       = [zeros(1,6) p(7:9) 0 0 0];
p(1)    = -mean(1:rsP.dim(1))*p(7);
p(2)    = -mean(1:rsP.dim(2))*p(8);
p(3)    = -mean(1:rsP.dim(3))*p(9);
rsP.mat = spm_matrix(p); clear p;
[x,y,z] = ndgrid(1:rsP.dim(1),1:rsP.dim(2),1:rsP.dim(3));
xyz     = [x(:) y(:) z(:) ones(numel(x),1)];
txyz    = ((inv(dispP.mat)*rsP.mat)*xyz')';
Bx      = spm_dctmtx(rsP.dim(1),ds.order(1));
By      = spm_dctmtx(rsP.dim(2),ds.order(2));
Bz      = spm_dctmtx(rsP.dim(3),ds.order(3));
nx      = rsP.dim(1); mx = ds.order(1);
ny      = rsP.dim(2); my = ds.order(2);
nz      = rsP.dim(3); mz = ds.order(3);
nof     = numel(ds.fot) + size(ds.sot,1);
for i=1:nof
   dispP.dat   = reshape(def(:,i),dispP.dim(1:3));
   c           = spm_bsplinc(dispP,ds.hold);
   field       = spm_bsplins(c,txyz(:,1),txyz(:,2),txyz(:,3),ds.hold);
   indx        = isnan(field);
   wgt         = ones(size(field));
   wgt(indx)   = 0;
   field(indx) = 0;
   AtA         = zeros(mx*my*mz);
   Aty         = zeros(mx*my*mz,1);
   for sl = 1:nz
      indx  = (sl-1)*nx*ny+1:sl*nx*ny;
      tmp1  = reshape(field(indx).*wgt(indx),nx,ny);
      tmp2  = reshape(wgt(indx),nx,ny);
      spm_krutil(Bz(sl,:)'*Bz(sl,:),spm_krutil(tmp2,Bx,By,1),AtA);
%      AtA   = AtA + kron(Bz(sl,:)'*Bz(sl,:),spm_krutil(tmp2,Bx,By,1));
      Aty   = Aty + kron(Bz(sl,:)',spm_krutil(tmp1,Bx,By,0));
   end
   beta((i-1)*prod(ds.order)+1:i*prod(ds.order)) = AtA\Aty;
end
spm_uw_show('EndRefit');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function msk = get_mask(P,xyz,ds,dm)
%
% Create a mask to avoid regions where data doesnt exist
% for all scans. This mask is slightly less Q&D than that
% of version 1 of Unwarp. It checks where data exist for
% all scans with present movement parameters and given
% (optionally) the static field. It creates a mask with
% non-zero values only for those voxels, and then does a
% 3D erode to preempt effects of re-estimated movement
% parameters and movement-by-susceptibility effects.
%
spm_uw_show('MaskStart',length(P));
msk = true(length(xyz),1);
for i=1:length(P)
   txyz = xyz * (P(i).mat\ds.M)';
   tmsk = (txyz(:,1)>=1 & txyz(:,1)<=P(1).dim(1) &...
           txyz(:,2)>=1 & txyz(:,2)<=P(1).dim(2) &...
           txyz(:,3)>=1 & txyz(:,3)<=P(1).dim(3));
   msk  = msk & tmsk;
   spm_uw_show('MaskUpdate',i);
end
%
% Include static field in mask estimation
% if one has been supplied.
%
if isfield(ds,'sfP') && ~isempty(ds.sfP)
   txyz = xyz * (ds.sfP.mat\ds.M)';
   tmsk = (txyz(:,1)>=1 & txyz(:,1)<=ds.sfP.dim(1) &...
           txyz(:,2)>=1 & txyz(:,2)<=ds.sfP.dim(2) &...
           txyz(:,3)>=1 & txyz(:,3)<=ds.sfP.dim(3));
   msk  = msk & tmsk;
end
msk = erode_msk(msk,dm);
spm_uw_show('MaskEnd');

% maskP                 = P(1);
% [mypath,myname,myext] = fileparts(maskP.fname);
% maskP.fname           = fullfile(mypath,['mask' myext]);
% maskP.dim             = [nx ny nz];
% maskP.dt              = P(1).dt;
% maskP                 = spm_write_vol(maskP,reshape(msk,[nx ny nz]));
return;
%_______________________________________________________________________

%_______________________________________________________________________
function omsk = erode_msk(msk,dim)
omsk = zeros(dim+[4 4 4]);
omsk(3:end-2,3:end-2,3:end-2) = reshape(msk,dim);
omsk = spm_erode(omsk(:),dim+[4 4 4]);
omsk = reshape(omsk,dim+[4 4 4]);
omsk = omsk(3:end-2,3:end-2,3:end-2);
omsk = omsk(:);
return
%_______________________________________________________________________

