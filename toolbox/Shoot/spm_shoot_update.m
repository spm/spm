function [u0,ll1,ll2,grad_norm] = spm_shoot_update(g,f,u0,phi,dt,prm, bs_args,scale)
% Shooting Of Diffeomorphisms (Spawn Of Dartel).
% FORMAT u0 = spm_shoot_update(g,f,u0,phi,dt,prm, bs_args)
% g        - template
% f        - individual
% u0       - initial velocity
% phi      - deformation
% dt       - Jacobian determinants
% prm      - Parameters of differential operator
% bs_args  - interpolation settings
% scale    - scaling of the update step
%
% u0       - updated initial velocity
% ll1      - matching part of objective function
% ll2      - regularisation part of objective function
% grad_norm - Norm of the 1st derivatives
%
% The easiest way to figure out what this function does is to read the code.
%________________________________________________________
% (c) Wellcome Trust Centre for NeuroImaging (2009)

% John Ashburner
% $Id: spm_shoot_update.m 5829 2014-01-06 20:02:09Z john $

if nargin<8, scale = 1.0; end
scale = max(min(scale,1.0),0.0);

d = [size(g{1}),1,1];
d = d(1:3);

if isempty(u0)
    u0   = zeros([d,3],'single');
end

[ll1,b,A] = mnom_derivs(g,f,phi,dt, bs_args);

m0       = spm_diffeo('vel2mom',u0,prm);
ll2      = 0.5*sum(sum(sum(sum(m0.*u0))));
var1     = sum(sum(sum(sum(b.^2))));
b        = b + m0;
var2     = sum(sum(sum(sum(b.^2))));
grad_norm = sqrt(var2/prod(d));
fprintf('%-10.5g %-10.5g %-10.5g %-10.5g %-10.5g\n',...
                            ll1/prod(d), ll2/prod(d), (ll1+ll2)/prod(d),...
                            var2/(var1+eps), grad_norm);
u0      = u0 - scale*spm_diffeo('fmg',A, b, [prm 3 2]);
clear A b
%=======================================================================

%=======================================================================
function [ll,b,A] = mnom_derivs(g,f,phi,dt, bs_args)
% Compute log-likelihood, first and second derivatives for multi-nomial matching
% FORMAT [ll,b,A] = mnom_derivs(g,f,phi,dt,bs_args)
% g       - cell array of template B-spline coefficients
% f       - cell array of individual subject data
% phi     - deformation field
% dt      - Jacobian determinants
% bs_args - B-spline arguments for sampling g.  Defaults to [2 2 2  1 1 1] if
%           not supplied.
%
% ll      - log-likelihood
% b       - first derivatives
% A       - (approx) second derivatives (Fisher information)
%
%_______________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

if nargin<4,
    bs_args = [2 2 2  1 1 1];
end

d  = [size(g{1}),1,1];
d  = d(1:3);

ll = 0;
b  = zeros([d,3],'single');
A  = zeros([d,6],'single');

[id{1},id{2},id{3}] = ndgrid(1:d(1),1:d(2),ceil(bs_args(3)/2)+1);

for z=1:d(3),
    ind = z+(-ceil(bs_args(3)/2):ceil(bs_args(3)/2));
    if bs_args(6), ind = rem(ind-1+d(3),d(3))+1; end

    f1  = cell(size(f));
    for k=1:numel(g),
        % Note the fudge with the indexing because spm_bsplins only works for double.
        [g1,d1,d2,d3] = spm_bsplins(double(g{k}(:,:,ind)),id{:},bs_args);
        slice(k) = struct('mu',exp(g1),'d1',d1,'d2',d2,'d3',d3);
        if isempty(phi)
            f1{k} = f{k}(:,:,z);
        else
            f1{k} = spm_diffeo('samp',f{k},phi(:,:,z,:)).*dt(:,:,z);
        end
    end
    s = zeros(d(1:2));
    for k=1:numel(g), s           = s + slice(k).mu; end
    for k=1:numel(g), slice(k).mu = slice(k).mu./s;  end

    b(:,:,z,:)  = 0;
    A(:,:,z,:)  = 0;
    for k=1:numel(g),
        tmp      = f1{k}.*log(slice(k).mu);
        ll       = ll - sum(tmp(:));
        tmp      = f1{k} - slice(k).mu.*dt(:,:,z);

        b(:,:,z,1) = b(:,:,z,1) + tmp.*slice(k).d1;
        b(:,:,z,2) = b(:,:,z,2) + tmp.*slice(k).d2;
        b(:,:,z,3) = b(:,:,z,3) + tmp.*slice(k).d3;
        for k1=1:numel(g),
            if k1~=k,
                tmp = -slice(k).mu.*slice(k1).mu.*dt(:,:,z);
            else
                tmp = max(slice(k).mu.*(1-slice(k1).mu),0).*dt(:,:,z);
            end
            A(:,:,z,1) = A(:,:,z,1) + tmp.*slice(k).d1.*slice(k1).d1;
            A(:,:,z,2) = A(:,:,z,2) + tmp.*slice(k).d2.*slice(k1).d2;
            A(:,:,z,3) = A(:,:,z,3) + tmp.*slice(k).d3.*slice(k1).d3;
            A(:,:,z,4) = A(:,:,z,4) + tmp.*slice(k).d1.*slice(k1).d2;
            A(:,:,z,5) = A(:,:,z,5) + tmp.*slice(k).d1.*slice(k1).d3;
            A(:,:,z,6) = A(:,:,z,6) + tmp.*slice(k).d2.*slice(k1).d3;
        end
    end
end

