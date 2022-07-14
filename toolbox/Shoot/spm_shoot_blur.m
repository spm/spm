function [sig,a] = spm_shoot_blur(t,prm,its,sig)
% A function for blurring ("smoothing") tissue probability maps
% FORMAT [sig,a_new] = spm_shoot_blur(t,prm,its,sig)
%     t   - sufficient statistics
%     prm - regularisation parameters (1,1,1, 0.01,0.02,1)
%     its - max no. iterations (12)
%     sig - optional starting estimates (ignored for now)
%
%     sig - "smoothed" average
%     a   - parameters
%
% The core of this procedure is described in:
%     John Ashburner & Karl J. Friston.
%     "Computing Average Shaped Tissue Probability Templates"
%     Neuroimage. 2009 Apr 1;45(2):333-41.
%
% However, there is an additional modification such that the the null space
% of the parameters is rotated out.
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2009-2022 Wellcome Centre for Human Neuroimaging


d   = [size(t),1,1,1];
if nargin<3, its = 16;                         end % Maximum no. iterations
if nargin<2, prm = [1,1,1, 0.01 0.02 1];       end % Default regularisation
rits = [1 1]; % No. cycles and no. relaxation iterations

W    = zeros([d(1:3) round(((d(4)-1)*d(4))/2)],'single'); % 2nd derivatives
gr   = zeros([d(1:3),d(4)-1],'single');                   % 1st derivatives

% Re-organise sufficient statistics to a form that is easier to work with
t    = max(t,eps);
s    = sum(t,4);
for k=1:d(4)
    t(:,:,:,k) = t(:,:,:,k)./s;
end
maxs   = max(s(:)); % Used for scaling the regularisation
%prm(4) = prm(4)+maxs*d(4)*1e-6;

% Only d(4)-1 fields need to be estimated because sum(a,4) = 0.  This matrix
% is used to rotate out the null space
R     = null(ones(1,d(4)));

% Initial starting estimates (if sig is passed)
a     = zeros([d(1:3),d(4)-1],'single');
if false && nargin>=4
    % This is disabled because starting estimates of 0 are much more stable
    % compared to initialising with estimates that might have greater
    % magnitudes than the values at the optimum. If this occurs, then
    % Gauss-Newton becomes very unstable. Note that the Hessian when
    % estimates are all zero is equal to Bohning's bound.
    for z=1:d(3) % Loop over planes
        sz = sig(:,:,z,:);
        sz = min(max(sz,0),1);
        sz(~isfinite(sz)) = 1/d(4);
        sz = squeeze(log(double(sz*(1-d(4)*1e-3)+1e-3)));
        for j1=1:(d(4)-1)
            az = zeros(d(1:2));
            for j2=1:d(4)
                az = az + R(j2,j1)*sz(:,:,j2); % Note the rotation
            end
            a(:,:,z,j1) = az;
        end
        clear sz az
    end
end

for i=1:its

    ll  = 0;
    for z=1:d(3) % Loop over planes

        % Compute softmax and -ve log likelihood for this plane
        [sig,ll1] = sftmax(a(:,:,z,:),R,t(:,:,z,:),s(:,:,z,:));
        sig       = double(reshape(sig,[d(1:2),d(4)]));
        ll        = ll - ll1;

        % Compute first derivatives (d(4)-1) x 1 
        grz = sig - double(reshape(t(:,:,z,:),[d(1:2),d(4)]));
        for j1=1:(d(4)-1)
            gr(:,:,z,j1) = 0;
            for j2=1:d(4)
                gr(:,:,z,j1) = gr(:,:,z,j1) + R(j2,j1)*grz(:,:,j2); % Note the rotation
            end
            gr(:,:,z,j1) = gr(:,:,z,j1).*s(:,:,z);
        end

        % Compute d(4) x d(4) matrix of second derivatives at each voxel.
        % These should be positive definate, but rounding errors may prevent this.
        % Regularisation is included to enforce +ve definateness.
        wz = zeros([d(1:2),d(4),d(4)]);
        for j1=1:d(4)
            wz(:,:,j1,j1) =   (1-sig(:,:,j1)).*sig(:,:,j1).*s(:,:,z);
            for j2=1:(j1-1)
                wz(:,:,j1,j2) = -sig(:,:,j1) .*sig(:,:,j2).*s(:,:,z);
                wz(:,:,j2,j1) = wz(:,:,j1,j2);
            end
        end

        % First step of rotating 2nd derivatives to (d(4)-1) x (d(4)-1)
        % by R'*W*R
        wz1 = zeros([d(1:2),d(4),d(4)-1]);
        for j1=1:d(4)
            for j2=1:(d(4)-1)
                tmp = zeros(d(1:2));
                for j3=1:d(4)
                    tmp = tmp + wz(:,:,j1,j3)*R(j3,j2);
                end
                wz1(:,:,j1,j2) = tmp;
            end
        end

        % Second step of rotating 2nd derivatives to (d(4)-1) x (d(4)-1)
        % by R'*W*R
        wz = zeros([d(1:2),d(4)-1,d(4)-1]);
        for j1=1:(d(4)-1)
            for j2=1:(d(4)-1)
                tmp = zeros(d(1:2));
                for j3=1:d(4)
                    tmp = tmp + R(j3,j1)*wz1(:,:,j3,j2);
                end
                wz(:,:,j1,j2) = tmp;
            end
        end

        % First pull out the diagonal of the 2nd derivs
        for j1=1:d(4)-1
            W(:,:,z,j1) = wz(:,:,j1,j1);% + maxs*sqrt(eps('single'))*d(4)^2;
        end

        % Then pull out the off diagonal parts (note that matrices are symmetric)
        jj = d(4);
        for j1=1:d(4)-1
           for j2=(j1+1):(d(4)-1)
               W(:,:,z,jj) = wz(:,:,j2,j1);
               jj = jj+1;
           end
        end
    end

    % ss1 and ss2 are for examining how close the 1st derivatives are to zero.
    % At convergence, the derivatives from the likelihood term should match those
    % from the prior (regularisation) term.
    ss1 = sum(sum(sum(sum(gr.^2))));
    gr1 = spm_field('vel2mom',a,prm);     % 1st derivative of the prior term
    ll1 = 0.5*sum(sum(sum(sum(gr1.*a)))); % -ve log probability of the prior term
    gr  = gr + gr1;                       % Combine the derivatives of the two terms
    ss2 = sum(sum(sum(sum(gr.^2))));      % This should approach zero at convergence
    mx  = max(max(max(sum(gr.^2,4))));
    reg = double(1e-3*sqrt(mx)*d(4));

    fprintf('%2d %8.4f %8.4f %8.4f %g\n', i, ll/prod(d(1:3)),ll1/prod(d(1:3)), (ll+ll1)/prod(d(1:3)), (ss2)/prod(d(1:3)));

    a   = a - spm_field(W,gr,[prm(1:3) prm(4)+reg prm(5:6) rits]); % Gauss-Newton update

    if ss2/ss1<1e-4, break; end        % Converged?
end

sig = sftmax(a,R);


%==========================================================================

%==========================================================================
function [sig,ll] = sftmax(a,R,t,s)
% Softmax function

d     = [size(a) 1 1 1];
sig   = zeros([d(1:3),d(4)+1],'single');
ll    = 0;

for j=1:size(a,3) % Loop over planes

    % Rotate the null-space back in to the data
    aj  = double(reshape(a(:,:,j,:),[d(1:2),d(4)]));
    sj  = zeros([d(1:2),d(4)+1]);
    for j1=1:d(4)+1
        sj(:,:,j1) = 0;
        for j2=1:d(4)
            sj(:,:,j1) = sj(:,:,j1) + R(j1,j2)*aj(:,:,j2);
        end
    end

    % Compute softmax
    if nargin>=3
        tj = double(reshape(t(:,:,j,:),[d(1:2),d(4)+1]));
        st = double(s(:,:,j));
        ll = ll + sum(sum(sum(sj.*tj,3).*st,2),1);
    end
    mx = max(sj,[],3);
   %sj = sj - mx;
    sj = bsxfun(@minus,sj,mx);
    ej = exp(sj)+eps;
    s  = sum(ej,3);
    if nargin>=3
        ll = ll - sum(sum(st.*(mx+log(s)),2),1);
    end
    for i=1:d(4)+1
        sig(:,:,j,i) = single(ej(:,:,i)./s);
    end
end
