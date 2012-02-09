

%Scans   = spm_select(Inf,'nifti');
%Scans=[
%'/home/john/Documents/spm8/toolbox/Longitudinal/fM00223_057.img'
%pyramid(level).scan(i).f'/home/john/Documents/spm8/toolbox/Longitudinal/fM00223_095.img'];
Scans = [...
'/home/john/Documents/spm8/toolbox/Longitudinal/Ged_serial_MRI/grp2_time1_BYRBR_2007-04-25.nii'
'/home/john/Documents/spm8/toolbox/Longitudinal/Ged_serial_MRI/grp2_time2_BYRBR_2008-03-07.nii'];
Scans = [...
'/home/john/Documents/spm8/toolbox/Longitudinal/Ged_serial_MRI/grp2_time1_STEDA_1998-06-23.nii'
'/home/john/Documents/spm8/toolbox/Longitudinal/Ged_serial_MRI/grp2_time2_STEDA_2000-01-20.nii'];

Scans=['oas1.img';'oas5.img'];
Scans = ['/home/john/Longitudinal/Patrick/timepoint1/c01-0.img'
         '/home/john/Longitudinal/Patrick/timepoint3/c01-6.img'];

Scans = ['/home/john/Longitudinal/Patrick/timepoint1/p01-0.img'
         '/home/john/Longitudinal/Patrick/timepoint3/p01-6.img'];

Scans = ['/home/john/IXI/orig/IXI020-Guys-0700-T1.nii'
         '/home/john/IXI/orig/IXI021-Guys-0703-T1.nii'];
%Scans = strvcat('IXI016-Guys-0697-T1.nii','fl_IXI016-Guys-0697-T1.nii');
%Scans = strvcat('../S/single_subj_T1.nii','../S/fl_single_subj_T1.nii');

%Scans = strvcat('/home/john/Longitudinal/Duezel/av31_T1_prae.img','/home/john/Longitudinal/Duezel/av31_T1_T3.img');


Nii = nifti(Scans);
%Nii(1).dat.scl_slope = 0.2;

prec = 1/mean(noise_estimate(Scans).^2);

bparam = [0 0 1e6];
wparam = [1e-2 1e-2 20 4 40];
sparam = [4 1];
ord    = [3 3 3 0 0 0];

B     = zeros(4,4,6);
B(1,4,1)=1;
B(2,4,2)=1;
B(3,4,3)=1;
B([1,2],[1,2],4)=[0 1;-1 0];
B([3,1],[3,1],5)=[0 1;-1 0];
B([2,3],[2,3],6)=[0 1;-1 0];


d = [0 0 0];
for i=1:numel(Nii),
    dm = [size(Nii(i).dat) 1];
    d  = max(d, dm(1:3));
end
d = prod(d-2)^(1/3);
clear pyramid
pyramid(max(ceil(log2(d)-log2(4)),1)) = struct('d',[1 1 1],'mat',eye(4),'scan',[]);
for i=numel(Nii):-1:1,
    pyramid(1).scan(i).f   = single(Nii(i).dat(:,:,:,1,1));
    pyramid(1).scan(i).mat = Nii(i).mat;
end

for level = 2:numel(pyramid),
    for i=numel(Nii):-1:1,
        pyramid(level).scan(i).f   = shoot3('restrict',pyramid(level-1).scan(i).f);
        pyramid(level).scan(i).f(~isfinite(pyramid(level).scan(i).f)) = 0;
        s1 = [size(pyramid(level-1).scan(i).f) 1];
        s2 = [size(pyramid(level  ).scan(i).f) 1];
        s  = s1(1:3)./s2(1:3);
        pyramid(level).scan(i).mat = pyramid(level-1).scan(i).mat*[diag(s), (1-s(:))*0.5; 0 0 0 1]; %%%%%%%%%%%%%
    end
end

for level=1:numel(pyramid),
    for i=1:numel(Nii)
        pyramid(level).scan(i).f = shoot3('bsplinc',pyramid(level).scan(i).f,ord);
    end
end


Mat0 = cat(3,pyramid(1).scan.mat);
dims = zeros(size(Matrices,3),3);
for i=1:size(dims,1),
    dims(i,:) = Nii(i).dat.dim(1:3);
end
[pyramid(1).mat,pyramid(1).d] = compute_avg_mat(Mat0,dims);

for level=2:numel(pyramid),
    pyramid(level).d   = ceil(pyramid(level-1).d/2);
    s                  = pyramid(level-1).d./pyramid(level).d;
    pyramid(level).mat = pyramid(level-1).mat*[diag(s), (1-s(:))*0.5; 0 0 0 1];
end


LL=[];

for level=numel(pyramid):-1:1,

    scan      = pyramid(level).scan;
    M_avg     = pyramid(level).mat;
    d         = pyramid(level).d;
    vx        = sqrt(sum(pyramid(level).mat(1:3,1:3).^2));

    if level==numel(pyramid),
        clear param
        for i=numel(Nii):-1:1,
            param(i) = struct('bias',zeros(d,'single'),'eb',0,'R',eye(4),'r',zeros(6,1),'s2',1,'vr',1,...
                              'v0',zeros([d 3],'single'),'ev',0,'y',identity(d),...
                              'J',repmat(reshape(eye(3,'single'),[1 1 1 3 3]),[d(1:3),1,1]));
        end
    else
        for i=1:numel(Nii),
            param(i).bias = shoot3('resize',param(i).bias,pyramid(level).d);
            bmom          = optimN_mex('vel2mom', param(i).bias, [vx bparam]);
            param(i).eb   = sum(bmom(:).*param(i).bias(:));

            param(i).v0   = shoot3('resize',param(i).v0,pyramid(level).d);
            for i1=1:3,
                s = pyramid(level).d(i1)/pyramid(level+1).d(i1);
                param(i).v0(:,:,:,i1) = param(i).v0(:,:,:,i1)*s;
            end
            m0          = shoot3('vel2mom',param(i).v0,[vx wparam]);
            param(i).ev = sum(sum(sum(sum(m0.*param(i).v0))));
        end
    end


    for iter=1:2*2^(level-1),


        % Compute deformations from initial velocities
        for i=1:numel(param),
            [param(i).y,param(i).J] = spm_shoot3d(param(i).v0,[vx wparam],sparam);
        end

        %fprintf('\nRigid body:');
        if true,
            [mu,vr,D] = compute_mean(pyramid(level), param, ord);
            ll = 0;
            for i=1:numel(param),
                param(i).vr = vr(i);
                ll          = ll - 0.5*prec*numel(scan(i).f)*param(i).vr - 0.5*param(i).eb - 0.5*param(i).ev;
            end
            fprintf(' %.5e', ll);
            LL = [LL ll/prod(d)];

            for i=1:numel(scan),
                [R,dR]    = spm_dexpm(param(i).r,B);
                M         = scan(i).mat\R*M_avg;
                [x1a,x2a] = ndgrid(1:d(1),1:d(2));

                Hess = zeros(12);
                gra  = zeros(12,1);
                for m=1:d(3)
                    dt   = shoot3('det',param(i).J(:,:,m,:,:));
                    E    = exp(param(i).bias(:,:,m));
                    f    = shoot3('bsplins',scan(i).f,transform_warp(M,param(i).y(:,:,m,:)),ord);
                    b    = f-mu(:,:,m).*E;

                    msk  = isfinite(f);
                    E    = E(msk);
                    b    = b(msk);
                    dt   = dt(msk);
                    x1   = x1a(msk);
                    x2   = x2a(msk);
                    d1   = D{1}(:,:,m);d1 = d1(msk).*E;
                    d2   = D{2}(:,:,m);d2 = d2(msk).*E;
                    d3   = D{3}(:,:,m);d3 = d3(msk).*E;

                    A    = [x1(:).*d1(:) x1(:).*d2(:) x1(:).*d3(:) ...
                            x2(:).*d1(:) x2(:).*d2(:) x2(:).*d3(:) ...
                                 m*d1(:)      m*d2(:)      m*d3(:) ...
                                   d1(:)        d2(:)        d3(:)];

                    Hess = Hess + double(A'*bsxfun(@times,A,dt));
                    gra  = gra  + double(A'*(dt.*b));
                end

                dA = zeros(12,6);
                for m=1:6,
                    tmp     = (R*M_avg)\dR(:,:,m)*M_avg;
                    dA(:,m) = reshape(tmp(1:3,:),12,1);
                end

                Hess       = dA'*Hess*dA*prec;
                gra        = dA'*gra*prec;
                param(i).r = param(i).r - Hess\gra;
                clear Y msk E A f dA
            end
            clear D

            % Mean correct the rigid-body transforms and compute matrices
            r_avg = mean(cat(2,param.r),2);
            for i=1:numel(param),
                param(i).r = param(i).r-r_avg;
                param(i).R = spm_dexpm(param(i).r,B);
            end
        end

        if true,
            [mu,vr] = compute_mean(pyramid(level), param, ord);
            ll = 0;
            for i=1:numel(param),
                param(i).vr = vr(i);
                ll          = ll - 0.5*prec*numel(scan(i).f)*param(i).vr - 0.5*param(i).eb - 0.5*param(i).ev;
            end
            fprintf(' %.5e', ll);
            LL = [LL ll/prod(d)];

            for i=1:numel(scan),
                M    = scan(i).mat\param(i).R*M_avg;
                gra  = zeros(d,'single');
                Hess = zeros(d,'single');

                for m=1:d(3)
                    dt          = shoot3('det',param(i).J(:,:,m,:,:));
                    f           = shoot3('bsplins',scan(i).f,transform_warp(M,param(i).y(:,:,m,:)),ord);
                    msk         = isfinite(f);
                    smu         = mu(:,:,m).*exp(param(i).bias(:,:,m));
                    f(~msk)     = 0;
                    smu(~msk)   = 0;
                    gra(:,:,m)  = smu.*(smu-f).*dt*prec;
                    Hess(:,:,m) = smu.*smu.*dt*prec;
                end

                gra           = gra + optimN_mex('vel2mom', param(i).bias, [vx bparam]);
                param(i).bias = param(i).bias - optimN_mex(Hess,gra,[vx bparam 2 2]); % Gauss-Newton update
                clear msk smu gra Hess f
            end

            % Zero-mean the bias parameter fields
            avg_bias = param(1).bias;
            for i=2:numel(param), avg_bias = avg_bias + param(i).bias; end
            avg_bias = avg_bias/numel(param);
            for i=1:numel(param),
                param(i).bias = param(i).bias - avg_bias;
                bmom          = optimN_mex('vel2mom', param(i).bias, [vx bparam]);
                param(i).eb   = sum(bmom(:).*param(i).bias(:));
            end
            clear avg_bias bmom

        end

        if true,
            % Compute deformations
            [mu,vr,D] = compute_mean(pyramid(level), param, ord);
            ll = 0;
            for i=1:numel(param),
                param(i).vr = vr(i);
                ll          = ll - 0.5*prec*numel(scan(i).f)*param(i).vr - 0.5*param(i).eb - 0.5*param(i).ev;
            end
            fprintf(' %.5e', ll);
            LL = [LL ll/prod(d)];

            for i=1:numel(scan),
                gra  = zeros([d,3],'single');
                Hess = zeros([d,6],'single');
                M    = scan(i).mat\param(i).R*M_avg;

                for m=1:d(3)
                    dt   = shoot3('det',param(i).J(:,:,m,:,:));
                    E    = exp(param(i).bias(:,:,m));
                    f    = shoot3('bsplins',scan(i).f,transform_warp(M,param(i).y(:,:,m,:)),ord);
                    b    = f-mu(:,:,m).*E;

                    msk           = ~isfinite(f);
                    b(msk)        = 0;
                    dt(msk)       = 0;
                    d1            = D{1}(:,:,m).*E;
                    d2            = D{2}(:,:,m).*E;
                    d3            = D{3}(:,:,m).*E;
                    gra(:,:,m,1)  = b.*d1.*dt;
                    gra(:,:,m,2)  = b.*d2.*dt;
                    gra(:,:,m,3)  = b.*d3.*dt;
                    Hess(:,:,m,1) = d1.*d1.*dt;
                    Hess(:,:,m,2) = d2.*d2.*dt;
                    Hess(:,:,m,3) = d3.*d3.*dt;
                    Hess(:,:,m,4) = d1.*d2.*dt;
                    Hess(:,:,m,5) = d1.*d3.*dt;
                    Hess(:,:,m,6) = d2.*d3.*dt;
                end
                Hess        = Hess*prec;
                gra         = gra*prec;

                gra         = gra + shoot3('vel2mom',param(i).v0,[vx wparam]);
                param(i).v0 = param(i).v0 - shoot3('fmg',Hess, gra, [vx wparam 2 2]);
            end
            avg_v0 = param(1).v0;
            for i=2:numel(param),
                avg_v0 = avg_v0 + param(i).v0;
            end
            avg_v0 = avg_v0/numel(param);
            for i=1:numel(param),
                param(i).v0 = param(i).v0 - avg_v0;
                m0          = shoot3('vel2mom',param(i).v0,[vx wparam]);
                param(i).ev = sum(sum(sum(sum(m0.*param(i).v0))));
            end
            clear avg_v0 m0

        end

        subplot(3,1,3); plot(LL,'.-'); drawnow;
        fprintf('\n');
    end

end

for i=1:numel(param),
    [param(i).y,param(i).J] = spm_shoot3d(param(i).v0,[vx wparam],sparam);
end


[pth,nam,ext] = fileparts(Nii(1).dat.fname);

Nii = nifti;
Nii.dat = file_array(['avg_' nam '.nii'],pyramid(1).d,'uint8',0,max(mu(:))/255,0);
Nii.mat = pyramid(1).mat;
Nii.descrip = 'Average';
create(Nii)
Nii.dat(:,:,:) = mu;

d  = (shoot3('det',param(2).J)./shoot3('det',param(1).J));
Nii = nifti;
Nii.dat = file_array(['Jacobian' nam '.nii'],pyramid(1).d,'float32',0,1,0);
Nii.mat = pyramid(1).mat;
Nii.descrip = 'Jacobian';
create(Nii)
Nii.dat(:,:,:) = d;


Nii = nifti;
Nii.dat = file_array(['div_' nam '.nii'],pyramid(1).d,'float32',0,1,0);
Nii.mat = pyramid(1).mat;
Nii.descrip = 'Velocity Divergence';
create(Nii)
Nii.dat(:,:,:) = shoot3('div',param(2).v0);


% 1st derivatives of bias correction objective function
% diff('1/2*(f-mu*exp(b))^2','b')
% gives: -mu*exp(b)*(f - mu*exp(b))
%
% 2nd derivatives
% diff(diff('1/2*(f-mu*exp(b))^2','b'),'b')
% gives: mu^2*exp(2*b) - mu*exp(b)*(f - mu*exp(b))
% Note that the second term is ignored in order to have positive values
% At the solution, we would expect the 1st derivatives to be zero, and
% the second term is simply the 1st derivative.

% Updates for mean
% diff('1/2*(f1-mu*exp(b1))^2 + 1/2*(f2-mu*exp(b2))^2','mu')
% gives: - exp(b1)*(f1 - mu*exp(b1)) - exp(b2)*(f2 - mu*exp(b2))
%
% solve('- exp(b1)*(f1 - mu*exp(b1)) - exp(b2)*(f2 - mu*exp(b2)) = 0','mu')
% gives: (f1*exp(b1) + f2*exp(b2))/(exp(2*b1) + exp(2*b2))

%diff('1/2*(f-mu(phi(x,pw1,pw2))*exp(bias(x,pb1,pb2)))^2','pw1')

%-exp(bias)*D(mu)(phi)*diff(phi, pw1)*(f - exp(bias)*mu(phi))


%exp(2*bias)*diff(phi, pw1)*diff(phi, pw2)*D(mu)(phi)^2 - exp(bias)*(f - mu(phi)*exp(bias))*diff(diff(phi, pw1), pw2)*D(mu)(phi) - exp(bias)*(D@@2)(mu)(phi)*(f - mu(phi)*exp(bias))*diff(phi, pw1)*diff(phi, pw2)

if false
nm=2;
d1=randn(nm,1);
d2=randn(nm,1);
d3=randn(nm,1);
d = [d1 d2 d3];
M=randn(3);
Jm0=randn(nm,3,3);
Jm=Jm0;
Jm      = reshape(reshape(permute(Jm,[1 3 2]),nm*3,3)*M(1:3,1:3),[nm 3 3]);
Dm{1} = (Jm(:,1,1).*d1 + Jm(:,1,2).*d2 + Jm(:,1,3).*d3);
Dm{2} = (Jm(:,2,1).*d1 + Jm(:,2,2).*d2 + Jm(:,2,3).*d3);
Dm{3} = (Jm(:,3,1).*d1 + Jm(:,3,2).*d2 + Jm(:,3,3).*d3);
cat(2,Dm{:})

(M'*squeeze(Jm0(1,:,:)))'*d(1,:)'

vx = sqrt(sum(pyramid(1).mat(1:3,1:3).^2));
dm = pyramid(1).d;
for i=1:16
pl=i*10+40;
subplot(4,4,i);
ld = log(dt(:,:,pl));
image((1:dm(1))*vx(1),(1:dm(2))*vx(2),(ld'+0.1)*64/0.2); axis image xy off
end


end

