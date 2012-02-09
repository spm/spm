function [mu,vr,D] = compute_mean(data, param, ord)
d     = data.d;
M_avg = data.mat;
scan  = data.scan;

mu = zeros(d,'single');
if nargout>=3,
    D  = {zeros(d,'single'),zeros(d,'single'),zeros(d,'single')};
end
%fprintf('Computing mean:');
ss0 = zeros(numel(scan),1);
ss1 = zeros(numel(scan),1);
for m=1:d(3),
    mum = zeros(d(1:2),'single');
    mgm = zeros(d(1:2),'single');
    if nargout>=3,
        Dm  = {zeros(d(1:2),'single'),zeros(d(1:2),'single'),zeros(d(1:2),'single')};
    end
    F  = cell(1,numel(scan));
    Dt = cell(1,numel(scan));
    Bf = cell(1,numel(scan));
    Msk= cell(1,numel(scan));
    Dr = cell(3,numel(scan));
    for i=1:numel(scan),
        M  = scan(i).mat\param(i).R*M_avg;

        y = transform_warp(M,param(i).y(:,:,m,:));
        if nargout>=2
            [F{i},Dr{1,i},Dr{2,i},Dr{3,i}]  = shoot3('bsplins',scan(i).f,y,ord);
        else
            F{i}  = shoot3('bsplins',scan(i).f,y,ord);
        end
        Msk{i}  = isfinite(F{i});
        Bf{i}   = exp(param(i).bias(:,:,m));
        Dt{i}   = shoot3('det',param(i).J(:,:,m,:,:));
    end
    mum = zeros(d(1:2),'single');
    mgm = zeros(d(1:2),'single');
    if nargout>=3,
        Dm  = {zeros(d(1:2),'single'),zeros(d(1:2),'single'),zeros(d(1:2),'single')};
    end

    for i=1:numel(scan),
        msk      = Msk{i};
        f        = F{i};
        E        = Bf{i};
        dt       = Dt{i};
        mum(msk) = mum(msk) + f(msk).*E(msk).*dt(msk);
        mgm(msk) = mgm(msk) + E(msk).*E(msk).*dt(msk);

        if nargout>=3
            M       = scan(i).mat\param(i).R*M_avg;
            nm      = sum(msk(:));
            Jm      = reshape(param(i).J(:,:,m,:,:),[d(1)*d(2),3,3]);
            Jm      = reshape(reshape(permute(Jm(msk,:,:),[1 2 3]),nm*3,3)*M(1:3,1:3),[nm 3 3]);
            d1      = Dr{1,i}(msk);
            d2      = Dr{2,i}(msk);
            d3      = Dr{3,i}(msk);

            Dm{1}(msk) = Dm{1}(msk) + (Jm(:,1,1).*d1 + Jm(:,2,1).*d2 + Jm(:,3,1).*d3).*E(msk).*dt(msk);
            Dm{2}(msk) = Dm{2}(msk) + (Jm(:,1,2).*d1 + Jm(:,2,2).*d2 + Jm(:,3,2).*d3).*E(msk).*dt(msk);
            Dm{3}(msk) = Dm{3}(msk) + (Jm(:,1,3).*d1 + Jm(:,2,3).*d2 + Jm(:,3,3).*d3).*E(msk).*dt(msk);

            clear d1 d2 d3
        end
    end
    mu(:,:,m) = mum./(mgm+eps);

    pl = ceil(size(mu,3)/2);
    if m==pl
        vx = sqrt(sum(data.mat(1:3,1:3).^2));
        dm = data.d;
        sc = {(1:dm(1))*vx(1),(1:dm(2))*vx(2)};
        subplot(3,2,1); imagesc(sc{:},mu(:,:,m)'); axis image xy off
        subplot(3,2,2); dt = shoot3('det',param(1).J(:,:,m,:,:)); imagesc(sc{:},dt'); axis image xy off
        subplot(3,2,3); imagesc(sc{:},Bf{1}'); axis image xy off
        subplot(3,2,4); imagesc(sc{:},(F{1}-mu(:,:,m).*Bf{1})'); axis image xy off
        drawnow;
    end
    for i=1:numel(scan),
        msk      = Msk{i};
        f        = F{i}(msk);
        E        = Bf{i}(msk);
        dt       = Dt{i}(msk);
        mum      = mu(:,:,m);
        mum      = mum(msk);
        ss0(i)   = ss0(i) + sum(dt);
        ss1(i)   = ss1(i) + sum((f-mum.*E).^2.*dt);
    end

    if nargout>=3,
        D{1}(:,:,m) = Dm{1}./(mgm+eps);
        D{2}(:,:,m) = Dm{2}./(mgm+eps);
        D{3}(:,:,m) = Dm{3}./(mgm+eps);
    end

end
vr = ss1./ss0;

