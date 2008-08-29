function [U,f1,dJ,iphi] = spm_dartel_variable(U,f,g,param,it1)

    dm = size(U);
    k1 = dm(4);
    dm = dm(1:3);
    param(2:4) = param(2:4)*k1;

    % Generate a sequence of warped templates
    g1    = cell(1,k1);
    g1{1} = g;
    if it1 || k1==1,
        for i=2:k1, g1{i} = g; end
    else
        g = log(max(cat(4,g,1-sum(g,4)),1e-12));
        for i = 1:(k1-1),
            v       = squeeze(single(U(:,:,:,i,:)));
            phi0    = dartel3('Exp',v,[param(8) 1 0]); drawnow
            clear v
            if i==1,
                phi = phi0;
            else
                phi = dartel3('comp',phi0,phi); drawnow
            end
            clear phi0

            g1{i+1} = exp(dartel3('samp',g,phi)); drawnow
            sg      = sum(g1{i+1},4);
            for n=1:size(g1{i+1},4)
                g1{i+1}(:,:,:,n) = g1{i+1}(:,:,:,n)./sg;
            end
            clear sg

        end;
        clear phi
    end

    er2  = 0;
    len2 = 0;
    len1 = 0;
    f1   = f;
    for i = k1:-1:1,
        v         = squeeze(single(U(:,:,:,i,:)));
        [v,ll]    = dartel3(v,f1,g1{i},param); drawnow
        fprintf('%g\t%g\t%g\n', ll);

        er2  = er2 +ll(1)/k1;
        len2 = len2+ll(2)/k1;
        len1 = len1+sqrt(ll(2)/k1);

        [phi0,dJ0]   = dartel3('Exp',v,[param(8) -1 1]); drawnow
        U(:,:,:,i,:) = reshape(v,[dm 1 3]);
        clear v
        if i==k1,
             iphi = phi0;
             dJ   = dJ0;
        else
            [iphi,dJ] = dartel3('comp',phi0,iphi,dJ0,dJ); drawnow
        end
        clear phi0 dJ0
        f1 = dartel3('samp',f,iphi); drawnow
    end;
    fprintf('*** %g %g %g %g ***\n', er2, len2, er2+len2, len1);

