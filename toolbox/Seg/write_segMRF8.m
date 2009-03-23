function write_segMRF8(p,opts)
% Write out VBM preprocessed data
%____________________________________________________________________________
% Copyright (C) 2008 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: write_segMRF8.m 2921 2009-03-23 17:59:50Z guillaume $


if nargin==1,
    opts = struct('biascor',0,'def',1,'tissues',ones(6,3),'cleanup',0);
end;
if numel(p)>0,
    b0  = spm_load_priors8(p(1).tpm);
end;
for i=1:numel(p),
    preproc_apply(p(i),opts,b0);
end;
return;
%=======================================================================

%=======================================================================
function preproc_apply(p,opts,b0)

sopts = opts.tissues;

[pth,nam,ext]=fileparts(p.image(1).fname);
bsol = p.Tbias;
d    = p.image(1).dim(1:3);

[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3  = 1:d(3);
d3  = [size(bsol{1}) 1];
bB3 = spm_dctmtx(d(3),d3(3),x3);
bB2 = spm_dctmtx(d(2),d3(2),x2(1,:)');
bB1 = spm_dctmtx(d(1),d3(1),x1(:,1));

prm  = [3 3 3 0 0 0];
C    = cell(1,3);
C{1} = spm_bsplinc(p.Twarp(:,:,:,1),prm);
C{2} = spm_bsplinc(p.Twarp(:,:,:,2),prm);
C{3} = spm_bsplinc(p.Twarp(:,:,:,3),prm);


mg  = p.mg;
mn  = p.mn;
vr  = p.vr;
Kb  = max(p.lkp);

for k1=1:size(sopts,1),
    dat{k1} = zeros(d(1:3),'uint8');
    if sopts(k1,3),
        Vt     = struct('fname',   fullfile(pth,['crf', num2str(k1), nam, ext]),...
                        'dim',     p.image(1).dim,...
                        'dt',      [spm_type('uint8') spm_platform('bigend')],...
                        'pinfo',   [1/255 0 0]',...
                        'mat',     p.image(1).mat,...
                        'n',       [1 1],...
                        'descrip', ['Tissue class ' num2str(k1)]);
        Vt     = spm_create_vol(Vt);
        VO(k1) = Vt;
    end;
end;

N   = numel(p.image);
lkp = p.lkp;

M = p.tpm(1).mat\p.Affine*p.image(1).mat;

Phi = (eye(Kb)-1/Kb)*1.0;

for it=1:3,
    spm_progress_bar('init',length(x3),[nam ' (iter ' num2str(it) ')'],'Planes completed');

    for z=1:length(x3),

        % Bias corrected image
        cr = cell(1,N);
        for n=1:N,
            f          = spm_sample_vol(p.image(n),x1,x2,o*x3(z),0);
            cr{n}      = exp(transf(bB1,bB2,bB3(z,:),bsol{n})).*f;
        end

        msk        = (f==0) | ~isfinite(f);
        q1         = likelihoods(cr,[],mg,mn,vr);
        q1         = reshape(q1,[d(1:2),numel(mg)]);
        [t1,t2,t3] = defs(C,z,p.MT,prm,x1,x2,x3,M);
        q          = zeros([d(1:2) Kb]);

        b          = spm_sample_priors8(b0,t1,t2,t3);
        if it>1,
            for k1=1:Kb,
                s{k1} = convn(double(dat{k1}(:,:,z)),[0 1 0; 1 0 1; 0 1 0],'same');
                if z>1,    s{k1} = s{k1} + double(dat{k1}(:,:,z-1)); end
                if z<d(3), s{k1} = s{k1} + double(dat{k1}(:,:,z+1)); end
                s{k1} = s{k1}/(6*255);
            end
            sb = zeros(d(1:2));
            for k1=1:Kb,
                tmp = zeros(d(1:2));
                for k2=1:Kb,
                    tmp = tmp + Phi(k1,k2)*s{k2};
                end
                b{k1} = b{k1}.*exp(tmp);
                sb    = sb + b{k1};
            end
            for k1=1:Kb,
                b{k1} = b{k1}./sb;
            end
        end
        for k1=1:Kb,
            q(:,:,k1) = sum(q1(:,:,lkp==k1),3).*b{k1};
        end

        sq = sum(q,3);
        sw = warning('off','all');
        for k1=1:Kb,
            tmp            = q(:,:,k1);
            tmp(msk)       = 0;
            tmp            = tmp./sq;
            dat{k1}(:,:,z) = uint8(round(255 * tmp));
        end
        warning(sw)
        spm_progress_bar('set',z);
    end
    spm_progress_bar('clear');
end

for z=1:length(x3),
    for k1=1:size(sopts,1),
        tmp = double(dat{k1}(:,:,z))/255;
        spm_write_plane(VO(k1),tmp,z);
    end;
end;

return;
%=======================================================================

%=======================================================================
function [x1,y1,z1] = defs(sol,z,MT,prm,x0,y0,z0,M)
iMT = inv(MT);
x1  = x0*iMT(1,1)+iMT(1,4);
y1  = y0*iMT(2,2)+iMT(2,4);
z1  = (z0(z)*iMT(3,3)+iMT(3,4))*ones(size(x1));
x1a = x0    + spm_bsplins(sol{1},x1,y1,z1,prm);
y1a = y0    + spm_bsplins(sol{2},x1,y1,z1,prm);
z1a = z0(z) + spm_bsplins(sol{3},x1,y1,z1,prm);
x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);
return;
%=======================================================================

%=======================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t = zeros(size(B1,1),size(B2,1),size(B3,1));
end;
return;
%=======================================================================

%=======================================================================
function p = likelihoods(f,bf,mg,mn,vr)
K  = numel(mg);
N  = numel(f);
M  = numel(f{1});
cr = zeros(M,N);
for n=1:N,
    if isempty(bf),
        cr(:,n) = double(f{n}(:));
    else
        cr(:,n) = double(f{n}(:).*bf{n}(:));
    end
end
p  = ones(numel(f{1}),K);
for k=1:K,
    amp    = mg(k)/sqrt((2*pi)^N * det(vr(:,:,k)));
    d      = cr - repmat(mn(:,k)',M,1);
    p(:,k) = amp * exp(-0.5* sum(d.*(d/vr(:,:,k)),2));
end
p = p + 1024*eps;
%=======================================================================

%=======================================================================

