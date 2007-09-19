function [sig,a] = spm_dartel_smooth(t,lam,its,vx)
% A function for smoothing tissue probability maps
% FORMAT sig = spm_dartel_smooth(t,lam,its,vx)
%
%________________________________________________________
% (c) Wellcome Centre for NeuroImaging (2007)

% 
if nargin<4, vx  = [1 1 1]; end;
if nargin<3, its = 16;      end;
if nargin<2, lam = 1;       end;

a   = zeros(size(t),'single');
d   = size(t);
W   = zeros([d(1:3) round((d(4)*(d(4)+1))/2)],'single');
s   = sum(t,4);
for k=1:d(4), t(:,:,:,k) = t(:,:,:,k)./s; end
for i=1:its,
    sig = softmax(a);
    gr  = sig - t;
    for k=1:d(4), gr(:,:,:,k) = gr(:,:,:,k).*s; end
    gr  = gr + optimN('vel2mom',a,[2 vx lam 0 0]);
    jj  = d(4)+1;
    for j1=1:d(4),
        W(:,:,:,j1) = sig(:,:,:,j1).*(1-sig(:,:,:,j1)).*s;
        for j2=1:j1-1,
            W(:,:,:,jj) = -sig(:,:,:,j1).*sig(:,:,:,j2).*s;
            jj = jj+1;
        end
    end
    a   = a - optimN(W,gr,[2 vx lam 0 lam*1e-3 1 1]);
    %a  = a - mean(a(:)); % unstable
    a   = a - sum(sum(sum(sum(a))))/numel(a);
    fprintf('%d  %g\n', i, sum(gr(:).^2));
    drawnow 
end;
sig = softmax(a);
return;

function sig = softmax(a)
sig = zeros(size(a),'single');
for j=1:size(a,3),
    aj = double(squeeze(a(:,:,j,:)));
    aj = min(max(aj-mean(aj(:)),-700),700);
    sj = exp(aj);
    s  = sum(sj,3);
    for i=1:size(a,4),
        sig(:,:,j,i) = single(sj(:,:,i)./s);
    end
end;

