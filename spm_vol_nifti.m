function V = spm_vol_nifti(fname,n)
if nargin<2,  n = [1 1];      end;
if ischar(n), n = str2num(n); end;
N = nifti(fname);
n = [n 1 1];
n = n(1:2);
dm = [N.dat.dim 1 1 1 1];
if any(n>dm(4:5)), V = []; return; end;

dt = struct(N.dat);
dt = double([dt.dtype dt.be]);
%d  = find(dt=='-');
%t  = spm_type(lower(dt(1:(d-1))));
%if spm_platform('bigend') ~= strcmp('BE',dt((d+1):end)),
%    t = t*256;
%end;

if isfield(N.extras,'mat') && size(N.extras.mat,3)>=n(1) && sum(sum(N.extras.mat(:,:,n(1))))~=0,
    mat = N.extras.mat(:,:,n(1));
else
    mat = N.mat;
end;

off = (n(1)-1+dm(4)*(n(2)-1))*ceil(spm_type(dt(1),'bits')*dm(1)*dm(2)/8)*dm(3) + N.dat.offset;
V   = struct('fname', N.dat.fname,...
             'mat',   mat,...
             'dim',  dm(1:3),...
             'dt',   dt,...
             'pinfo',[N.dat.scl_slope N.dat.scl_inter off]',...
             'n',    n,...
             'descrip', N.descrip,...
             'private',N);
