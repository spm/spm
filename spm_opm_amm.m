function [mfD,Yinds] = spm_opm_amm(S)
% models brain signal and interference as a set of geometrically adaptive
% multipole moments
% FORMAT D = spm_opm_amm(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object                                - Default: no Default
%   S.li             -  internal harmonic order   - Default: 9
%   S.le             -  external harmonic order   - Default: 2
%   S.window        - temporal window size (s)        - 10
%   S.prefix        - prefix to filename          - Default 'm'
%   S.corrLim       - correlation limit          - Default 1
%   S.plotSpheroid  - flag to plot spheroid      - Default 1
% Output:
%   D               - denoised MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright  Tim Tierney


%-Set default values
%--------------------------------------------------------------------------
errorMsg = 'an MEEG object must be supplied.';
if ~isfield(S, 'D'),             error(errorMsg); end
if ~isfield(S, 'li'),            S.li = 9; end
if ~isfield(S, 'le'),            S.le = 2; end
if ~isfield(S, 'corrLim'),       S.corrLim = 1; end
if ~isfield(S, 'window'),        S.window = 10; end
if ~isfield(S, 'skip'),          S.skip = 0; end
if ~isfield(S, 'chunkSize'),     S.chunkSize = 512; end
if ~isfield(S, 'prefix'),        S.prefix = 'm'; end
if ~isfield(S, 'reducerank'),    S.reducerank = 0; end
if ~isfield(S, 'plotSpheroid'),  S.plotSpheroid = 0; end

%-Get design matrix
%--------------------------------------------------------------------------
s = sensors(S.D,'MEG');
if isempty(s)==1
    error('Could not find sensor positions')
end

%-Get usable channels
%--------------------------------------------------------------------------
chaninds = indchantype(S.D,'MEG');
badinds = badchannels(S.D);
usedinds = setdiff(chaninds,badinds);
usedLabs= chanlabels(S.D,usedinds);

%usedLabs = intersect(usedLabs,s.label);
[~,sinds] = spm_match_str(usedLabs,s.label);
usedLabs = s.label(sinds);
Yinds = indchannel(S.D,usedLabs);    
    
%-fit the ellipsoid
%--------------------------------------------------------------------------
v = s.chanpos(sinds,:);
n = s.chanori(sinds,:);
vrange = abs((max(v)-min(v)));
[~,ind]=max(vrange);
if ind==1
[ o, r]=spheroid_fit(v,1);
end

if ind==2
[ o, r ]=spheroid_fit(v,2);
end

if ind==3
[ o, r]=spheroid_fit(v,3);
end

if (ind~=2)
error('Y is not longest axis.... fix please')
end

inside = (v(:,1)-o(1)).^2/r(1)^2+(v(:,2)-o(2)).^2/r(2)^2+(v(:,3)-o(3)).^2/r(3)^2;
c = sum(inside>1);
stepsize = max(r*.005);

while c~=size(v,1)
  rt = r-stepsize;
  inside = (v(:,1)-o(1)).^2/rt(1)^2+(v(:,2)-o(2)).^2/rt(2)^2+(v(:,3)-o(3)).^2/rt(3)^2;
  cc = sum(inside>1);
  if(cc>=c)
    r = r-stepsize;  
    c = cc;
  end 
end

if S.plotSpheroid
  figure()
  plot3(v(:,1),v(:,2),v(:,3),'.k')
  daspect([1,1,1])
  hold on
  [X,Y,Z]=ellipsoid(o(1),o(2),o(3),r(1),r(2),r(3),10);
  plot3(X(:),Y(:),Z(:),'.')
  daspect([1,1,1])
end
%-construct the projectors
%--------------------------------------------------------------------------
a = max(r);
b = min(r);

vtest = double(bsxfun(@minus,v,o'));
external = spm_epharm(vtest,n,a,b,S.le);
inelipse  = spm_ipharm(vtest,n,a,b,S.li);

if(S.reducerank)
  [Q,s] = svd(inelipse,'econ'); 
  Ve = cumsum(diag(s))/sum(diag(s));
  inelipse = Q(:,Ve<S.reducerank);
end

Pout =external*pinv(external);
M = eye(size(external,1))-Pout;
Pin =M*inelipse*pinv(M*inelipse)*M;

fprintf('%-40s: %30s\n','Created Design Matrix',spm('time'));

%-Get Data indices
%--------------------------------------------------------------------------


if (size(Yinds,1)~=size(Pin,1))
    error('data size ~= number of sensors with orientation information');
end
    
%-create output dataset object
%--------------------------------------------------------------------------
fprintf('Creating output dataset\n'); 
outname = fullfile(path(S.D),[S.prefix fname(S.D)]);
mfD = clone(S.D,outname);
mfD.save();


%-Update forward modelling information
%--------------------------------------------------------------------------
%fprintf('%-40s: %30s\n','Updating Sensor Information',spm('time'));
grad = mfD.sensors('MEG');
tmpTra= eye(size(grad.coilori,1));
tmpTra(sinds,sinds)=Pin;
grad.tra                = tmpTra*grad.tra;
grad.balance.previous   = grad.balance.current;

params = ['amm: ' mat2str([S.li,S.le,[o',a,b]],6)];

grad.balance.current    = params;
mfD = sensors(mfD,'MEG',grad);
if isfield(mfD,'inv')
    if isfield(mfD.inv{1},'gainmat')
        fprintf(['Clearing current forward model, please recalculate '...
            'with spm_eeg_lgainmat\n']);
        mfD.inv{1} = rmfield(mfD.inv{1},'gainmat');
    end
    if isfield(mfD.inv{1},'datareg')
        mfD.inv{1}.datareg.sensors = grad;
    end
    if isfield(mfD.inv{1},'forward')
        voltype = mfD.inv{1}.forward.voltype;
        mfD.inv{1}.forward = [];
        mfD.inv{1}.forward.voltype = voltype;
        mfD = spm_eeg_inv_forward(mfD,1);
    end
end
mfD.save();



%-canonical correlations
%--------------------------------------------------------------------------
if(S.skip>0)
  chunks = S.skip*mfD.fsample:S.window*mfD.fsample:size(mfD,2); 
else
  chunks = 1:S.window*mfD.fsample:size(mfD,2);
end

if chunks(end)<size(mfD,2)
    chunks(end+1)=size(mfD,2);
end



for i = 1:(length(chunks)-1)
    inds = chunks(i):(chunks(i+1)-1);
    Ytemp = S.D(:,inds,1);
    Y=Ytemp(Yinds,:);
    inner = Pin*Y;
    Ytemp(Yinds,:)=inner;
    
    if S.corrLim<1
        outer = Pout*Y;
        inter = Y-inner-outer;
        Oinner = orth(inner');
        Ointer = orth(inter');
        C = Oinner'*Ointer;
        [~,Sc,Z] = svd(C);
        noise = Ointer*Z;
        s= diag(Sc);
        noisevec= noise(:,1:sum(s>S.corrLim));
        Beta = noisevec'*inner';
        mod = noisevec*Beta;
        binnew = inner-mod';
        Ytemp(Yinds,:)=binnew;
    end
    mfD(:,inds,:)= Ytemp;
end
mfD.save();

%-Complete
%--------------------------------------------------------------------------


fprintf('%-40s: %30s\n','Completed',spm('time'));

end

function [ o, r] = spheroid_fit( X, ax )
%  [o, r] = ellipsoid_fit( X, ax);
%
% Parameters:
%  X   - Coordinates  n x 3 matrix
%  ax  - numeric indicating longer axis
%
% Output:
% o   - origin
% r   - radii


x =X(:,1);
y =X(:,2);
z =X(:,3);
on = ones(size(x,1),1);
b = x.^2 + y.^2 + z.^2;
if ax==1
    A = [ y.^2 + z.^2 - 2*x.^2, 2*x,2*y,2*z,on];
    beta = pinv(A) *b;
    v(1) = -2 * beta(1) - 1;
    v(2) = beta(1) - 1;
    v(3) = beta(1) - 1;
    v = [ v(1) v(2) v(3) 0 0 0 beta( 2 : 5 )' ];
end

if ax==2
    A = [ x.^2 + z.^2 - 2*y.^2, 2*x,2*y,2*z,on];
    beta = pinv(A)*b;
     v(1) = beta(1)-1;
    v(2) = -2*beta(1)-1;
    v(3) = beta(1)-1;
    v = [ v(1) v(2) v(3) 0 0 0 beta( 2 : 5 )' ];
end

if ax==3
    A = [ x.^2 + y.^2 - 2*z.^2, 2*x,2*y,2*z,on];
    beta = pinv(A) *b; 
    v = beta;
    v(1) = beta(1) - 1;
    v(2) = beta(1) - 1;
    v(3) = -2*beta(1) - 1;
    v = [ v(1) v(2) v(3) 0 0 0 beta( 2 : 5 )' ];
end

A = [ v(1) v(4) v(5) v(7); ...
      v(4) v(2) v(6) v(8); ...
      v(5) v(6) v(3) v(9); ...
      v(7) v(8) v(9) v(10) ];
    

o = -A( 1:3, 1:3 ) \ v( 7:9 )';
T = eye( 4 );
T( 4, 1:3 ) = o';
R = T * A * T';
[ vec, s ] = eig( R( 1:3, 1:3 ) / -R( 4, 4 ) );
r = sqrt( 1 ./ diag( abs( s ) ) );
sgns = sign( diag( s ) );
r = r .* sgns;
r =vec*r;

end
