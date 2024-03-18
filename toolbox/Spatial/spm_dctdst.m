function varargout = spm_dctdst(varargin)
% Function pointers to various forms of sin and cosine transforms etc
% FORMAT fun = spm_dctdst
% fun - a structure with function handles
%
% Multidimensional transforms
% FORMAT G = fun.function(F)
% where function can be:
% dct2n  - Multidimensional type II discrete cosine transform
% idct2n - Multidimensional inverse type II discrete cosine transform
% dst1n  - Multidimensional type I discrete sin transform
% idst1n - Multidimensional inverse type I discrete sin transform
% dst2n  - Multidimensional type II discrete sin transform
% idst2n - Multidimensional inverse type II discrete sin transform
%
% One dimensional transforms of columns
% FORMAT G = fun.function(F)
% where function can be:
% dct2   - Type II discrete cosine transform
% idct2  - Inverse type II discrete cosine transform
% dst1   - Type I discrete sin transform
% idst1  - Inverse type I discrete sin transform
% dst2   - Type II discrete sin transform
% idst2  - Inverse type II discrete sin transform
%
% FORMAT A = fun.permute2mat(B,dim)
% B      - Multidimensional array
% dim    - Dimension to put into the columns
% A      - Matrix of re-arranged values
%
% FORMAT B = fun.permute2vol(A,dim,d)
% A      - Matrix of re-arranged values
% dim    - Dimension of multidimensional array
%          corresponding with columns of A
% d      - Dimensions of multidimensional array
% B      - Multidimensional array
%
%__________________________________________________________________________
%
% Code works only for real data. Note that it is still a work in progress,
% so is likely to change considerably.
% Some functions remain undocumented for now.
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging

[varargout{1:nargout}] = spm_subfun(localfunctions,varargin{:});
end


%==========================================================================
% For testing
%==========================================================================
function unused
N  = 5;
f0 = [3; -1; 2; -2; 1];
x  = randn(N,1);

% Work with DCT
% Convolve with f
ff0 = [flipud(f0(2:end)); f0];
tmp = convn([flipud(x); x; flipud(x)],ff0,'same');
tmp = tmp(N+(1:N));

X   = dct2(x);
f   = [f0; zeros(N-length(f0),1)];
F   = fft([f; 0; flipud(f(2:end))]);
F   = real(F(1:N)); % DCT-I
FX  = F.*X;
fx  = idct2(FX);
fprintf('DCT-II:\n');
disp([tmp fx])

% Work with DST
% Convolve with f
ff0 = [flipud(f0(2:end)); f0];
tmp = convn([-flipud(x); x; -flipud(x)], ff0,'same');
tmp = tmp(N+(1:N));

X   = dst2(x);
f   = [f0; zeros(N-length(f0),1)];
F   = fft([f; 0; flipud(f(2:end))]);
F   = real(F((1:N)+1)); % ??
FX  = F.*X;
fx  = idst2(FX);
fprintf('DST-II:\n');
disp([tmp fx])

N1 = 64;
N2 = 64;
ux = randn(N1,N2);
uy = randn(N1,N2);

f  = zeros(N1,N2,'single');
f(1,1)   =  4;
f(1,2)   = -1;
f(2,1)   = -1;
f(end,1) = -1;
f(1,end) = -1;

f  = gpuArray(f);
ux = gpuArray(ux);
uy = gpuArray(uy);

Fx = multidim({@dfta1,@dfts1},f);
Fx = (0.1*Fx.^2 + 2*Fx + 0.000001);

Fy = multidim({@dfts1,@dfta1},f);
Fy = (0.1*Fy.^2 + 2*Fy + 0.000001);

vx = multidim({@idst2,@idct2},multidim({@dst2,@dct2},ux)./Fx);
vy = multidim({@idct2,@idst2},multidim({@dct2,@dst2},uy)./Fy);

subplot(4,2,1); imagesc(ux);  axis image; title('u_x');
subplot(4,2,2); imagesc(uy);  axis image; title('u_y');
subplot(4,2,3); imagesc(vx);  axis image; title('v_x');
subplot(4,2,4); imagesc(vy);  axis image; title('v_y');
drawnow

[idx,idy] = ndgrid(1:N1,1:N2);

px = idx+vx;
py = idy+vy;
subplot(2,1,2);
plotdef(px,py)
axis image
end


%==========================================================================
% For testing
%==========================================================================
function plotdef(px,py)
d  = size(px);

hold off
tx = [px(:,1)      px(:,1)];
ty = [zeros(d(1),1) py(:,1)];
plotdef1(tx,ty,'r-','r-');
hold on

tx = [px(:,end) px(:,end)];
ty = [py(:,end) ones(d(1),1)*(d(2)+1)];
plotdef1(tx,ty,'r-','r-');

tx = [zeros(1,d(2)) ; px(1,:)];
ty = [py(1,:) ; py(1,:)];
plotdef1(tx,ty,'r-','r-');

tx = [px(end,:) ; ones(1,d(2))*(d(1)+1)];
ty = [py(end,:) ; py(end,:)];
plotdef1(tx,ty,'r-','r-');

plotdef1(px,py,'b-','k-');

hold off
end


%==========================================================================
% For testing
%==========================================================================
function plotdef1(px,py,cx,cy)
    sk = 1;
    plot(px(1:sk:end,:), py(1:sk:end,:), cx,...
         px(:,1:sk:end)',py(:,1:sk:end)',cy,'LineWidth',1);
end


%==========================================================================
% Multidimensional type II DCT
%==========================================================================
function x = dct2n(varargin)
    x = multidim(@dct2,varargin{:});
end

%==========================================================================
% Multidimensional inverse type II DCT
%==========================================================================
function x = idct2n(varargin)
    x = multidim(@idct2,varargin{:});
end

%==========================================================================
% Multidimensional type I DST
%==========================================================================
function x = dst1n(varargin)
    x = multidim(@dst1,varargin{:});
end

%==========================================================================
% Multidimensional inverse type I DST
%==========================================================================
function x = idst1n(varargin)
    x = multidim(@idst1,varargin{:});
end

%==========================================================================
% Multidimensional type II DST
%==========================================================================
function x = dst2n(varargin)
    x = multidim(@dst2,varargin{:});
end

%==========================================================================
% Multidimensional inverse type II DST
%==========================================================================
function x = idst2n(varargin)
    x = multidim(@idst2,varargin{:});
end

%==========================================================================
% Multidimensional work in progress
%==========================================================================
function x = dftsn(varargin)
    x = multidim(@dfts1,varargin{:});
end

%==========================================================================
% Multidimensional work in progress
%==========================================================================
function x = dftan(varargin)
    x = multidim(@dfta1,varargin{:});
end

%==========================================================================
% Multidimensional work in progress
%==========================================================================
function x = dftpadn(varargin)
    x = real(multidim(@dftpad,varargin{:}));
end

%==========================================================================
% Do multidimensional transforms
%==========================================================================
function x = multidim(fun,x,varargin)
    % FFT in the 2nd dimension of a 3D array is slow
    % so permuting the data first.
    if nargin<3, dims = find(size(x)>1);
    else,        dims = varargin{1};
    end
    d   = size(x);
    if ~iscell(fun)
        for dim=dims
            x = permute2vol(fun(permute2mat(x,dim)),dim,d);
        end
    else
        if length(fun) ~= length(dims)
            error('Incompatible dimensions.');
        end
        for i = 1:length(dims)
            dim  = dims(i);
            funi = fun{i};
            x    = permute2vol(funi(permute2mat(x,dim)),dim,d);
        end
    end
end


%==========================================================================
% Permute a multidimensional array into a matrix
%==========================================================================
function X = permute2mat(X,dim)
    % Convert an N-D array into a matrix, with elements
    % in the dim dimension as columns.
    perm = [dim, 1:(dim-1) (dim+1):ndims(X)];
    d    = size(X);
    dp   = d(perm);
    if d~=1, X = permute(X, perm); end
    X    = reshape(X,[dp(1) prod(dp(2:end))]);
end


%==========================================================================
% Permute a matrix into a multidimensional array
%==========================================================================
function X = permute2vol(X,dim,d)
    % Convert a matrix into an array of size d
    % with columns put into the dim dimension.
    perm = [dim, 1:(dim-1) (dim+1):numel(d)];
    dp   = d(perm);
    X    = reshape(X,dp);
    iperm(perm) = 1:numel(d);
    if dim~= 1, X = permute(X,iperm); end
end


%==========================================================================
% Work in progress
%==========================================================================
function X = dfts1(X)
    % FFT of matrix columns, after padding the middle with
    % zeros.  Assumes X is symmetric, so the results are real.
    X2  = padft(X);
    X   = real(X2(1:size(X,1),:));
end


%==========================================================================
% Work in progress
%==========================================================================
function X = dfta1(X)
    % FFT of matrix columns, after padding the middle with
    % zeros.  Assumes X is anti-symmetric, so the results are imaginary.
    X2  = padft(X);
    X   = X2((1:size(X,1))+1,:);
end


%==========================================================================
% Work in progress
%==========================================================================
function X = dftpad(X)
    % FFT of matrix columns, after padding the middle with
    % zeros.  Assumes X is anti-symmetric, so the results are imaginary.
    X2  = padft(X);
    X   = X2(1:size(X,1),:);
end

%==========================================================================
%
%==========================================================================
function X2 = padft(X,N2)
    N   = size(X,1);
    if nargin<2, N2 = 2*N; end
    ind = [1:ceil(N/2) ((1-floor(N/2)):0)+N2];
    X2  = zeros(N2,size(X,2),'like',X);
    X2(ind,:) = X;
    X2  = fft(X2,[],1);
end


%==========================================================================
% 1D type II DCT
%==========================================================================
function X = dct2(X)
    % Discrete cosine transform (type-II) of matrix columns.
    N   = size(X,1);
    t   = (pi/(2*N))*(0:(N-ones('like',X)))';
   %w   = exp(-1i*t)/sqrt(2*N); w(1) =  w(1)/sqrt(2);
    wr  =  cos(t)/sqrt(2*N);   wr(1) = wr(1)/sqrt(2);
    wi  = -sin(t)/sqrt(2*N);   wi(1) = wi(1)/sqrt(2);
    if rem(N,2)==1
        ind = [1:N N:-1:1];
        X   = X(ind,:);
        X   = fft(X, [], 1);
        X   = X(1:N,:);
    else
        ind = [1:2:N N:-2:2];
        wr  = 2*wr;
        wi  = 2*wi;
        X   = X(ind,:);
        X   = fft(X, [], 1);
    end
    X   = real(X).*wr - imag(X).*wi;
end


%==========================================================================
% 1D type II IDCT
%==========================================================================
function X = idct2(X)
    % Inverse discrete cosine transform (type-II) of matrix columns.
    N = size(X,1);
    t = (pi/(2*N))*(0:(N-ones('like',X)))';
    w = exp(1i*t)*sqrt(2*N);
    if rem(N,2)==1
        % Odd
        ind  = [1:N 1 N:-1:2];
        ind2 = 1:N;
        w(1) = w(1)*sqrt(2);
        w    = [w; 0; -1i*w(2:N)];
        X    = X(ind,:);
        X    = X.*w;
        X    = ifft(X, [], 1, 'symmetric');
        X    = X(ind2,:);
    else
        % Even
        ind2 = [1:(N/2) ; N:-1:(N/2+1)];
        ind2 = ind2(:);
        w(1) = w(1)/sqrt(2);
        X    = ifft(X.*w, [], 1);
        X    = X(ind2,:);
    end
    X = real(X);
end


%==========================================================================
% 1D type I DST
%==========================================================================
function X = dst1(X)
    % Discrete sine transform (Type-I) of matrix columns.
    % Output matches that from MATLAB's PDE toolbox
    N         = size(X,1);
    X2        = zeros((N+1)*2,size(X,2),'like',X);
    in1       = 2:(N+1);
    in2       = (N+1)*2:-1:(N+3);
    X2(in1,:) =  X;
    X2(in2,:) = -X;
    X2        = fft(X2,[],1);
    X         = -imag(X2(2:(N+1),:))/2;
   %X         = -imag(X2([N+1,2:N],:))/2;
end


%==========================================================================
% 1D type I IDST
%==========================================================================
function X = idst1(X)
    % Inverse discrete sine transform (Type-I) of matrix columns.
    N = size(X,1);
    X = dst1(X)*(2/(N+1));
end


%==========================================================================
% 1D type II DST
%==========================================================================
function X = dst2(X)
    % Discrete sine transform (type-II) of matrix columns.
    N    = size(X,1);
    t    = (pi/(2*N))*(1:(N*ones('like',X)))';
    wr   = -sin(t)/sqrt(2*N);
    wi   =  cos(t)/sqrt(2*N);

    ind1 = 1:N;
    ind2 = (2*N):-1:(N+1);
    X2   = zeros(N*2,size(X,2),'like',X);
    X2(ind1,:) =  X;
    X2(ind2,:) = -X;
    X2   = fft(X2,[],1);
    X    = X2(2:(N+1),:);
    X    = real(X).*wr + imag(X).*wi;
end


%==========================================================================
% 1D type II IDST
%==========================================================================
function X = idst2(X)
    % Inverse discrete sine transform (type-II) of matrix columns.
    N    = size(X,1);
    t    = (pi/(2*N))*(1:(N*ones('like',X)))';
    wr   = -sin(t)*sqrt(2*N);
    wi   =  cos(t)*sqrt(2*N);

    ind1 = 2:(N+1);
    ind2 = (2*N):-1:(N+2);
    X2   = zeros(N*2,size(X,2),'like',X);
    X2(ind1, :) = X.*wr + 1i*(X.*wi);
    X2(ind2, :) = conj(X2(2:N,:));
    X2   = ifft(X2,[],1);
    X    = real(X2(1:N,:));
end

