function t = spm_sepmul3d(B1,B2,B3,T)
% Compute t = kron(B3,kron(B2,B1))*T(:)
%
% FORMAT t = spm_sepmul3d(B1,B2,B3,T)
% B1 - x-dim basis functions [nx kx]
% B2 - y-dim basis functions [ny ky]
% B3 - z-dim basis functions [nz kz]
% T  - parameters encoding of the field [kx ky kz]
% t  - Reconstructed field [nx ny nz]
%
% If T is a vector, then so is the output
%
% Note that DCT basis functions are usually used,
% but other forms are available.  For example,
% sparse B-spline basis functions could be used.
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2003-2023 Wellcome Centre for Human Neuroimaging

d  = [size(B1); size(B2); size(B3)];
t  = zeros(d(:,1)','like',T);
t1 = reshape(T, d(1,2)*d(2,2),d(3,2))*B3';
for k=1:d(3,1)
    t(:,:,k) = B1*reshape(t1(:,k),d(1:2,2)')*B2';
end
if  iscolumn(T), t = t(:);
elseif isrow(T), t = t(:)'; end

