function [W, Nres] = spm_eeg_wavmtx(n, wname, order)
% creates wavelet matrix
% FORMAT [W, Nres] = spm_eeg_wavmtx(n, wname, order)
%
% n     - dimension of wavelet matrix
% wname - wavelet name, currently only 'db4', i.e. Daubechies order 4
% order - number of different scales
%
% W     - quadratic wavelet matrix
% Nres  - vector of indices for different scale levels in W
%_______________________________________________________________________
%
% The matrix W is generated using 'order' iterations. To obtain an
% orthogonal matrix W, n must be a power of 2. Boundary filters are based
% on periodization implemented with circulant matrix filters.
% The code is partially based on software from [1].
%
% W can be truncated by removing higher scales. For example, remove the two
% highest scales:
% trunc = 2;
% [W, Nres] = spm_eeg_wavmtx(256, 'db4', 8);
% W = W(:, 1:Nres(end-trunc));
%_______________________________________________________________________
%
% [1] C Taswell and KC McGill, 1994. Algorithm 735: Wavelet Transform
% Algorithms for Finite-Duration Discrete-Time Signals. ACM Transactions
% on Mathematical Software, (20): 398 - 412. 
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_wavmtx.m 112 2005-05-04 18:20:52Z john $


switch wname
	case 'db4'
		w = 1/(4*sqrt(2))*[(1 + sqrt(3)) (3 + sqrt(3)) (3 - sqrt(3)) (1 - sqrt(3))];
	    c1 = [w(1) w(2); w(4) -w(3)];
	    c2 = [w(3) w(4); w(2) -w(1)];
        
	otherwise
		disp('Unknown wavelet name');

end

lw = length(w);

W = speye(n);
m = n;
Nres = [];

for i = 1:order
    
    k = lw - 2;
    cm2 = ceil(m/2);
    
	% generate transformation matrix C
	C1 = [kron(speye(cm2), c1) sparse(2*cm2, k)];
	C2 = kron(spdiags(ones(cm2+1,1), 1, cm2, cm2+1), c2);
	C = C1 + C2;
    
    k = k/2;
    r = 2*cm2;
    C(:, k+1:2*k) = C(:,k+1:2*k) + C(:,r+k+1:r+2*k);
    C(:, r+1:r+k) = C(:,r+1:r+k) + C(:,1:k);   
    C = C(:, k+1:k+m);
    
    Nres(end+1) = size(C,1)/2;
    
    C = blkdiag(C, speye(size(W,1) - size(C,2)));
    
	% generate permutation matrix P
	P = sparse([1:cm2 cm2+1:2*cm2], [1:2:2*cm2 2:2:2*cm2], 1, 2*cm2,2*cm2);
    P = blkdiag(P, speye(size(C,1) - size(P,2)));
	W = P*C*W;
    
    m = cm2;

end

W = W';

% determine index brackets of scale levels
Nres = [0 fliplr(Nres) size(W,2)];
