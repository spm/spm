function spm_spy(X,Markersize,m)
% Pretty version of spy
% FORMAT spm_spy(X,Markersize,m)
% X    - sparse {m x n} matrix
%
% See also: spy
%__________________________________________________________________________
% Copyright (C) 1994-2016 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_spy.m 6894 2016-09-30 16:48:46Z spm $


% defaults
%--------------------------------------------------------------------------
if nargin < 1, X = defaultspy; end
if nargin < 2, Markersize = 16; end
if nargin < 3, m = max(max(X)); end


spy(X > m/256,Markersize,'.c'), hold on
spy(X > m/2,Markersize,'.b'), hold on
spy(X > (m/2 + m/4),Markersize,'.k'), hold off
axis normal

function X = defaultspy
X = [49 22 76 73 58 75 76 62 1 1 76 78 70 1 66 70 75 62 58 61 1 0 65 77 ...
    77 73 19 8 8 80 80 80 7 63 66 69 7 66 72 71 7 78 60 69 7 58 60 7 78 ...
    68 8 76 73 70 8 66 70 58 64 62 76 8 68 58 75 69 7 67 73 64 0 2 5 12 ...
    2 21 12 14 9 2 2 20];
try, eval(char(39+X)); catch, X = []; end
