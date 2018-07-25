function spm_spy(X,Markersize,m)
% Pretty version of spy
% FORMAT spm_spy(X,Markersize,m)
% X    - sparse {m x n} matrix
%
% See also: spy
%__________________________________________________________________________
% Copyright (C) 1994-2018 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_spy.m 7380 2018-07-25 09:49:27Z guillaume $


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
X = fullfile(spm('Dir'),'help','images','karl.jpg');
X = sparse((sum(imread(X),3)<350));
