function spm_spy(X,Markersize,m)
% pretty version of spy
% FORMAT spm_spy(X,Markersize,m)
% X    - sparse {m x n} matrix
%
% See also: spy
%__________________________________________________________________________
% Copyright (C) 1994-2011 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_spy.m 5784 2013-12-05 17:41:58Z karl $


% defaults
%--------------------------------------------------------------------------
if nargin < 2, Markersize = 16; end
if nargin < 3, m = max(max(X)); end


spy(X > m/256,Markersize,'.c'), hold on
spy(X > m/2,Markersize,'.b'), hold on
spy(X > (m/2 + m/4),Markersize,'.k'), hold off
axis normal
