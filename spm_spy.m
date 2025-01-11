function spm_spy(X,Markersize,m)
% Pretty version of spy
% FORMAT spm_spy(X,Markersize,m)
% X    - sparse {m x n} matrix
%
% See also: spy
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 1994-2022 Wellcome Centre for Human Neuroimaging


% defaults
%--------------------------------------------------------------------------
if nargin < 1, X = defaultspy;  end
if nargin < 2, Markersize = 16; end
if nargin < 3, m = max(max(X)); end

hold off
for p = flip(1:4)
    [i,j] = find(X > m/p);
    plot(j,i,'.','Markersize',Markersize,'Color',[1,1,1]*(p - 1)/4)
    hold on
end
hold off
axis ij
[x,y] = size(X);
set(gca,'XLim',[1,y],'YLim',[1,x])


% test routine
%--------------------------------------------------------------------------
function X = defaultspy
X = fullfile(spm('Dir'),'help','images','karl.jpg');
X = sparse((sum(imread(X),3) < 350));


