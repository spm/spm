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
if nargin < 1, X = defaultspy; end
if nargin < 2, Markersize = 16; end
if nargin < 3, m = max(max(X)); end

spy(X > m/256,Markersize,'.c'), hold on
spy(X > m/2,Markersize,'.b'), hold on
spy(X > (m/2 + m/4),Markersize,'.k'), hold off
axis normal
if nargin < 1 && rand < 0.05, try, life(X), end, end


function X = defaultspy
X = fullfile(spm('Dir'),'help','images','karl.jpg');
X = sparse((sum(imread(X),3)<350));


function life(X)
h = get(gca,'Children'); delete(h(2:3)); h = h(1);
set(h,'Marker','o','MarkerSize',2); pause(2);
while true
    [m,n] = size(X);
    a = circshift(1:n,1); b = circshift(1:n,-1);
    c = circshift(1:m,1); d = circshift(1:m,-1);
    Y = X(:,a) + X(:,b) + X(c,:) + X(d,:) + ...
        X(c,a) + X(d,b) + X(c,b) + X(d,a);
    X = (X & (Y == 2)) | (Y == 3);
    [i,j] = find(X);
    set(h,'XData',j,'YData',i);
    drawnow;
end
