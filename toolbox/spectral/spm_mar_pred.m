function [y,y_pred] = spm_mar_pred (X,mar)
% Get predictions from MAR model
% FORMAT [y,y_pred] = spm_mar_pred (X,mar)
%
% X              T-by-d matrix containing d-variate time series0)
%
% mar            see spm_mar.m for data structure
%
% y              Target values
% y_pred         Predicted values
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_mar_pred.m 1143 2008-02-07 19:33:33Z spm $

d=size(X,2);    % dimension of time series
N=size(X,1);    % length of time series
p=mar.p;

% Embedding of multiple time series 
% giving x=[(x1(t-p) x2(t-p) .. xd(t-p)) (x1(t-p+1) x2(t-p+1)..xd(t-p+1)) ...
%           (x1(t-1) x2(t-1) .. xd(t-1))] on each row

x=[];
for i=1:p,
  tmpx=X(i:N-p+i,:);
  x=[x,tmpx];
end

% Reorder columns of x 
% giving x=[(x1(t-1) x2(t-1) .. xd(t-1)) (x1(t-2) x2(t-2)..xd(t-2)) ...
%           (x1(t-p) x2(t-p) .. xd(t-p))] on each row
for i=1:p
  start=(i-1)*d+1;
  stop=start+d-1;
  chunk(i).x=x(:,[start:1:stop]);
end
x=[];
for i=p:-1:1,
  x=[x chunk(i).x];
end
% and remove last row
Nrows=size(x,1);
x=x(1:Nrows-1,:);

y=X([p+1:1:N],:);

y_pred = x*mar.wmean;