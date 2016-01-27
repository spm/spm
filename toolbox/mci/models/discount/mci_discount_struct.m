function [M,U] = mci_discount_struct (Nobs)
% Set up data structures for discounting model
% FORMAT [M,U] = mci_discount_struct (Nobs)
%
% Nobs      number of data points
%
% M         model structure
% U         U.X is the design matrix
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: mci_discount_struct.m 6697 2016-01-27 14:57:28Z spm $

% Number of data points
try, T=Nobs; catch, T=100; end

% Set rewards
rmax=70; % maximum reward (pounds)
U.r1=floor(rand(T,1)*24+16); % low
U.r2=floor(U.r1+rand(T,1).*(rmax-U.r1)); % high

% Set delays
tmax=75; % maximum delay (weeks)
U.t1=floor(rand(T,1)*54+1); % short
U.t2=floor(U.t1+rand(T,1).*(tmax-U.t1)); % long

M.L='mci_discount_like';
M.dL='mci_discount_deriv';
M.IS='mci_discount_gen';

Np=2;
M.pE=zeros(Np,1);
M.pC=diag([1 1]);
M.l=1;
M.T=T;
M.N=T;


