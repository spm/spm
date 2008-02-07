function [x,t,xi] = spm_boxcars(T,fs,len)
% Generate boxcar variable
% FORMAT [x,t,xi] = spm_boxcars(T,fs,len)
%
% T         Length of time series (secs)
% fs        Sampling rate, (Hz)
% len       Length of top of boxcar (secs)
%
% x         Event stream (1-event, 0-no event) (samples)
% t         time index (secs) eg. for plot(t,x) 
% xi        Sample numbers of events  (samples)
%
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_boxcars.m 1143 2008-02-07 19:33:33Z spm $

N=T*fs;
t=[1/fs:1/fs:T];

L=len*fs;
box=[zeros(L,1);ones(L,1)];
x=zeros(N,1);
xx=repmat(box,floor(N/length(box)),1);
x(1:length(xx))=xx;

xi=find(x==1);



