function res = indsample(this, t)
% Method for getting the sample closest to some time point
% FORMAT res = indsample(this, t)
% this       - MEEG object
% t          - vector of time points in seconds
%
% res        - vector of sample indices matching time points
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: indsample.m 2921 2009-03-23 17:59:50Z guillaume $

res = NaN(1,length(t));
if this.Nsamples > 0
    T = time(this);
    for i = 1:length(t)
        [m,res(i)] = min(abs(T-t(i)));
        if m > (1/this.Fsample)
            warning('Could not find an index matching the requested time %d sec', t(i));
            res(i) = NaN;
        end
    end
end
