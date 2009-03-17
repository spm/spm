function res = indsample(this, t)
% Method for getting the sample closest to some time point.
% The time is specified in sec.
% FORMAT res = indsample(this, t)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: indsample.m 2889 2009-03-17 12:02:04Z vladimir $

if this.Nsamples>0
    T = time(this);
    for i = 1:length(t)
        [m,res(i)] = min(abs(T-t(i)));
        if m > (1/this.Fsample)
            warning('Could not find an index matching the requested time %d sec', t(i));
            res(i) = nan;
        end
    end
end
