function res = indsample(obj, t)
% Method for getting the sample closest to some time point.
% The time is specified in sec.
% FORMAT res = indsample(obj, t)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: indsample.m 2418 2008-10-30 10:33:39Z stefan $

if obj.Nsamples>0
    T = time(obj);
    [m,res] = min(abs(T-t));
    
end
