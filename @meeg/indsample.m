function res = indsample(obj, t)
% Method for getting the sample closest to some time point
% FORMAT res = indsample(obj, t)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: indsample.m 1373 2008-04-11 14:24:03Z spm $

if obj.Nsamples>0
    T = time(obj);
    [m,res] = min(abs(T-t));
    
end
