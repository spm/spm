function res = units(obj, ind)
% Method for getting all units, over channels
% FORMAT res = units(obj, ind)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel

if obj.Nsamples>0
        res = {obj.channels.units};
else
    res = [];
end

if nargin > 1
    res = res(ind);
end

res = strvcat(res);
