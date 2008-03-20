function obj = putdata(obj, Ichannels, Isamples, Itrials, data)
% Writes data to file_array 
% FORMAT res = selectdata(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: selectdata.m 1125 2008-01-30 12:12:18Z vladimir $

if obj.Nsamples == 0
    res = [];
    return;
end

if nargin<5
    error('Insufficient number of arguments for data selection');
end

obj.data.y(Ichannels, Isamples, Itrials) = data;