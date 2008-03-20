function res = badchannels(obj)
% Method for getting index vector of bad channels 
% FORMAT res = badchannels(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel

res = find(cat(1,obj.channels(:).bad));
