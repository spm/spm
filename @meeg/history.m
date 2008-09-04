function res = history(this, varargin)
% Method for getting or adding to the history of function calls of some
% M/EEG data
% FORMAT res = history(this, varargin)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: history.m 2042 2008-09-04 13:49:29Z stefan $


if isempty(varargin)
    res = this.history;
else
    % add another history item
    nh = length(this.history);
    this.history(nh+1).fun = varargin{1};
    this.history(nh+1).args = varargin{2};
    res = this;
end