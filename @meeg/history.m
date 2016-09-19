function res = history(this, varargin)
% Method for getting or adding to the history of function calls of some
% M/EEG data
% FORMAT res = history(this, varargin)
% _______________________________________________________________________
% Copyright (C) 2008-2016 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: history.m 6883 2016-09-19 13:42:05Z vladimir $


if isempty(varargin)
    res = this.history;
else
    % add another history item
    if length(varargin) > 2 % To enable reset
        nh = 0;
        this.history = [];
    else
        nh = length(this.history);
    end

    if ischar(varargin{1})
        this.history(nh+1).fun = varargin{1};
        
        if isstruct(varargin{2}) && isfield(varargin{2}, 'D')
            if  isa(varargin{2}.D, 'meeg')
                varargin{2}.D = fullfile(varargin{2}.D);
            elseif isa(varargin{2}.D, 'cell')
                for i = 1:numel(varargin{2}.D)
                    if isa(varargin{2}.D{i}, 'meeg')
                        varargin{2}.D{i} = fullfile(varargin{2}.D{i});
                    end
                end
            end
        end

        this.history(nh+1).args = varargin{2};
                
        [dum, this.history(nh+1).ver] = spm('ver'); 
    elseif isstruct(varargin{1})
        this.history = varargin{1};
    end

    res = this;
end