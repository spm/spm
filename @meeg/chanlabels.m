function res = chanlabels(this, varargin)
% Method for getting/setting the channel labels
% FORMAT res = chanlabels(this, ind, label)
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: chanlabels.m 5933 2014-03-28 13:22:28Z vladimir $

if this.montage.Mind == 0
    if nargin == 3
        ind = varargin{1};
        label = varargin{2};

        if iscell(label) && length(label)>1
            if isnumeric(ind) && length(ind)~=length(label)
                error('Indices and values do not match');
            end

            if length(label)>1
                for i = 1:length(label)
                    for j = (i+1):length(label)
                        if strcmp(label{i}, label{j})
                            error('All labels must be different');
                        end
                    end
                end
            end

        end
    end

    res = getset(this, 'channels', 'label', varargin{:});
else
% case with an online montage applied
    if nargin == 3
        ind = varargin{1};
        label = varargin{2};

        if iscell(label) && length(label)>1
            if isnumeric(ind) && length(ind)~=length(label)
                error('Indices and values do not match');
            end

            if length(label)>1
                for i = 1:length(label)
                    for j = (i+1):length(label)
                        if strcmp(label{i}, label{j})
                            error('All labels must be different');
                        end
                    end
                end
            end

        end
        
        this.montage.M(this.montage.Mind) = getset(this.montage.M(this.montage.Mind), 'channels', 'label', varargin{:});
        res = this;
    else
        res = getset(this.montage.M(this.montage.Mind), 'channels', 'label', varargin{:});
    end
    
end