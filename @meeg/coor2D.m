function [res, plotind] = coor2D(this, ind, val, mindist)
% x and y coordinates of channels in 2D plane
% FORMAT coor2D(this)
%__________________________________________________________________________

% Vladimir Litvak, Laurence Hunt
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


megind = indchantype(this, {'MEG', 'MEGPLANAR', 'MEGCOMB'});
eegind = indchantype(this, {'EEG'});
otherind = setdiff(1:nchannels(this), [megind eegind]);

if nargin==1 || isempty(ind)
    if nargin<3 || (size(val, 2)<nchannels(this))
        if ~isempty(megind)
            ind = megind;
        elseif ~isempty(eegind)
            ind = eegind;
        else
            ind = 1:nchannels(this);
        end
    else
        ind = 1:nchannels(this);
    end
elseif ischar(ind)
    switch upper(ind)
        case 'MEG'
            ind = megind;
        case 'EEG'
            ind = eegind;
        case ':'
            ind = 1:nchannels(this);
        otherwise
            ind = otherind;
    end
end

if nargin < 3 || isempty(val)
    if ~isempty(intersect(ind, megind))
        
        if this.montage.Mind==0
            if ~any(cellfun('isempty', {this.channels(megind).X_plot2D this.channels(megind).Y_plot2D}))
                meg_xy = [this.channels(megind).X_plot2D; this.channels(megind).Y_plot2D];
            elseif all(cellfun('isempty', {this.channels(megind).X_plot2D this.channels(megind).Y_plot2D}))
                meg_xy = grid(length(megind));
            else
                error('Either all or none of MEG channels should have 2D coordinates defined.');
            end
        else
            if ~any(cellfun('isempty', {this.montage.M(this.montage.Mind).channels(megind).X_plot2D ...
                    this.montage.M(this.montage.Mind).channels(megind).Y_plot2D}))
                meg_xy = [this.montage.M(this.montage.Mind).channels(megind).X_plot2D; ...
                    this.montage.M(this.montage.Mind).channels(megind).Y_plot2D];
            elseif all(cellfun('isempty', {this.montage.M(this.montage.Mind).channels(megind).X_plot2D...
                    this.montage.M(this.montage.Mind).channels(megind).Y_plot2D}))
                meg_xy = grid(length(megind));
            else
                error('Either all or none of EEG channels should have 2D coordinates defined.');
            end
        end
    end
    
    if ~isempty(intersect(ind, eegind))
        if this.montage.Mind==0
            if ~any(cellfun('isempty', {this.channels(eegind).X_plot2D this.channels(eegind).Y_plot2D}))
                eeg_xy = [this.channels(eegind).X_plot2D; this.channels(eegind).Y_plot2D];
            elseif all(cellfun('isempty', {this.channels(eegind).X_plot2D this.channels(eegind).Y_plot2D}))
                eeg_xy = grid(length(eegind));
            else
                error('Either all or none of EEG channels should have 2D coordinates defined.');
            end
        else
            if ~any(cellfun('isempty', {this.montage.M(this.montage.Mind).channels(eegind).X_plot2D ...
                    this.montage.M(this.montage.Mind).channels(eegind).Y_plot2D}))
                eeg_xy = [this.montage.M(this.montage.Mind).channels(eegind).X_plot2D; ...
                    this.montage.M(this.montage.Mind).channels(eegind).Y_plot2D];
            elseif all(cellfun('isempty', {this.montage.M(this.montage.Mind).channels(eegind).X_plot2D...
                    this.montage.M(this.montage.Mind).channels(eegind).Y_plot2D}))
                eeg_xy = grid(length(eegind));
            else
                error('Either all or none of EEG channels should have 2D coordinates defined.');
            end
        end
    end
    
    if ~isempty(intersect(ind, otherind))
        other_xy = grid(length(otherind));
    end
    
    xy = zeros(2, length(ind));
    plotind = zeros(1, length(ind));
    for i = 1:length(ind)
        [found, loc] = ismember(ind(i), megind);
        if found
            xy(:, i) = meg_xy(:, loc);
            plotind(i) = 1;
        else
            [found, loc] = ismember(ind(i), eegind);
            if found
                xy(:, i) = eeg_xy(:, loc);
                plotind(i) = 2;
            else
                [found, loc] = ismember(ind(i), otherind);
                if found
                    xy(:, i) = other_xy(:, loc);
                    plotind(i) = 3;
                end
            end
        end
    end
    
    if nargin > 3 && ~isempty(mindist)
        xy = shiftxy(xy,mindist);
    end
    
    res = xy;
else
    if this.montage.Mind==0
        this = getset(this, 'channels', 'X_plot2D', ind, val(1, :));
        this = getset(this, 'channels', 'Y_plot2D', ind, val(2, :));
    else
        this.montage.M(this.montage.Mind) = getset(this.montage.M(this.montage.Mind), 'channels', 'X_plot2D', ind, val(1, :));
        this.montage.M(this.montage.Mind) = getset(this.montage.M(this.montage.Mind), 'channels', 'Y_plot2D', ind, val(2, :));
    end
    res = this;
end


function xy = grid(n)

ncol = ceil(sqrt(n));
x = 0:(1/(ncol+1)):1;
x = 0.9*x+0.05;
x = x(2:(end-1));
y = fliplr(x);
[X, Y] = meshgrid(x, y);
xy = [X(1:n); Y(1:n)];


function xy = shiftxy(xy,mindist)

x = xy(1,:);
y = xy(2,:);

l=1;
i=1; %filler
mindist = mindist/0.999; % limits the number of loops
while (~isempty(i) && l<50)
    xdiff = repmat(x,length(x),1) - repmat(x',1,length(x));
    ydiff = repmat(y,length(y),1) - repmat(y',1,length(y));
    xydist= sqrt(xdiff.^2 + ydiff.^2); %euclidean distance between all sensor pairs
    
    [i,j] = find(xydist<mindist*0.999);
    rm=(i<=j); i(rm)=[]; j(rm)=[]; %only look at i>j
    
    for m = 1:length(i);
        if (xydist(i(m),j(m)) == 0)
            diffvec = [mindist./sqrt(2) mindist./sqrt(2)];
        else
            xydiff = [xdiff(i(m),j(m)) ydiff(i(m),j(m))];
            diffvec = xydiff.*mindist./xydist(i(m),j(m)) - xydiff;
        end
        x(i(m)) = x(i(m)) - diffvec(1)/2;
        y(i(m)) = y(i(m)) - diffvec(2)/2;
        x(j(m)) = x(j(m)) + diffvec(1)/2;
        y(j(m)) = y(j(m)) + diffvec(2)/2;
    end
    l = l+1;
end

xy = [x; y];
