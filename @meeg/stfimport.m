function this = stfimport(this, stfname)
% Use an SPM5 stf file to assign X_plot2D and Y_plot2D
% FORMAT ctfimport(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: stfimport.m 1255 2008-03-27 19:44:05Z vladimir $

if nargin < 2
    P = spm_select(1, '\.mat$', 'Select sensor template file', [], fullfile(spm('dir'), 'EEGtemplates'));
else
    [tmp1,tmp2,ext] = fileparts(stfname);
    if ~strcmp(ext, '.mat')
        stfname = [stfname '.mat'];
    end
    
    P = fullfile(spm('dir'), 'EEGtemplates', stfname);
end

stf = load(P); % must contain Cpos, Cnames

% identify channels in template file and copy their coordinates
chan = chanlabels(this);
for i = 1:nchannels(this)
    index = [];
    for j = 1:stf.Nchannels
        if ~isempty(find(strcmpi(chan{i}, stf.Cnames{j})))
            index = [index j];
            this = chantype(this, i, 'EEG');
        end
    end

    if isempty(index)
        warning(sprintf('No channel named %s found in channel template file.', chan{i}));
    else
        X = stf.Cpos(1, index(1));
        [this.channels(i).X_plot2D] = X;

        Y = stf.Cpos(2, index(1));
        [this.channels(i).Y_plot2D] = Y;
    end
end

[res this] = checkmeeg(struct(this));
this = meeg(this);
