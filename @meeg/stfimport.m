function obj = stfimport(obj, stfname)
% Use an SPM5 stf file to assign X_plot2D and Y_plot2D
% FORMAT ctfimport(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: stfimport.m 1236 2008-03-20 18:15:33Z stefan $

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
chan = chanlabels(obj);
for i = 1:nchannels(obj)
    index = [];
    for j = 1:stf.Nchannels
        if ~isempty(find(strcmpi(deblank(chan(i,:)), stf.Cnames{j})))
            index = [index j];
        end
    end

    if isempty(index)
        warning(sprintf('No channel named %s found in channel template file.', deblank(chan(i,:))));
    else
        X = stf.Cpos(1, index(1));
        [obj.channels(i).X_plot2D] = X;

        Y = stf.Cpos(2, index(1));
        [obj.channels(i).Y_plot2D] = Y;

    end
end



[res obj] = checkmeeg(struct(obj));
obj = meeg(obj);
