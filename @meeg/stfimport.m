function obj = stfimport(obj, stfname)
% Use an SPM5 stf file to assign X_plot2D and Y_plot2D
% FORMAT ctfimport(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id $

if nargin < 2
    P = spm_select(1, '\.mat$', 'Select sensor template file', [], fullfile(spm('dir'), 'EEGtemplates'));
else
    if ~strcmp(stfname((end-3):end), '.mat')
        stfname = [stfname '.mat'];
    end
    
    P = fullfile(spm('dir'), 'EEGtemplates', stfname);
end

stf=load(P); % must contain Cpos, Cnames

[sel1, sel2] = match_str(chanlabels(obj), stf.Cnames);

X = num2cell(stf.Cpos(1, sel2));
[obj.channels(sel1).X_plot2D] = deal(X{:});

Y = num2cell(stf.Cpos(2, sel2));
[obj.channels(sel1).Y_plot2D] = deal(Y{:});

[res obj] = checkmeeg(struct(obj));
obj = meeg(obj);
