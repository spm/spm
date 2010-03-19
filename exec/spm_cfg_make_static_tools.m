function spm_cfg_make_static_tools
% Collect all toolboxes in the current SPM installation and create code to
% initialise their batch configuration in a compiled SPM version. 
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Volkmar Glauche
% $Id: spm_cfg_make_static_tools.m 3788 2010-03-19 15:58:33Z volkmar $ 

fid = fopen(fullfile(spm('dir'),'config','spm_cfg_static_tools.m'),'w');
fprintf(fid,'function values = spm_cfg_static_tools\n');
fprintf(fid,...
    '%% Static listing of all batch configuration files in the SPM toolbox folder\n');
% create code to insert toolbox config
%-Toolbox autodetection
%-Get the list of toolbox directories
tbxdir = fullfile(spm('Dir'),'toolbox');
d  = dir(tbxdir); d = {d([d.isdir]).name};
dd = regexp(d,'^\.');
%(Beware, regexp returns an array if input cell array is of dim 0 or 1)
if ~iscell(dd), dd = {dd}; end
d  = {'' d{cellfun('isempty',dd)}};
ft = {};
ftc = {};
%-Look for '*_cfg_*.m' or '*_config_*.m' files in these directories
for i=1:length(d)
    d2 = fullfile(tbxdir,d{i});
    di = dir(d2); di = {di(~[di.isdir]).name};
    f2 = regexp(di,'.*_cfg_.*\.m$');
    if ~iscell(f2), f2 = {f2}; end
    fi = {di{~cellfun('isempty',f2)}};
    if ~isempty(fi)
        ft = [ft(:); fi(:)];
    else
        % try *_config_*.m files, if toolbox does not have '*_cfg_*.m' files
        f2 = regexp(di,'.*_config_.*\.m$');
        if ~iscell(f2), f2 = {f2}; end
        fi = {di{~cellfun('isempty',f2)}};
        ftc = [ftc(:); fi(:)];
    end;
end
if ~isempty(ft)||~isempty(ftc)
    if isempty(ft)
        ftstr = '';
    else
        ft = cellfun(@(cft)strtok(cft,'.'),ft,'UniformOutput',false);
        ftstr  = sprintf('%s ', ft{:});
    end
    if isempty(ftc)
        ftcstr = '';
    else
        ftc = cellfun(@(cftc)strtok(cftc,'.'),ftc,'UniformOutput',false);
        ftcstr = sprintf('cfg_struct2cfg(%s) ', ftc{:});
    end
    fprintf(fid,'values = {%s %s};\n', ftstr, ftcstr);
end
fclose(fid);