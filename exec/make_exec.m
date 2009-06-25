function make_exec
% SPM can be compiled using Matlab 7.
% This will generate a standalone program, which can be run
% outside Matlab, and therefore does not use up a Matlab license.
% The executable needs to be dynamically linked with runtime libraries
% bundled with the Matlab package.  Before executing, ensure that your
% LD_LIBRARY_PATH is set appropriately, so that $MATLAB/bin/$arch
% is included. LD_LIBRARY_PATH or PATH should also contain the directory
% containing the *.ctf files
%
% Software compiled with later MATLAB versions may need 
% LD_LIBRARY_PATH=$MATLAB/bin/$arch:$MATLAB/sys/os/$arch
%
% Note that when compiling, it is important to be careful with your
% startup.m file.  See the following link for more information:
% http://www.mathworks.com/support/solutions/data/1-QXFMQ.html?1-QXFMQ
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: make_exec.m 3226 2009-06-25 18:13:11Z volkmar $ 

%=======================================================================
%-Files to include explicitly
%=======================================================================
% Matlab compiler will include all files that are referenced explicitly
% in the compiled code. If there is a file missing (e.g. data file,
% (f)eval'ed file), add it to the 'includefiles' list with its full
% path. If you want to add a directory and all of its contents, add it to
% the 'includedirs' list with its full path.
% By default, all files in spm('dir') are included and all subdirectories
% except matlabbatch, config and exec. Selected files from these
% directories will be included if they are referenced from any compiled
% function. Including spm('dir') instead would cause unwanted side
% effects during batch initialisation.

%-List of files
%-----------------------------------------------------------------------
[includefiles includedirs] = cfg_getfile('FPList',spm('dir'),'.*');

%-Clean up list of directories
%-----------------------------------------------------------------------
includedirs = includedirs(~(strcmp(includedirs,fullfile(spm('dir'),'config'))| ...
                            strcmp(includedirs,fullfile(spm('dir'),'exec'))| ...
                            strcmp(includedirs,fullfile(spm('dir'),'matlabbatch'))));

includefiles = [includefiles; includedirs];
%-Add '-a' switch for each file to include
%-----------------------------------------------------------------------
sw           = cell(size(includefiles));
[sw{:}]      = deal('-a');
tmp          = [sw includefiles]';
includefiles = tmp(:);

%=======================================================================
%-Configuration management
%=======================================================================
% 1) locate all currently used batch configs
% 2) copy them to fullfile(spm('dir'),'exec') with unique names
% 3) create fullfile(spm('dir'),'exec','cfg_master.m') which calls the
%    configs at runtime 
% compiled spm_jobman will execute this file and add the applications

%-Locate batch configs and copy them
%-----------------------------------------------------------------------
apps = which('cfg_mlbatch_appcfg','-all');
cfgfiles = cell(1,2*numel(apps)+2);
for k = 1:numel(apps)
    cfgfiles{2*k-1} = '-a';
    cfgfiles{2*k}   = fullfile(spm('dir'),'exec',sprintf('cfg_mlbatch_appcfg_%d.m',k));
    copyfile(apps{k}, cfgfiles{2*k});
end

%-Create code for cfg_master.m
%-----------------------------------------------------------------------
cfgfiles{end-1} ='-a';
cfgfiles{end}   = fullfile(spm('dir'),'exec','cfg_master.m');
fid             = fopen(cfgfiles{end},'w');
fprintf(fid,'function cfg_master\n');
for k = 1:numel(apps)
    fprintf(fid,'[cfg, def] = cfg_mlbatch_appcfg_%d;\n', k);
    if strcmp(apps{k},...
              fullfile(spm('dir'),'config','cfg_mlbatch_appcfg.m'))
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
                ft = {ft{:} fi{:}};
            else
                % try *_config_*.m files, if toolbox does not have '*_cfg_*.m' files
                f2 = regexp(di,'.*_config_.*\.m$');
                if ~iscell(f2), f2 = {f2}; end
                fi = {di{~cellfun('isempty',f2)}};
                ftc = {ftc{:} fi{:}};
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
            % assume that tools are the last thing in SPM config
            fprintf(fid,'cfg.values{end}.values = {%s %s};\n', ftstr, ftcstr);
        end
    end
    fprintf(fid,'cfg_util(''addapp'', cfg, def);\n');
end
fclose(fid);

%=======================================================================
%-Compile
%=======================================================================
mcc('-m','-v','-o',['spm_'    computer],'exec_spm.m'   ,'spm_load.m', '-I',spm('Dir'),includefiles{:},cfgfiles{:})
mcc('-m','-v','-o',['jobman_' computer],'exec_jobman.m','spm_load.m', '-I',spm('Dir'),includefiles{:},cfgfiles{:})
