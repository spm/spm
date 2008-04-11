function spm_tbx_config2cfg(c0)
% Convert SPM5 toolbox configuration to Matlabbatch
% FORMAT spm_tbx_config2cfg(c)
% Input:
% c      - SPM5 toolbox configuration structure
% Output: (written to disk in the current working directory)
% tbx_cfg_<toolboxtag>.m - Code to generate a Matlabbatch configuration
%                          tree similar to the SPM5 configuration struct
% tbx_def_<toolboxtag>.m - Code to set toolbox defaults.
%
% Both files should be placed in the root directory of the toolbox
% instead of the old config file. They will be picked up during spm
% initialisation by spm_cfg.m/spm_def.m.
% The full subscript reference path will become 
% spmjobs{}.tools{}.<toolboxtag> in the configuration tree and
% spmjobs.tools.<toolboxtag> in the internal defaults structure used to
% initialise the configuration.
%
% CAVE: No code is generated for subfunctions that are present in the old
% config file. This code has to be transferred manually. A transition
% from .vfiles to .vout callbacks is strongly encouraged. This requires
% - computation functions to return a single results variable (any kind
%   of MATLAB variable allowed, but struct or cell preferred to combine
%   multiple results)
% - a .vout callback to describe subscript references into this output
%   variable for each individual output.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Volkmar Glauche
% $Id: spm_tbx_config2cfg.m 1374 2008-04-11 14:36:35Z volkmar $

% Convert to cfg_ tree. This will produce warnings if some elements could
% not be converted properly.
c = cfg_struct2cfg(c0);

% generate code for configuration
[cstr tag] = gencode(c,'',{},'',cfg_tropts(cfg_findspec, 0, Inf, 0, Inf, true));
cfgname = sprintf('tbx_cfg_%s', tag);
fid = fopen(sprintf('%s.m', cfgname),'w');
fprintf(fid, 'function %s = %s\n', tag, cfgname);
fprintf(fid, ...
        '%% MATLABBATCH Configuration file for toolbox ''%s''\n', c.name);
fprintf(fid, ...
        '%% This code has been automatically generated.\n', c.name);
fprintf(fid, '%s\n', cstr{:});
fclose(fid);

% harvest defaults
[u1 d] = harvest(c,c,true,false);
% generate code for defaults
dstr = gencode(d, tag);
defname = sprintf('tbx_def_%s', tag);
fid = fopen(sprintf('%s.m', defname),'w');
fprintf(fid, 'function varargout = %s\n', defname);
fprintf(fid, ...
        '%% MATLABBATCH Defaults file for toolbox ''%s''\n', c.name);
fprintf(fid, ...
        '%% This code has been automatically generated.\n', c.name);
fprintf(fid, '%s\n', dstr{:});
fprintf(fid, 'varargout{1} = %s;\n', tag);
fprintf(fid, 'varargout{2} = ''%s'';\n', tag);
fclose(fid);
