% SPM5 (c) 1991,1994-2003,2005
% Statistical Parametric Mapping           -                        SPM5
%_______________________________________________________________________
%  ___  ____  __  __
% / __)(  _ \(  \/  )  
% \__ \ )___/ )    (   Statistical Parametric Mapping
% (___/(__)  (_/\/\_)  SPM - http://www.fil.ion.ucl.ac.uk/spm/
%_______________________________________________________________________
%
% This Contents.m file holds the version ID for this release of Matlab,
% and contains a manifest of the included functions and their version numbers.
%
% SPM5 is written for Matlab 6.5.1 or Matlab 7.0.1 under UNIX and Windows
% ( Compiled binaries of external MEX functions are provided for:       )
% (                   Solaris2, Linux and Windows                       )
%
% See spm.man for details of this release.
% See README.txt for information on installation and getting started.
% See spm_motd.man for last minute release details.
%
% SPM (being the collection of files given in the manifest below) is
% free but copyright software, distributed under the terms of the GNU
% General Public Licence as published by the Free Software Foundation
% (either version 2, as given in file spm_LICENCE.man, or at your option,
% any later version). Further details on "copyleft" can be found at
% http://www.gnu.org/copyleft/.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% $Id: Contents.m 123 2005-05-06 12:15:13Z john $

%=======================================================================
% PROGRAMMERS NOTE:
% This <Contents.m> is the contents file for SPM, used by spm('Ver') to
% recover the version number and copyright information.
%   Line1: Version (first word) & copyright information (rest of line).
%   Line2: One line description
% Matlab's ver also uses <Contents.m> files to identify toolbox versions.
%   Line1: Toolbox Description
%   Line2: Version xxx dd-mmm-yyyy
%=======================================================================
%%% Extract LastRevision number from 'Id' tag in SPM files
% mext = {'.m','.c','.h','.man','.xml', ...       %- source code
%   	  '.dll','.mexmac','.mexsol','.mexglx'};  %- MEX
% spmdir = spm('Dir');
% d = dir(fullfile(spmdir,'*'));
% f = {d(~[d.isdir]).name};
% d = {d([d.isdir]).name}; d = {d{~ismember(d,{'.' '..'})}};
% L = max(cellfun('length',f));
% pat = '\$Id: (\S+) (\d+) ([0-9-]+) ([0-9:]+Z) (\w+) \$';
% for i=1:length(f)
%     [pathstr, name, ext] = fileparts(f{i});
%     if ismember(ext,mext)
%   	  fid = fopen(fullfile(spmdir,f{i}),'r');
%   	  if fid == -1, continue; end
%   	  V = 'none';
%   	  while 1
%   		  tline = fgetl(fid);
%   		  if ~ischar(tline), break, end
%   		  tok = regexp(tline, pat, 'tokens');
%   		  if ~isempty(tok), V = tok{1}{2}; break; end
%   	  end
%   	  fclose(fid);
%   	  fprintf('%% %s %s%s\n',f{i},blanks(L-length(f{i})),V);
%     end
% end
%=======================================================================
