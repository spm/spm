function varargout = spm_matx(matfile,varargin)
% Utility to extract variables from a mat-file
% FORMAT [v1,v2,...] = spm_matx(matfile,n1,n2,n3,...)
%
% matfile      - name of Matlab mat-file
% n1,n2,n3,... - names of variables to be extracted from mat-file (strings)
% v2,v3,v3,... - variables extracted from mat-file
%_______________________________________________________________________
%
% This utility enables limited extraction (and renaming) of variables
% from a Matlab mat-file. It works by loading the mat-file in the
% functions workspace, and sending the named variables back as output
% arguments.
%
% E.g. : Matlab mat-file 'test.mat' in the current workign directory
% contains variables X, Y and Z, amongst other things. You already have
% an X, Y and Z variable, and don't wish to load the other variables
% from the mat-file. Using spm_matx, you can obtain X Y and Z as 'X1',
% 'Y1', and 'Z1' by:
%           [X1,Y1,Z1] = spm_matx('test','X','Y','Z')
%
%
% NB: Unfortunately this function has to load the *entire* mat-file
% into memory. To load individual data structures from mat-files
% requires an external C routine using the mat* functions of the mat.h
% Matlab C API. (See the MatLab Application Program Interface Guide",
% section 5 for details.) Since this is rather involved, and this
% functionality is intended for small mat-files, this m-file function
% implementation is more practicable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_matx.m 1143 2008-02-07 19:33:33Z spm $



%-Check arguments
%-----------------------------------------------------------------------
if nargin<2, varargout={}; return, end
matfile = [spm_str_manip(matfile,'s'),'.mat'];
if exist(matfile,'file')~=2, error(['invalid mat-file: ',matfile]), end


%-Note mat-file variables into structure M (feature introduced in v5.2)
%-----------------------------------------------------------------------
M = load(matfile);


%-Loop through varargin putting named variables into output arguments
%-----------------------------------------------------------------------
varargout = cell(1,nargout);
for i=1:min(max(1,nargout),nargin-1)
    if isfield(M,varargin{i})
        varargout{i} = getfield(M,varargin{i});
    else
        warning(['variable "',varargin{i},'" not found in mat-file: ',...
                matfile])
    end
end
