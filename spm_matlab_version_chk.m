function [status, fieldsUsed] = spm_matlab_version_chk(chk,tbx)
% Check a version number against a Toolbox version
% FORMAT [status, fieldsUsed] = spm_matlab_version_chk(chk,tbx)
% chk        - Version number to be checked {string}
% tbx        - Name of toolbox to check [Default: 'MATLAB']
%
% status     - Defines the outcome of the comparison
%              -1: Toolbox version is earlier than the user supplied version
%               0: Toolbox and user versions are the same
%               1: Toolbox version is later than the user supplied version
%                  Think of it this way, the sign of status is determined
%                  by MATLAB_TOOLBOX_VERSION - USER_VERSION (i.e., THE 
%                  VERSION YOU INPUT).
% fieldsUsed - deprecated [Returns {}]
%__________________________________________________________________________
%
% This function is deprecated, use SPM_CHECK_VERSION instead.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Darren Gitelman
% $Id: spm_matlab_version_chk.m 4145 2010-12-23 15:18:30Z guillaume $

%persistent runonce
%if isempty(runonce)
%    warning(['spm_matlab_version_check is deprecated. ',...
%        'Use spm_check_version instead.']);
%    runonce = 1;
%end

if nargin < 1, error('Please provide a version number to be checked.'); end
if nargin < 2, tbx = 'MATLAB'; end

status = spm_check_version(tbx,chk);

if nargout > 1, fieldsUsed = {}; end
