function [hReg,SPM,VOL,xX,xSDM] = spm_results
% Display and analysis of regional effects (Grandfathered)
% FORMAT spm_results
%_______________________________________________________________________
%
% spm_results functionality has been absorbed into spm_results_ui.
%
% Note that spm_results_ui is a function, but the results section still
% relies on having the SPM results summary structures available in the
% base workspace. Therefore, calls to spm_results_ui must collect the
% appropriate return arguments [hReg,SPM,VOL,xX,xSDM].
%
% For interim backwards compatability, spm_results is maintained as a
% gateway to spm_results_ui. If called without at least four return
% arguments, the required variables are assigned in the base
% workspace.
%
%___________________________________________________________________________
% %W% Andrew Holmes %E%


%-Print warning of obsolescence
%-----------------------------------------------------------------------
warning('spm_results is grandfathered, use spm_results_ui instead')


%-Pass on to spm_results_ui & return necessary arguments
%-----------------------------------------------------------------------
varargout = cell(1,5);
[varargout{1:5}] = spm_results_ui;


%-If not 5 output arguments, then use assignin('Base',...
%-----------------------------------------------------------------------
if nargout<5, assignin('base','xSDM',	varargout{5}), end
if nargout<4, assignin('base','xX',	varargout{4}), end
if nargout<3, assignin('base','VOL',	varargout{3}), end
if nargout<2, assignin('base','SPM',	varargout{2}), end
if nargout<1, assignin('base','hReg',	varargout{1}), end
