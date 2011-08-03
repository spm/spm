function spm_chi2_plot(action,varargin)
% Display a plot showing convergence of an optimisation routine.
% FORMAT spm_chi2_plot('Init',title,ylabel,xlabel)
% Initialise the plot in the 'Interactive' window.
%
% FORMAT spm_chi2_plot('Set',value)
% Update the plot.
%
% FORMAT spm_chi2_plot('Clear')
% Clear the 'Interactive' window.
%__________________________________________________________________________
%
% This function is deprecated, use SPM_PLOT_CONVERGENCE instead.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_chi2_plot.m 4418 2011-08-03 12:00:13Z guillaume $

persistent runonce
if isempty(runonce)
    warning(['spm_chi2_plot is deprecated. ',...
        'Use spm_plot_convergence instead.']);
    runonce = 1;
end

if ~nargin, action = 'Init'; end
spm_plot_convergence(action,varargin{:});
