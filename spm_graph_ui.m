function spm_graph_ui(varargin)
% User interface for editing graph attributes
% FORMAT spm_graph_ui(h)
% h - axis handle
%_______________________________________________________________________
% %W% Karl Friston %E%

%-Print warning of obsolescence
%-----------------------------------------------------------------------
warning('spm_graph_ui is grandfathered, use spm_results_ui(''PlotUi'',... instead')

%-Pass on arguments to spm_DesMtx
%-----------------------------------------------------------------------
spm_results_ui('PlotUi',varargin{:});

