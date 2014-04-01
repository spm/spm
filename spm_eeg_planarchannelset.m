function [planar] = spm_eeg_planarchannelset(input)

% FUNCTION that defines the planar gradiometer channel combinations
% The output cell-array contains the horizontal label, vertical label
% and the label after combining the two.
%
% Use as
%   [planar] =  spm_eeg_planarchannelset(data)
%
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% $Id: spm_eeg_planarchannelset.m 5936 2014-04-01 09:40:26Z vladimir $

if isa(input, 'cell')    && ~isempty(input) && isa(input{1}, 'char')
    data.label = input;
else
    data = input;
end

planar = ft_senslabel(lower(ft_senstype(data)), 'output', 'planarcombined');
