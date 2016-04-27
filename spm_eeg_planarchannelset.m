function planar = spm_eeg_planarchannelset(data)
% Define the planar gradiometer channel combinations
% FORMAT planar = spm_eeg_planarchannelset(data)
%
% The output cell array contains the horizontal label, vertical label and
% the label after combining the two.
%__________________________________________________________________________
% Copyright (C) 2013-2015 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_planarchannelset.m 6781 2016-04-27 10:01:38Z vladimir $


if isa(data, 'cell') && ~isempty(data) && isa(data{1}, 'char')
    input.label = data;
else
    input = data;
end

planar = ft_senslabel(lower(ft_senstype(input)), 'output', 'planarcombined');
