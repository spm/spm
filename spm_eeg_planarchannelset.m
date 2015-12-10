function planar = spm_eeg_planarchannelset(data)
% Define the planar gradiometer channel combinations
% FORMAT planar = spm_eeg_planarchannelset(data)
%
% The output cell array contains the horizontal label, vertical label and
% the label after combining the two.
%__________________________________________________________________________
% Copyright (C) 2013-2015 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_planarchannelset.m 6638 2015-12-10 15:43:57Z guillaume $


if isa(data, 'cell') && ~isempty(data) && isa(data{1}, 'char')
    data = struct('label',data);
end

planar = ft_senslabel(lower(ft_senstype(data)), 'output', 'planarcombined');
