function planar = spm_eeg_planarchannelset(data)
% Define the planar gradiometer channel combinations
% FORMAT planar = spm_eeg_planarchannelset(data)
%
% The output cell array contains the horizontal label, vertical label and
% the label after combining the two.
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2013-2022 Wellcome Centre for Human Neuroimaging


if isa(data, 'cell') && ~isempty(data) && isa(data{1}, 'char')
    input.label = data;
else
    input = data;
end

senstype = lower(ft_senstype(input));

k = strfind(senstype, '_combined');
if ~isempty(k)
    senstype(k:(k+length('_combined')-1)) = [];
end        

try
    planar = ft_senslabel(senstype, 'output', 'planarcombined');
catch
    planar = ft_senslabel([senstype '_planar'], 'output', 'planarcombined');
end
