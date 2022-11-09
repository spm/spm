function s = spm_mg_switch(V)
% Switching output
% FORMAT s = spm_mg_switch(V)
%
% Switching output s: determined by voltage (V) dependent magnesium
% blockade parameters as per Durstewitz, Seamans & Sejnowski 2000.
%__________________________________________________________________________

% Rosalyn Moran
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


s = 1.50265./(1 + 0.33*exp(-0.06.*V));
