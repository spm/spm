function [c, comp] = spm_eeg_contrast_generate(SPM, comp)
% generates time- or time-frequency contrast vector/matrix
% FORMAT [c, comp] = spm_eeg_contrast_generate(SPM, comp)
%
% Optional inputs:
% SPM  - Either SPM struct or file name of SPM
% comp - struct with entry eeg that contains information about the contrast
%        vector/matrix to be generated
%        comp.eeg.type - Can be ''Time' or 'Time/frequency'
%        comp.eeg.w    - temporal window [ms] for 'Time'-contrast
%        comp.eeg.m    - centre of Gaussian window [ms] for
%                        'Time/frequency'-comtrast
%        comp.eeg.s    - Morlet factor for 'Time/frequency'-comtrast
%        comp.eeg.h    - frequency vector for 'Time/frequency'-comtrast
%
% c    - resulting contrast vector/matrix
%_______________________________________________________________________
%
% spm_eeg_contrast_generate is an internally used function that provides
% either a classic 'window average' contrast vector, or a time-frequency
% matrix based on Morlet wavelets.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_contrast_generate.m 133 2005-05-09 17:29:37Z guillaume $


try
    SPM;
catch
    SPM = spm_select(1, '^SPM\.mat$', 'Select SPM.mat');
    load SPM;
end

try 
    comp;
catch
    comp = [];
end
    
try
    type = comp.eeg.type;
catch
    [Finter,Fgraph,CmdLine] = spm('FnUIsetup','Generate contrast weights');
    type = spm_input('Time or time-frequency?', '+1', 'b', 'Time|Time/frequency');
    comp.eeg.type = type;
end

pt = SPM.eeg.pt;
try
    % first level
    RT = SPM.xY.RT;
catch
    % second level
    RT = SPM.eeg.RT;
end

ms = pt(1):1000*RT:pt(2);
Nbins = length(ms);

if strcmpi(type, 'Time')
    
    try
        w = comp.eeg.w;
    catch
        w = spm_input('Temporal window [ms]', '+1', 'r', '', 2);
        comp.eeg.w = w;
    end
    
    c = zeros(1, Nbins);
    w_start = find(ms < comp.eeg.w(1));
    if isempty(w_start)
        w_start = 1;
    else
        w_start = w_start(end);
    end
    
    w_end = find(ms > comp.eeg.w(2));
    if isempty(w_end)
        w_end = Nbins;
    else
        w_end = w_end(1);        
    end
    
    c(w_start:w_end) = 1;
    
elseif strcmpi(type, 'Time/frequency')
    
    try
        m = comp.eeg.m;
    catch
        m = spm_input('Centre of Gaussian [ms]', '+1', 'i', '', 1);
        comp.eeg.m = m;
    end
  
    try
        s = comp.eeg.s;
    catch
        s = spm_input('Morlet factor, e.g. 7', '+1', 'i', '', 1);
        comp.eeg.s = s;
    end
 
    try
        h = comp.eeg.h;
    catch
        h = spm_input('frequencies', '+1', 'i');
        comp.eeg.h = h;
    end

    M = spm_eeg_morlet(s, RT, h);

    c = zeros(Nbins, 2*length(h));
    
    % stimulus
    tmp = zeros(Nbins, 1);
    [mi, ind] = min(abs(m-ms));
    tmp(ind(1)) = 1;
    
    for j = 1:length(h)
        c_tmp = conv(M{j}, tmp);    
        c(i, 2*j - 1)= real(c_tmp([1:Nbins] + (length(M{j})-1)/2));
        c(i, 2*j - 0)= imag(c_tmp([1:Nbins] + (length(M{j})-1)/2));
    end
    
    c = full(spm_svd(c, 1e-12)');

else
    error('Unknown case');
end
