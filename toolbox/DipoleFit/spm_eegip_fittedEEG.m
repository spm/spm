function [Vfit,SV] = spm_eegip_fittedEEG(sdip,ind_s,varargin)
% FORMAT function [Vfit,SV] = spm_eegip_fittedEEG(sdip,ind_s,avg,el_set,type,maxphi)
% or
% FORMAT [Vfit,SV] = fittedEEG(sdip,ind_s)
%
% Calculate the fitted EEG potential from the source
% dipoles. And if required the interpolated scalp map.
%
% IN:
%   - sdip  :  fitted source dipoles structure.
%             or output from summarised results.
%   - ind_s : index of the source set(s) to be used
%             By default: 1:sdip.n_seeds
%   - avg   : data can be averaged over t_wind and only one image produced (avg=1)
%             or one image is created per time bin (avg=0, default), SV(1...Ntb).
%   - el_set: electrodes set, an index refering to the sets defined
%              in spm_EEGodel, or the electrodes coordinates (on a shere)
%   - types : 'simple' EEG interpolation (1, default) or Laplacian (2)
%   - maxphi: Maximum phi angle used for the interpolation and display (defaut 2*pi/3)
% OUT:
%   - Vfit  : cell array of fitted EEG values. One cell per 
%             source set, over the time window.
%   - SV    : cell array of interpolated scalp EEG map.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips,
% $Id: spm_eegip_fittedEEG.m 143 2005-05-11 17:13:13Z christophe $


calcSV = 0;
dispSV = 0;

if nargin<1
    Pdip = spm_get(1,'*dip*.mat','Select dipole file');
    load(Pdip);
    if ~exist('sdip') & exist('result')
        sdip = result;
    end
    calcSV = 1;
    dispSV = 1;
end

if nargin<2 
    if sdip.n_seeds>1
        % Which source set ?
        [ind_s] = spm_input(['Source set: [',num2str(sdip.n_seeds),'] ?'],'+1','e',1) ;
    else
        ind_s = 1
    end
elseif isempty(ind_s)
    % buil them all
    ind_s = 1:sdip.n_seeds
end
if nargout>1
    calcSV = 1;
end

Nind_s = length(ind_s);
Vfit = cell(1,Nind_s);
if calcSV
    SV = cell(1,Nind_s);
end
for ii=1:Nind_s
    ind_s_i = ind_s(ii);
    Vfit{ii} = sdip.L{ind_s_i}*sdip.j{ind_s_i};
    if calcSV
        SV{ii} = spm_eeg_CreateSV(Vfit{ii},[],varargin{:});
    end
    if dispSV
        spm_eeg_DrawSV(SV{ii},[]);
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [SV_kk,h] = drawSV(SV{16})
% 
% 
% dd = data(:,sdip.wind_be(1):sdip.wind_be(2))
% SVdd=CreateSV(dd)
% drawSV(SVdd)
% 
% di = dd-Vfit{16};
% SVdi=CreateSV(di)
% drawSV(SVdi)
% 
