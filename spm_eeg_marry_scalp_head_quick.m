function [col_scalp]=marry_scalp_head_quick(V,C,V3,elec_verts,dv,id);
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% James Kilner
% $Id: spm_eeg_marry_scalp_head_quick.m 112 2005-05-04 18:20:52Z john $

col_scalp=zeros(length(dv),1);
P = fullfile(spm('dir'), 'EEGtemplates');
load(fullfile(P, 'params2.mat'));
col_scalp(col_scalp==0)=NaN;
allid(allid(:,1)==0,:)='';
alltds(alltds(:,1)==0,:)='';
alli(alli(:,1)==0,:)='';

col_scalp(allid)=sum((alltds(:,:).*C(alli(:,1:10))),2);
