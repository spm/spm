function [P] = mci_lds_par2lat (Pt,M)
% Convert parmas to latent params
% FORMAT [P] = mci_lds_par2lat (Pt,M)
%
% Pt    params 
% M     model struct
%
% P     params (latent)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_lds_par2lat.m 6548 2015-09-11 12:39:47Z will $

a=Pt(1:M.d)/M.a_typical;
a=log(a);
b=Pt(M.d+1:end);

P=[a(:);b(:)];