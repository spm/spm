function [P] = mci_lds_par2lat (Pt,M)
% Convert parmas to latent params
% FORMAT [P] = mci_lds_par2lat (Pt,M)
%
% Pt    params 
% M     model struct
%
% P     params (latent)

a=Pt(1:M.d)/M.a_typical;
a=log(a);
b=Pt(M.d+1:end);

P=[a(:);b(:)];