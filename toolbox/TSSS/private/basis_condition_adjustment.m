% [SN_new,sss_indices] = basis_condition_adjustment(SN,dim_in,cond_threshold)
%
% Improve the condition number of the SSS basis by dropping unnecessary
% vectors from the basis. Input argument SN is the original normalized SSS
% basis, Lin is the order number of the internal basis, and cond_threshold
% is the acceptable condition number of the basis. 
%
% Copyright (c) 2016, Elekta Oy
% ---------------------------------------
% 
% Redistribution and use of the Software in source and binary forms, with or without 
% modification, are permitted for non-commercial use.
% 
% The Software is provided "as is" without warranties of any kind, either express or
% implied including, without limitation, warranties that the Software is free of defects,
% merchantable, fit for a particular purpose. Developer/user agrees to bear the entire risk 
% in connection with its use and distribution of any and all parts of the Software under this license.
% 

%
function [SN_new,sss_indices, n] = basis_condition_adjustment(SN, dim_in, cond_threshold)

%dim_in = (Lin + 1)^2 - 1;
n_basis = size(SN,2);

sss_cond = cond(SN);
sss_indices_in = 1:dim_in;
SNout = SN(:,dim_in+1:end);

while sss_cond > cond_threshold
    for j = 1:length(sss_indices_in)
        SN2 = [SN(:,setdiff(sss_indices_in,sss_indices_in(j))) SNout];
        c(j) = cond(SN2);
    end
    [sss_cond,drop] = min(c);
    sss_indices_in = setdiff(sss_indices_in,sss_indices_in(drop));
end
SN_new = [SN(:,sss_indices_in) SNout];
sss_indices = [sss_indices_in dim_in+1:n_basis];
n = length(sss_indices_in);

