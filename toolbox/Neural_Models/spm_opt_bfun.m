function spm_opt_bfun


d       = round(get(gca,'CurrentPoint'));
i       = d(1);
pE      = get(gcf,'UserData');
P       = spm_fieldindices(pE,i);
j       = spm_fieldindices(pE,P);
j       = find(j == i);
[i j k] = ind2sub(size(spm_cat(getfield(pE,P))),j);


fprintf('\n%s(%i,%i,%i)\n',P,i,j,k);

return
