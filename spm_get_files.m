function Q = spm_get_files(dir,fil,uni)

[P,tmp] = spm_list_files(dir,fil);
for k = 1:size(P,1)
  Q(k,:) = fullfile(dir,P(k,:));
end

