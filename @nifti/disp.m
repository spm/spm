function disp(obj)
% Disp a NIFTI-1 object
% _______________________________________________________________________
% @(#)disp.m	1.1 John Ashburner 04/11/26

sz = size(obj);
fprintf('NIFTI object: ');
if length(sz)>4,
    fprintf('%d-D\n',length(sz));
else
    for i=1:(length(sz)-1),
        fprintf('%d-by-',sz(i));
    end;
    fprintf('%d\n',sz(end));
end;
if prod(sz)==1,
    display(structn(obj))
end;
return;
