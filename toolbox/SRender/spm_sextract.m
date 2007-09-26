function spm_sextract(job)
% Surface extraction
%_______________________________________________________________________
% Copyright (C) 2007 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

images = job.images;
Vi     = spm_vol(strvcat(images));
n      = numel(Vi);                %-#images
if n==0, error('no input images specified'), end

[pth,nam,ext] = fileparts(images{1});

for k=1:numel(job.surface),
    expression = job.surface(k).expression;
    thresh     = job.surface(k).thresh;

    y = zeros(Vi(1).dim(1:3),'single');
    spm_progress_bar('Init',Vi(1).dim(3),expression,'planes completed');
    for p  = 1:Vi(1).dim(3),
        B  = spm_matrix([0 0 -p 0 0 0 1 1 1]);
        im = cell(1,n);
        for i=1:n,
            M     = inv(B*inv(Vi(1).mat)*Vi(i).mat);
            im{i} = spm_slice_vol(Vi(i),M,Vi(1).dim(1:2),[0,NaN]);
        end
        try,
            y(:,:,p) = real(single(efun(im,expression)));
        catch,
            error(['Can not evaluate "' expression '".']);
        end
        spm_progress_bar('Set',p);
    end
    spm_smooth(y,y,[1.5,1.5,1.5]);

    [faces,vertices] = isosurface(y,thresh);
    % Swap around x and y because isosurface does for some
    % wierd and wonderful reason.
    Mat      = Vi(1).mat(1:3,:)*[0 1 0 0;1 0 0 0;0 0 1 0; 0 0 0 1];
    vertices = (Mat*[vertices' ; ones(1,size(vertices,1))])';
    matname  = fullfile(pth,sprintf('surf_%s_%.3d.mat',nam,k));
    if spm_matlab_version_chk('7.0') >=0,
        save(matname,'-V6','faces','vertices','expression','thresh','images');
    else
        save(matname,'faces','vertices','expression','thresh','images');
    end;
    spm_progress_bar('Clear');
end
return
%_______________________________________________________________________
%_______________________________________________________________________
function y = efun(im,f)
for i=1:numel(im),
    eval(['i' num2str(i) '= im{i};']);
end
y = eval(f);
return

