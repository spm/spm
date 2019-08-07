function out = spm_sextract(job)
% Surface extraction
%__________________________________________________________________________
% Copyright (C) 2008-2016 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_sextract.m 7652 2019-08-07 11:30:35Z john $


images = job.images;
Vi     = spm_vol(char(images));

out.Surface = {};

for k=1:numel(job.surface)
    expression = job.surface(k).expression;
    thresh     = job.surface(k).thresh;
    
    y = zeros(Vi(1).dim(1:3),'single');
    spm_progress_bar('Init',Vi(1).dim(3),expression,'planes completed');
    for p  = 1:Vi(1).dim(3)
        B  = spm_matrix([0 0 -p 0 0 0 1 1 1]);
        im = cell(1,numel(Vi));
        for i=1:numel(im)
            M     = inv(B*inv(Vi(1).mat)*Vi(i).mat);
            im{i} = spm_slice_vol(Vi(i),M,Vi(1).dim(1:2),[0,NaN]);
        end
        try
            y(:,:,p) = real(single(efun(im,expression)));
        catch
            error(['Can not evaluate "' expression '".']);
        end
        spm_progress_bar('Set',p);
    end
    spm_smooth(y,y,[1.5,1.5,1.5]);
    
    [faces,vertices] = isosurface(y,thresh);
    
    if isempty(vertices), error('No surface.'); end
    
    % Swap around x and y because isosurface does for some
    % wierd and wonderful reason.
    Mat      = Vi(1).mat(1:3,:)*[0 1 0 0;1 0 0 0;0 0 1 0; 0 0 0 1];
    vertices = (Mat*[vertices' ; ones(1,size(vertices,1))])';
    giiname  = spm_file(images{1},...
        'prefix','surf_',...
        'suffix',sprintf('_%.3d',k),...
        'ext','gii', 'number','');
    save(gifti(struct('faces',faces,'vertices',vertices)),giiname);
    out.SurfaceFile{k} = giiname;
    spm_progress_bar('Clear');
end

%==========================================================================
function y = efun(im,f)
for i=1:numel(im)
    eval(['i' num2str(i) '= im{i};']);
end
y = eval(f);

