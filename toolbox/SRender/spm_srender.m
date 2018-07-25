function spm_srender(job)
% A function for rendering surfaces
% FORMAT spm_srender(job)
% job - a job structure (see tbx_cfg_render.m)
%__________________________________________________________________________
% Copyright (C) 2008-2018 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_srender.m 7381 2018-07-25 10:27:54Z guillaume $


fg  = spm_figure('GetWin','Graphics');
ren = get(fg,'Renderer');
clf(fg);
set(fg,'Renderer','OpenGL');
ax = axes('Parent',fg,'DeleteFcn',['rotate3d off; set(gcf,''Renderer'',''' ren ''');']);

try
    set(fg,'CurrentAxes',ax);
    cameratoolbar(fg);
    drawnow;
end

for i=1:numel(job.Object)
     obj = job.Object(i);
     for j=1:numel(obj.SurfaceFile)
         FVo = struct(gifti(obj.SurfaceFile{j}));
         FV  = struct('faces',FVo.faces,'vertices',FVo.vertices);
         p   = patch(FV, 'Parent',ax,...
             'FaceColor', [obj.Color.Red,obj.Color.Green, obj.Color.Blue],...
             'FaceVertexCData',  [],...
             'EdgeColor',        'none',...
             'FaceLighting',     'gouraud',...
             'EdgeLighting',     'gouraud',...
             'SpecularStrength', obj.SpecularStrength,...
             'AmbientStrength',  obj.AmbientStrength,...
             'DiffuseStrength',  obj.DiffuseStrength,...
             'SpecularExponent', obj.SpecularExponent,...
             'FaceAlpha',        obj.FaceAlpha);
    end
end
for i=1:numel(job.Light)
    obj = job.Light(i);
    l   = light('Parent',ax,...
        'Position',obj.Position,...
        'Color',[obj.Color.Red,obj.Color.Green, obj.Color.Blue]);
end

%set(0,'CurrentFigure',fg);
set(fg,'CurrentAxes',ax);
axis image equal off;
drawnow;
