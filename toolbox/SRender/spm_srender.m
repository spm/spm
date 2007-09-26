function spm_srender(job)
% A function for rendering surfaces
% FORMAT spm_render(job)
% job - a job structure (see spm_config_render.m)
%_______________________________________________________________________
% Copyright (C) 2007 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$


fg  = spm_figure('GetWin','Graphics');
ren = get(fg,'Renderer');
clf(fg);
set(fg,'Renderer','OpenGL');
ax = axes('Parent',fg,'DeleteFcn',['rotate3d off; set(gcf,''Renderer'',''' ren ''');']);

for i=1:numel(job.Object),
     obj = job.Object(i);
     for j=1:numel(obj.SurfaceFile),
         FVo = load(obj.SurfaceFile{j});
         FV  = struct('faces',FVo.faces,'vertices',FVo.vertices);
         p  = patch(FV, 'Parent',ax,...
             'FaceColor', [obj.Color.Red,obj.Color.Green, obj.Color.Blue],...
             'FaceVertexCData', [],...
             'EdgeColor', 'none',...
             'FaceLighting', 'phong',...
             'SpecularStrength' ,obj.SpecularStrength,...
             'AmbientStrength', obj.AmbientStrength,...
             'DiffuseStrength', obj.DiffuseStrength,...
             'SpecularExponent', obj.SpecularExponent,...
             'FaceAlpha',obj.FaceAlpha);
    end
end
for i=1:numel(job.Light),
    obj = job.Light(i);
    l   = light('Parent',ax,...
        'Position',obj.Position,...
        'Color',[obj.Color.Red,obj.Color.Green, obj.Color.Blue]);
end

set(0,'CurrentFigure',fg);
set(fg,'CurrentAxes',ax);
axis image;
rotate3d on;
drawnow;

