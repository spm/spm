function render = tbx_cfg_render
% Configuration file for toolbox 'Rendering'
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: tbx_cfg_render.m 5010 2012-10-19 11:47:42Z john $

% ---------------------------------------------------------------------
% images Input Images
% ---------------------------------------------------------------------
images         = cfg_files;
images.tag     = 'images';
images.name    = 'Input Images';
images.help    = {'These are the images that are used by the calculator.  They are referred to as i1, i2, i3, etc in the order that they are specified.'};
images.filter  = 'image';
images.ufilter = '.*';
images.num     = [1 Inf];
% ---------------------------------------------------------------------
% expression Expression
% ---------------------------------------------------------------------
expression         = cfg_entry;
expression.tag     = 'expression';
expression.name    = 'Expression';
expression.help    = {
                      'Example expressions (f):'
                      '    * Mean of six images (select six images)'
                      '       f = ''(i1+i2+i3+i4+i5+i6)/6'''
                      '    * Make a binary mask image at threshold of 100'
                      '       f = ''i1>100'''
                      '    * Make a mask from one image and apply to another'
                      '       f = ''i2.*(i1>100)'''
                      '             - here the first image is used to make the mask, which is applied to the second image'
                      '    * Sum of n images'
                      '       f = ''i1 + i2 + i3 + i4 + i5 + ...'''
                      }';
expression.strtype = 's';
expression.num     = [2  Inf];
expression.val     = {'i1'};
% ---------------------------------------------------------------------
% thresh Surface isovalue(s)
% ---------------------------------------------------------------------
thresh         = cfg_entry;
thresh.tag     = 'thresh';
thresh.name    = 'Surface isovalue(s)';
thresh.help    = {'Enter the value at which isosurfaces through the resulting image is to be computed.'};
thresh.strtype = 'e';
thresh.num     = [1  1];
thresh.val     = {0.5};
% ---------------------------------------------------------------------
% surface Surface
% ---------------------------------------------------------------------
surface         = cfg_branch;
surface.tag     = 'surface';
surface.name    = 'Surface';
surface.val     = {expression thresh };
surface.help    = {'An expression and threshold for each of the surfaces to be generated.'};
% ---------------------------------------------------------------------
% Surfaces Surfaces
% ---------------------------------------------------------------------
Surfaces         = cfg_repeat;
Surfaces.tag     = 'Surfaces';
Surfaces.name    = 'Surfaces';
Surfaces.help    = {'Multiple surfaces can be created from the same image data.'};
Surfaces.values  = {surface };
Surfaces.num     = [0 Inf];
% ---------------------------------------------------------------------
% SExtract Surface Extraction
% ---------------------------------------------------------------------
SExtract         = cfg_exbranch;
SExtract.tag     = 'SExtract';
SExtract.name    = 'Surface Extraction';
SExtract.val     = {images Surfaces };
SExtract.help    = {'User-specified algebraic manipulations are performed on a set of images, with the result being used to generate a surface file. The user is prompted to supply images to work on and a number of expressions to evaluate, along with some thresholds. The expression should be a standard matlab expression, within which the images should be referred to as i1, i2, i3,... etc. An isosurface file is created from the results at the user-specified threshold.'};
SExtract.prog = @spm_local_sextract;
%SExtract.vfiles = @filessurf;
SExtract.vout = @vout_sextract;

% ---------------------------------------------------------------------
% SurfaceFile Surface File
% ---------------------------------------------------------------------
SurfaceFile         = cfg_files;
SurfaceFile.tag     = 'SurfaceFile';
SurfaceFile.name    = 'Surface File';
SurfaceFile.help    = {'Filename of the surf_*.gii file containing the rendering information. This can be generated via the surface extraction routine in SPM. Normally, a surface is extracted from grey and white matter tissue class images, but it is also possible to threshold e.g. an spmT image so that activations can be displayed.'};
SurfaceFile.filter = 'mesh';
SurfaceFile.ufilter = '.*';
SurfaceFile.num     = [1 1];
% ---------------------------------------------------------------------
% Red Red
% ---------------------------------------------------------------------
Red         = cfg_menu;
Red.tag     = 'Red';
Red.name    = 'Red';
Red.help    = {'The intensity of the red colouring (0 to 1).'};
Red.labels = {
              '0.0'
              '0.2'
              '0.4'
              '0.6'
              '0.8'
              '1.0'
              }';
Red.values = {
              0
              0.2
              0.4
              0.6
              0.8
              1
              }';
Red.val    = {1};
% ---------------------------------------------------------------------
% Green Green
% ---------------------------------------------------------------------
Green         = cfg_menu;
Green.tag     = 'Green';
Green.name    = 'Green';
Green.help    = {'The intensity of the green colouring (0 to 1).'};
Green.labels = {
                '0.0'
                '0.2'
                '0.4'
                '0.6'
                '0.8'
                '1.0'
                }';
Green.values = {
                0
                0.2
                0.4
                0.6
                0.8
                1
                }';
Green.val     = {1};
% ---------------------------------------------------------------------
% Blue Blue
% ---------------------------------------------------------------------
Blue         = cfg_menu;
Blue.tag     = 'Blue';
Blue.name    = 'Blue';
Blue.help    = {'The intensity of the blue colouring (0 to 1).'};
Blue.labels = {
               '0.0'
               '0.2'
               '0.4'
               '0.6'
               '0.8'
               '1.0'
               }';
Blue.values = {
               0
               0.2
               0.4
               0.6
               0.8
               1
               }';
Blue.val     = {1};
% ---------------------------------------------------------------------
% Color Color
% ---------------------------------------------------------------------
Color         = cfg_branch;
Color.tag     = 'Color';
Color.name    = 'Color';
Color.val     = {Red Green Blue };
Color.help    = {'Specify the colour using a mixture of red, green and blue. For example, white is specified by 1,1,1, black is by 0,0,0 and purple by 1,0,1.'};
% ---------------------------------------------------------------------
% DiffuseStrength Diffuse Strength
% ---------------------------------------------------------------------
DiffuseStrength         = cfg_menu;
DiffuseStrength.tag     = 'DiffuseStrength';
DiffuseStrength.name    = 'Diffuse Strength';
DiffuseStrength.help    = {'The strength with which the object diffusely reflects light. Mat surfaces reflect light diffusely, whereas shiny surfaces reflect speculatively.'};
DiffuseStrength.labels = {
                          '0.0'
                          '0.2'
                          '0.4'
                          '0.6'
                          '0.8'
                          '1.0'
                          }';
DiffuseStrength.values = {
                          0
                          0.2
                          0.4
                          0.6
                          0.8
                          1
                          }';
DiffuseStrength.val    = {0.8};
% ---------------------------------------------------------------------
% AmbientStrength Ambient Strength
% ---------------------------------------------------------------------
AmbientStrength         = cfg_menu;
AmbientStrength.tag     = 'AmbientStrength';
AmbientStrength.name    = 'Ambient Strength';
AmbientStrength.help    = {'The strength with which the object reflects ambient (non-directional) lighting.'};
AmbientStrength.labels = {
                          '0.0'
                          '0.2'
                          '0.4'
                          '0.6'
                          '0.8'
                          '1.0'
                          }';
AmbientStrength.values = {
                          0
                          0.2
                          0.4
                          0.6
                          0.8
                          1
                          }';
AmbientStrength.val    = {0.2};
% ---------------------------------------------------------------------
% SpecularStrength Specular Strength
% ---------------------------------------------------------------------
SpecularStrength         = cfg_menu;
SpecularStrength.tag     = 'SpecularStrength';
SpecularStrength.name    = 'Specular Strength';
SpecularStrength.help    = {'The strength with which the object specularly reflects light (i.e. how shiny it is). Mat surfaces reflect light diffusely, whereas shiny surfaces reflect speculatively.'};
SpecularStrength.labels = {
                           '0.0'
                           '0.2'
                           '0.4'
                           '0.6'
                           '0.8'
                           '1.0'
                           }';
SpecularStrength.values = {
                           0
                           0.2
                           0.4
                           0.6
                           0.8
                           1
                           }';
SpecularStrength.val    = {0.2};
% ---------------------------------------------------------------------
% SpecularExponent Specular Exponent
% ---------------------------------------------------------------------
SpecularExponent         = cfg_menu;
SpecularExponent.tag     = 'SpecularExponent';
SpecularExponent.name    = 'Specular Exponent';
SpecularExponent.help    = {'A parameter describing the specular reflectance behaviour. It relates to the size of the high-lights.'};
SpecularExponent.labels = {
                           '0.01'
                           '0.1'
                           '10'
                           '100'
                           }';
SpecularExponent.values = {
                           0.01
                           0.1
                           10
                           100
                           }';
SpecularExponent.val    = {10};
% ---------------------------------------------------------------------
% SpecularColorReflectance Specular Color Reflectance
% ---------------------------------------------------------------------
SpecularColorReflectance         = cfg_menu;
SpecularColorReflectance.tag     = 'SpecularColorReflectance';
SpecularColorReflectance.name    = 'Specular Color Reflectance';
SpecularColorReflectance.help    = {'Another parameter describing the specular reflectance behaviour.'};
SpecularColorReflectance.labels = {
                                   '0.0'
                                   '0.2'
                                   '0.4'
                                   '0.6'
                                   '0.8'
                                   '1.0'
                                   }';
SpecularColorReflectance.values = {
                                   0
                                   0.2
                                   0.4
                                   0.6
                                   0.8
                                   1
                                   }';
SpecularColorReflectance.val    = {0.8};
% ---------------------------------------------------------------------
% FaceAlpha Face Alpha
% ---------------------------------------------------------------------
FaceAlpha         = cfg_menu;
FaceAlpha.tag     = 'FaceAlpha';
FaceAlpha.name    = 'Face Alpha';
FaceAlpha.help    = {'The opaqueness of the surface.  A value of 1 means it is opaque, whereas a value of 0 means it is transparent.'};
FaceAlpha.labels = {
                    '0.0'
                    '0.2'
                    '0.4'
                    '0.6'
                    '0.8'
                    '1.0'
                    }';
FaceAlpha.values = {
                    0
                    0.2
                    0.4
                    0.6
                    0.8
                    1
                    }';
FaceAlpha.val    = {1};
% ---------------------------------------------------------------------
% Object Object
% ---------------------------------------------------------------------
Object         = cfg_branch;
Object.tag     = 'Object';
Object.name    = 'Object';
Object.val     = {SurfaceFile Color DiffuseStrength AmbientStrength SpecularStrength SpecularExponent SpecularColorReflectance FaceAlpha };
Object.help    = {'Each object is a surface (from a surf_*.gii file), which may have a number of light-reflecting qualities, such as colour and shinyness.'};
% ---------------------------------------------------------------------
% Objects Objects
% ---------------------------------------------------------------------
Objects         = cfg_repeat;
Objects.tag     = 'Objects';
Objects.name    = 'Objects';
Objects.help    = {'Several surface objects can be displayed together in different colours and with different reflective properties.'};
Objects.values  = {Object };
Objects.num     = [0 Inf];
% ---------------------------------------------------------------------
% Position Position
% ---------------------------------------------------------------------
Position         = cfg_entry;
Position.tag     = 'Position';
Position.name    = 'Position';
Position.help    = {'The position of the light in 3D.'};
Position.strtype = 'e';
Position.num     = [1  3];
Position.val     = {[100 100 100]};
% ---------------------------------------------------------------------
% Red Red
% ---------------------------------------------------------------------
Red         = cfg_menu;
Red.tag     = 'Red';
Red.name    = 'Red';
Red.help    = {'The intensity of the red colouring (0 to 1).'};
Red.labels = {
              '0.0'
              '0.2'
              '0.4'
              '0.6'
              '0.8'
              '1.0'
              }';
Red.values = {
              0
              0.2
              0.4
              0.6
              0.8
              1
              }';
Red.val    = {1};
% ---------------------------------------------------------------------
% Green Green
% ---------------------------------------------------------------------
Green         = cfg_menu;
Green.tag     = 'Green';
Green.name    = 'Green';
Green.help    = {'The intensity of the green colouring (0 to 1).'};
Green.labels = {
                '0.0'
                '0.2'
                '0.4'
                '0.6'
                '0.8'
                '1.0'
                }';
Green.values = {
                0
                0.2
                0.4
                0.6
                0.8
                1
                }';
Green.val    = {1};
% ---------------------------------------------------------------------
% Blue Blue
% ---------------------------------------------------------------------
Blue         = cfg_menu;
Blue.tag     = 'Blue';
Blue.name    = 'Blue';
Blue.help    = {'The intensity of the blue colouring (0 to 1).'};
Blue.labels = {
               '0.0'
               '0.2'
               '0.4'
               '0.6'
               '0.8'
               '1.0'
               }';
Blue.values = {
               0
               0.2
               0.4
               0.6
               0.8
               1
               }';
Blue.val    = {1};
% ---------------------------------------------------------------------
% Color Color
% ---------------------------------------------------------------------
Color         = cfg_branch;
Color.tag     = 'Color';
Color.name    = 'Color';
Color.val     = {Red Green Blue };
Color.help    = {'Specify the colour using a mixture of red, green and blue. For example, white is specified by 1,1,1, black is by 0,0,0 and purple by 1,0,1.'};
% ---------------------------------------------------------------------
% Light Light
% ---------------------------------------------------------------------
Light         = cfg_branch;
Light.tag     = 'Light';
Light.name    = 'Light';
Light.val     = {Position Color };
Light.help    = {'Specification of a light source in terms of position and colour.'};
% ---------------------------------------------------------------------
% Lights Lights
% ---------------------------------------------------------------------
Lights         = cfg_repeat;
Lights.tag     = 'Lights';
Lights.name    = 'Lights';
Lights.help    = {'There should be at least one light specified so that the objects can be clearly seen.'};
Lights.values  = {Light };
Lights.num     = [0 Inf];
% ---------------------------------------------------------------------
% SRender Surface Rendering
% ---------------------------------------------------------------------
SRender         = cfg_exbranch;
SRender.tag     = 'SRender';
SRender.name    = 'Surface Rendering';
SRender.val     = {Objects Lights };
SRender.help    = {'This utility is for visualising surfaces.  Surfaces first need to be extracted and saved in surf_*.gii files using the surface extraction routine.'};
SRender.prog = @spm_local_srender;
% ---------------------------------------------------------------------
% render Rendering
% ---------------------------------------------------------------------
render         = cfg_choice;
render.tag     = 'render';
render.name    = 'Rendering';
render.help    = {'This is a toolbox that provides a limited range of surface rendering options. The idea is to first extract surfaces from image data, which are saved in rend_*.mat files. These can then be loaded and displayed as surfaces. Note that OpenGL rendering is used, which can be problematic on some computers. The tools are limited - and they do what they do.'};
render.values  = {SExtract SRender };
%render.num     = [0 Inf];

%======================================================================
function spm_local_srender(job)
if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','SRender')); end
spm_srender(job);

%======================================================================
function out=spm_local_sextract(job)
if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','SRender')); end
out=spm_sextract(job);

%======================================================================
function dep = vout_sextract(job)
dep = cfg_dep;
for k=1:numel(job.surface),
    dep(k)            = cfg_dep;
    dep(k).sname      = ['Surface File ' num2str(k)];
    dep(k).src_output = substruct('.','SurfaceFile','()',{k});
%   dep(k).tgt_spec   = cfg_findspec({{'filter','.*\.gii$'}});
    dep(k).tgt_spec   = cfg_findspec({{'filter','mesh'}});
end

