function varargout=spm_XYZreg_Ex2(varargin)
% Example of Registry enabled XYZ GUI control / function
% FORMAT...
%_______________________________________________________________________
%
% Help goes here...
%
% Object must be indentifiable via a unique HandleGraphics object.
% In this code, this handle is called hMe.
%
% This HandleGraphics objects 'UserData' *must* be a structure.
% The structure must have a field called 'hReg', which stores the handle
% of the registry when linked, and is empty when not. Some utility features
% of spm_XYZreg will set/delete the handle directly...
%
% There must be a 'SetCoords' function for this object, with call:
%   spm_res_ui('SetCoords',xyz,hMe,hC)
% ...this can handle interna, co-ordinate setting (as in this example), but
% must also call the registry.
%
% The registry update function is:
%   spm_XYZreg('SetCoords',xyz,hReg,hMe);
% ...which must be called at all points where the local co-ordinates can be
% changed. It is robust to invalid or empty hReg.
%
% It's *vital* to specify caller handles (hC), so that the registry doesn't
% end up in an infinite loop of updating!
%
% Hey, if your function has multiple places where you can change the XYZ,
% you could use an ``internal'' registry locally, with the external registry
% as one of it's entries! (I think?)
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_XYZreg_Ex2.m 5219 2013-01-29 17:07:07Z spm $


%=======================================================================
switch lower(varargin{1}), case 'create'
%=======================================================================
% hMe = spm_XYZreg_Ex2('Create',M,D,xyz)
if nargin<4, xyz=[0;0;0]; else, xyz=varargin{4}; end
if nargin<3, error('Insufficient arguments'), end
D = varargin{3};
M = varargin{2};

xyz = spm_XYZreg('RoundCoords',xyz,M,D);

F = figure;

%-Create control:
% - note the UserData structure: hReg is necessary for cross-registration,
%   M, D, & xyz are used internally...
% - note how the DeleteFcn closes the objects figure...
hMe = uicontrol(F,'Style','Pushbutton',...
    'String',sprintf('[%.3f, %.3f, %.3f]',xyz),...
    'Position',[100 100 300 30],...
    'FontSize',14,'FontWeight','Bold',...
    'Callback',...
    'spm_XYZreg_Ex2(''SetCoords'',input(''Enter [x,y,z]'''' vector: ''),gcbo);',...
    'UserData',struct(...
        'hReg', [],...
        'M',    M,...
        'D',    D,...
        'xyz',  xyz ),...
    'DeleteFcn','delete(gcbf)');

varargout = {hMe};



%=======================================================================
case 'setcoords'    % Set co-ordinates
%=======================================================================
% [xyz,d] = spm_XYZreg_Ex2('SetCoords',xyz,hMe,hC)
if nargin<4, hC=0; else hC=varargin{4}; end
if nargin<3, error('Insufficient arguments'), else hMe=varargin{3}; end
%-Or possibly some clever code to guess the handle
if nargin<2, error('Set co-ords to what!'), else xyz=varargin{2}; end

UD = get(hMe,'UserData');

%-Check validity of coords only when called without a caller handle
%-----------------------------------------------------------------------
if hC<=0
    [xyz,d] = spm_XYZreg('RoundCoords',xyz,UD.M,UD.D);
    if d>0 && nargout<2, warning(sprintf(...
        '%s: Co-ords rounded to nearest voxel center: Discrepancy %.2f',...
        mfilename,d)), end
else
    d = [];
end

%-Update local XYZ information
%-----------------------------------------------------------------------
UD.xyz = xyz;
set(hMe,'UserData',UD)
set(hMe,'String',sprintf('[%.3f, %.3f, %.3f]',xyz))

%-Tell the registry, if we've not been called by the registry...
%-----------------------------------------------------------------------
if (~isempty(UD.hReg) && UD.hReg~=hC)
    spm_XYZreg('SetCoords',xyz,UD.hReg,hMe);
end

%-Return arguments
%-----------------------------------------------------------------------
varargout = {xyz,d};



%=======================================================================
otherwise
%=======================================================================
error('Unknown action string')

%=======================================================================
end
