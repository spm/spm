function varargout = spm_eeg_inv_vde(varargin)
% SPM_EEG_INV_VDE M-file for spm_eeg_inv_vde.fig
%      SPM_EEG_INV_VDE, by itself, creates a new SPM_EEG_INV_VDE or raises the existing
%      singleton*.
%
%      H = SPM_EEG_INV_VDE returns the handle to a new SPM_EEG_INV_VDE or the handle to
%      the existing singleton*.
%
%      SPM_EEG_INV_VDE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPM_EEG_INV_VDE.M with the given input arguments.
%
%      SPM_EEG_INV_VDE('Property','Value',...) creates a new SPM_EEG_INV_VDE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spm_eeg_inv_vde_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spm_eeg_inv_vde_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help spm_eeg_inv_vde

% Last Modified by GUIDE v2.5 05-Dec-2006 09:07:07

% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
% Jeremie Mattout
% $Id: spm_eeg_inv_vde.m 1039 2007-12-21 20:20:38Z karl $

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spm_eeg_inv_vde_OpeningFcn, ...
                   'gui_OutputFcn',  @spm_eeg_inv_vde_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before spm_eeg_inv_vde is made visible.
function spm_eeg_inv_vde_OpeningFcn(hObject, eventdata, handles, varargin)
handles.main_handles = varargin{1};
set(handles.Location,'Enable','on');
set(handles.Display,'Enable','off');
handles.output = hObject;
guidata(hObject,handles);

% --- Outputs from this function are returned to the command line.
function varargout = spm_eeg_inv_vde_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Executes on button press in Location.
function Location_Callback(hObject, eventdata, handles)
axes(handles.main_handles.sources_axes);
vert = handles.main_handles.vert;
handles.location = datacursormode(handles.main_handles.figure1);
set(handles.location,'Enable','on','DisplayStyle','datatip','SnapToDataVertex','off');
waitforbuttonpress;
vde = getCursorInfo(handles.location);
vde = vde.Position;
datacursormode off
CurrVert = sum(((ones(length(vert),1)*vde - vert).^2)');
handles.vde = find(CurrVert == min(CurrVert));
set(handles.Location,'Enable','off');
set(handles.Display,'Enable','on');
guidata(hObject,handles);


% --- Executes on button press in Size.
function Size_Callback(hObject, eventdata, handles)
axes(handles.main_handles.sources_axes);
vert = handles.main_handles.vert;
face = handles.main_handles.face;
Neighbours = FindNeighB(handles.vde,face);
amp = handles.main_handles.srcs_disp;
ampref  = amp(handles.vde(1));
ampnew  = amp(Neighbours);
amptemp = abs(ampnew - ampref);
newsrc = find(amptemp == min(amptemp));
handles.vde = [handles.vde Neighbours(newsrc)];
axes(handles.main_handles.sources_axes);
hold on;
handles.hpts = plot3(vert(handles.vde,1),vert(handles.vde,2),vert(handles.vde,3),'sk','MarkerFaceColor','k','MarkerSize',14);
guidata(hObject,handles);

% --- Executes on button press in Display
function Display_Callback(hObject, eventdata, handles)
figure;
if length(handles.vde) > 1
   TimeSeries = mean(handles.main_handles.srcs_data(handles.vde,:)); 
else
    TimeSeries = handles.main_handles.srcs_data(handles.vde,:);
end
D = handles.main_handles.D;
woi = D.inv{handles.main_handles.val}.inverse.woi;
Xstep = (woi(2) - woi(1))/(handles.main_handles.dimT - 1);
X     = woi(1):Xstep:woi(2);
plot(X,TimeSeries,'r','LineWidth',3);
set(handles.hpts,'Visible','off');
hold off;
close(handles.figure1);
axes(handles.main_handles.sources_axes)
cameramenu on

function Neighbours = FindNeighB(CurrVert,face)
Neighbours = [];
for i = 1:length(CurrVert)
    [Icv,Jcv] = find(face == CurrVert(i));
    CurrNeigh = [];
    for j = 1:length(Icv)
        CurrNeigh = [CurrNeigh face(Icv(j),:)];
    end
    CurrNeigh = unique(CurrNeigh);
    Neighbours = [Neighbours CurrNeigh];
end
Neighbours = setdiff(Neighbours,CurrVert);
return

function pushbutton1_Callback(hObject,eventdata,handles)