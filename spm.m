
% SPM startup
% FORMAT spm
%_______________________________________________________________________
% spm prompts for PET or fMRI and calls the apropriate user interface
% routine
%
% see also spm_ui.m and spm_fmri_ui.m
%
%_______________________________________________________________________
% %W% %E%

%-Clean up workspace, print welcome banner
%-----------------------------------------------------------------------
close all
clear all
clc
disp(' ___  ____  __  __                                                  ')
disp('/ __)(  _ \(  \/  )  Statistical Parametric Mapping                 ')
disp('\__ \ )___/ )    (   The Wellcome Department of Cognitive Neurology ')
disp('(___/(__)  (_/\/\_)  Version SPMdevel:96                            ')
fprintf('\n')

%-Open startup window, set window defaults
%-----------------------------------------------------------------------
whitebg(0,'w')

figure('Color',[1 1 1]*.8,'Name',' ',...
	'NumberTitle','off','Position',[360 200 400 280],'Resize','off');

set(0,'DefaultFigureColormap',gray);
set(0,'DefaultFigurePaperType','a4letter')

%-Frames and text
%-----------------------------------------------------------------------
uicontrol(1,'Style','Frame','Position',[10 120 380 140]);
uicontrol(1,'Style','Frame','Position',[10 020 380 090]);

c = 'STATISTICAL PARAMETRIC MAPPING  -  {SPM95}';
uicontrol(1,'Style','Text', 'Position',[10 220 380 40],...
'String',c,'ForegroundColor',[1 0 0])

c = 'The Wellcome Department of Cognitive Neurology,';
uicontrol(1,'Style','Text', 'Position',[10 200 380 16],'String',c)
c = 'The Institute of Neurology and';
uicontrol(1,'Style','Text', 'Position',[10 180 380 16],'String',c)
c = 'The MRC Cyclotron Unit,';
uicontrol(1,'Style','Text', 'Position',[10 160 380 16],'String',c)
c = 'Hammersmith Hospital, London UK';
uicontrol(1,'Style','Text', 'Position',[10 140 380 16],'String',c)

global MODALITY
MODALITY = 'unknown';

%-Objects with Callbacks - main spm_*_ui.m routines
%-----------------------------------------------------------------------
uicontrol(1,'String','PET and SPECT',         'Position',[040 066 150 30],...
'CallBack','global MODALITY ; MODALITY = ''PET''; spm_ui','Interruptible','yes','ForegroundColor',[0 1 1]);

uicontrol(1,'String','fMRI time-series',      'Position',[210 066 150 30],...
'CallBack','global MODALITY ; MODALITY = ''FMRI''; spm_fmri_ui','Interruptible','yes','ForegroundColor',[0 1 1]);

uicontrol(1,'String','notes and bibliography','Position',[040 030 320 30],...
'CallBack','spm_bib','Interruptible','yes','ForegroundColor',[0 1 1]);
