function [cfg, artifact] = artifact_manual(cfg);

% ARTIFACT_MANUAL allows the user to detect artifacts manually using visual
% inspection.
%
% You must specify the channels that should be displayed:
%
%   cfg.artfctdef.manual.channel        = cell-array with channels to be displayed.
% Be careful not to specify to much channels
% because the function then will be very slow.
%
% You can specify the following configuration options:
%
%   cfg.artfctdef.manual.pretrialtime   = time shown before trialstart (default 0)
%   cfg.artfctdef.manual.posttrialtime  = time shown after trialend    (default 0)
%   cfg.artfctdef.manual.fft	        = 'no' (default) or 'yes' turns on FFT window
%   cfg.artfctdef.manual.padding        = 'no' (default) or FFT-padding in seconds
%
% The FFT will be executed on the complete time interval including pre and post
% trial times.
%
%   cfg.artfctdef.manual.timeaxrelative = 'yes' (default) or 'no'.
%
% Set to yes defines the time axes relative to trialstart. Set to no defines it
% relative to the beginning of the experiment.
%
%   cfg.artfctdef.manual.blc       = 'no' (default) or 'yes' apply baseline correction
%   cfg.artfctdef.manual.bpfilter  = 'no' (default) or 'yes' apply bandpass filter
%   cfg.artfctdef.manual.bpfreq    = [0.3 30] in Hz
%   cfg.artfctdef.manual.bpfiltord = 2
%
% See also REJECTARTIFACT

% Undocumented local options:
% cfg.trl
% cfg.version
% cfg.artfctdef.manual.maxnumberofchannels = 20 (default)
%
% This function depends on PREPROC which has the following options:
% cfg.absdiff
% cfg.blc in ARTIFACT_MANUAL as cfg.artfctdef.manual.blc, documented
% cfg.blcwindow
% cfg.boxcar
% cfg.bpfilter in ARTIFACT_MANUAL as cfg.artfctdef.manual.bpfilter, documented
% cfg.bpfiltord in ARTIFACT_MANUAL as cfg.artfctdef.manual.bpfiltord, documented
% cfg.bpfilttype
% cfg.bpfreq in ARTIFACT_MANUAL as cfg.artfctdef.manual.bpfreq, documented
% cfg.derivative
% cfg.detrend
% cfg.dftfilter
% cfg.dftfreq
% cfg.hilbert
% cfg.hpfilter
% cfg.hpfiltord
% cfg.hpfilttype
% cfg.hpfreq
% cfg.implicitref
% cfg.lnfilter
% cfg.lnfiltord
% cfg.lnfreq
% cfg.lpfilter
% cfg.lpfiltord
% cfg.lpfilttype
% cfg.lpfreq
% cfg.medianfilter
% cfg.medianfiltord
% cfg.rectify
% cfg.refchannel
% cfg.reref

% Copyright (C) 2004, Geerten Kramer, FCDC
% version 23-11-2004
%
% $Log: artifact_manual.m,v $
% Revision 1.14  2006/11/29 09:06:36  roboos
% renamed all cfg options with "sgn" into "channel", added backward compatibility when required
% updated documentation, mainly in the artifact detection routines
%
% Revision 1.13  2006/06/14 12:43:51  roboos
% removed the documentation for cfg.lnfilttype, since that option is not supported by preproc
%
% Revision 1.12  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.11  2006/04/19 09:41:07  ingnie
% updated documentation
%
% Revision 1.10  2006/02/07 08:20:24  roboos
% fixed silly bug that was introduced along with dataset2files
%
% Revision 1.9  2006/01/31 13:49:29  jansch
% included dataset2files to ensure the presence of cfg.headerfile or cfg.datafile
% whenever needed
%
% Revision 1.8  2005/12/20 08:36:47  roboos
% add the artifact Nx2 matrix to the output configuration
% changed some indentation and white space, renamed a few variables
%
% Revision 1.7  2005/10/31 13:02:36  geekra
% - Made the field cfg.artfctdef.manual.channel required and removed the default
%   value 'all'. This is done because 'all' selects more channels than the
%   function can handle.
% - Build in a check on the number of channels selected, resulting in an error
%   if: 1. No channels are selected or 2. more then 20 channels are selected.
% - Added the non documented option cfg.artfctdef.manual.maxnumberofchannels
% - Changed some Dutch comment into English.
%
% Revision 1.6  2005/06/29 12:42:00  roboos
% added version to the output configuration
%
% Revision 1.5  2005/05/17 17:50:36  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.4  2005/03/02 14:28:44  roboos
% added backward compatibility support for xxx.sgn (copy into xxx.channel)
%
% Revision 1.3  2005/01/25 13:15:41  roboos
% added check for cfg.datatype=continuous and extended the call to read_fcdc_data with the boundary check for non-continuous data
%
% Revision 1.2  2004/12/20 12:22:56  roboos
% fixed multiple bugs, added general preproc function, cleaned up code
%
% Revision 1.1  2004/12/20 08:47:50  roboos
% initial implementation, based on Geertens version of 23 November 2004
%

% set default parameters if necessary.
if ~isfield(cfg, 'artfctdef'),                          cfg.artfctdef                            = [];       end
if ~isfield(cfg.artfctdef,'manual'),                    cfg.artfctdef.manual                     = [];       end
if ~isfield(cfg.artfctdef.manual,'zscale'),             cfg.artfctdef.manual.zscale              = 100;      end
if ~isfield(cfg.artfctdef.manual,'fft'),                cfg.artfctdef.manual.fft                 ='no';      end
if ~isfield(cfg.artfctdef.manual,'padding'),            cfg.artfctdef.manual.padding             ='no';      end
if ~isfield(cfg.artfctdef.manual,'pretrialtime'),       cfg.artfctdef.manual.pretrialtime        = 0;        end
if ~isfield(cfg.artfctdef.manual,'posttrialtime'),      cfg.artfctdef.manual.posttrialtime       = 0;        end
if ~isfield(cfg.artfctdef.manual,'timeaxrelative'),     cfg.artfctdef.manual.timeaxrelative      = 'yes';    end
if ~isfield(cfg.artfctdef.manual,'maxnumberofchannels'),cfg.artfctdef.manual.maxnumberofchannels = 20;       end

% for backward compatibility
if isfield(cfg.artfctdef.manual,'sgn')
  cfg.artfctdef.manual.channel = cfg.artfctdef.manual.sgn;
  cfg.artfctdef.manual         = rmfield(cfg.artfctdef.manual, 'sgn');
end

if ~isfield(cfg.artfctdef.manual,'channel'),
  % set an unusual default because all crashes the program.
  cfg.artfctdef.manual.channel = [];
end

if ~isfield(cfg, 'datatype') || ~strcmp(cfg.datatype, 'continuous')
  % datatype is unknown or not continuous, perform epoch boundary check
  iscontinuous = 0;
else
  % do not perform epoch boundary check, usefull for pseudo-continuous data
  iscontinuous = strcmp(cfg.datatype, 'continuous');
end

% read the header and do some preprocessing on the configuration
fprintf('Reading raw data...');
cfg = dataset2files(cfg);
hdr = read_fcdc_header(cfg.headerfile);
cfg.artfctdef.manual.channel=channelselection(cfg.artfctdef.manual.channel, hdr.label);
cfg.artfctdef.manual.trl=cfg.trl;
if(isempty(cfg.artfctdef.manual.channel))
	error(sprintf('\nNo channels selected for artifact_manual!\nSelect at least one channel in cfg.artfctdef.manual.channel'));
end;
channelindx=match_str(hdr.label, cfg.artfctdef.manual.channel);
ntrial=size(cfg.trl, 1);
nch=length(channelindx);
if(nch<1)
	error(sprintf('\nNo channels selected for artifact_manual!\nSelect at least one channel in cfg.artfctdef.manual.channel'));
elseif(nch>cfg.artfctdef.manual.maxnumberofchannels)
	error(sprintf('\nMore than %i channels selected in cfg.artfctdef.manual.channel',cfg.artfctdef.manual.maxnumberofchannels));
end

show=read_fcdc_data(cfg.datafile, hdr, 1, hdr.nTrials*hdr.nSamples, channelindx, iscontinuous);
show=show';

N=length(show);
S=floor(N./hdr.Fs); % In seconds
x=1:N;
x=x./hdr.Fs; % In seconds.

fprintf(' done.\n');
fprintf('Processing raw data for artifact_manual...');

% these elements are stored inside the figure so that the callback routines can modify them
dat.RejMarkList=zeros(length(cfg.trl),1);
dat.ResRejMarkList=zeros(length(cfg.trl),1);
dat.RejCount=0;
dat.ResRejCount=0;
dat.stop=0;
dat.trln=1;
dat.numtrl=length(cfg.trl);
dat.FFT=strcmp(cfg.artfctdef.manual.fft,'yes')|strcmp(cfg.artfctdef.manual.fft,'Yes')...
  |strcmp(cfg.artfctdef.manual.fft,'on')|strcmp(cfg.artfctdef.manual.fft,'On');
if(strcmp(cfg.artfctdef.manual.padding,'no')|strcmp(cfg.artfctdef.manual.padding,'No')...
    |strcmp(cfg.artfctdef.manual.padding,'Off')|strcmp(cfg.artfctdef.manual.padding,'off'))
  cfg.artfctdef.manual.padding=[];
end;
dat.IfRelOn=strcmp(cfg.artfctdef.manual.timeaxrelative,'yes')|strcmp(cfg.artfctdef.manual.timeaxrelative,'Yes')...
  |strcmp(cfg.artfctdef.manual.timeaxrelative,'on')|strcmp(cfg.artfctdef.manual.timeaxrelative,'On');

ival=cell(dat.numtrl,1);
dataX=cell(dat.numtrl,1);
dataY=cell(dat.numtrl,1);
dataFx=cell(dat.numtrl,1);
dataFy=cell(dat.numtrl,1);
dataXLim=cell(dat.numtrl,1);
dataYLim=cell(dat.numtrl,1);

% first we calculate everything and put it into varables.
for(i=1:dat.numtrl)
  begpadding = round(cfg.artfctdef.manual.pretrialtime.*hdr.Fs);
  endpadding = round(cfg.artfctdef.manual.posttrialtime.*hdr.Fs);
  ival{i}=(cfg.trl(i,1)-begpadding):(cfg.trl(i,2)+endpadding);
  dataY{i}=show(ival{i},:);

  % don't remove the padding
  [dataY{i}, dumlab, dumtime, dumcfg] = preproc(dataY{i}', cfg.artfctdef.manual.channel, hdr.Fs, cfg.artfctdef.manual, cfg.trl(i,3)-begpadding, 0, 0);
  dataY{i} = dataY{i}';

  if(~rem(i,fix(dat.numtrl/10)))fprintf('.');end;
end;
clear show; % Now, we don't need the complete dataset anymore...
for(i=1:dat.numtrl) % we go on with the calculations ....
  dataX{i}=x(ival{i})-dat.IfRelOn.*(x(ival{i}(1))+cfg.artfctdef.manual.pretrialtime);

  dataYLim{i}=[];
  for(j=1:nch)
    mini=min(dataY{i}(:,j));
    maxi=max(dataY{i}(:,j));

    dataYLim{i}=[dataYLim{i};[mini-.03.*(maxi-mini),maxi+.03.*(maxi-mini)]];
    dataXLim{i}=[min(dataX{i}),max(dataX{i})];
  end;

  if(dat.FFT) % if FFT is on, calculate the FFT of the trials too...
    tmp=[];
    for(k=1:nch) tmp=[tmp,dataY{i}(:,k)-mean(dataY{i}(:,k),1)];end;
    dataFy{i}=abs(fft(tmp,cfg.artfctdef.manual.padding,1));
    tmp=floor((max(ival{i})-min(ival{i}))/2);
    dataFx{i}=(1:tmp)./tmp.*hdr.Fs./2;
    dataFy{i}=dataFy{i}(1:length(dataFx{i}),:);
  end;
  if(~rem(i,fix(dat.numtrl/10)))fprintf('.');end;
end;

fprintf(' done.\n');

thisfig('fig_trace');% See the help of thisfig
set(fig_trace,'name','Traces');
if(dat.FFT);
  thisfig('fig_spec');
  set(fig_spec,'name','Spectra of traces');
end;
gt=[5 40 5 25 0 50 5 60 5 60 5 60 5 60 5 60 5 25 0 25];
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:1)),5,gt(2),18],'String','Done','Callback',@stop);
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:3)),5,gt(4),18],'String','<','Callback',@prevtrial);
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:5)),5,gt(6),18],'String','>','Callback',@nexttrial);
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:7)),5,gt(8),18],'String','Reject >','Callback',@reject_next);
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:9)),5,gt(10),18],'String','Reject','Callback',@reject);
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:11)),5,gt(12),18],'String','Accept','Callback',@undo);
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:13)),5,gt(14),18],'String','Accept All','Callback',@undoAll);
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:15)),5,gt(16),18],'String','Redo All','Callback',@redoAll);
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:17)),5,gt(18),18],'String','<<','Callback',@pptrial);
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:19)),5,gt(20),18],'String','>>','Callback',@nntrial);
set(fig_trace, 'KeyPressFcn',@keypress);
% set(get(fig_trace, 'children'),'KeyPressFcn',@keypress);
show_help;

thisfig('fig_trace'); % See the help of thisfig.

while(ishandle(fig_trace))

  if(dat.FFT) % Plots FFT when requested
    if(exist('HF')) % this is for remembering the zoom of the FFT-window
      if(ishandle(HF)) % while zapping trough the data.
        saveXLim=get(HF,'XLim');
        saveYLim=get(HF,'YLim');
      end;
    end;
    thisfig('fig_spec');
    plot(dataFx{dat.trln},dataFy{dat.trln});
    HF=get(fig_spec,'CurrentAxes');
    if(exist('saveXLim'))set(HF,'XLim',saveXLim);end;
    if(exist('saveYLim'))set(HF,'YLim',saveYLim);end;
    xlabel('Frequency [Hz]');
    ylabel('Power [A.U.]');
    thisfig('fig_trace');
  end;

  for(j=1:nch) % This loop fils the traces figure with the traces.
    H=subplot(nch,1,j);
    plot(dataX{dat.trln},dataY{dat.trln}(:,j),'b-');
    line([dataXLim{1}(1)+cfg.artfctdef.manual.pretrialtime dataXLim{1}(1)+cfg.artfctdef.manual.pretrialtime], dataYLim{dat.trln}(j,:), 'color', 'g');
    line([dataXLim{1}(2)-cfg.artfctdef.manual.posttrialtime dataXLim{1}(2)-cfg.artfctdef.manual.posttrialtime], dataYLim{dat.trln}(j,:), 'color', 'r');
    set(H,'YLim',dataYLim{dat.trln}(j,:) + [-eps +eps]);
    set(H,'XLim',dataXLim{dat.trln});
    xlabel('Time [s]','position',...
      [dataXLim{dat.trln}(2),dataYLim{dat.trln}(j,1)-(dataYLim{dat.trln}(j,2)-dataYLim{dat.trln}(j,1)).*.05]...
      ,'HorizontalAlignment','left');
    ylabel(sprintf('%s', cfg.artfctdef.manual.channel{j}));
  end;

  if(dat.RejMarkList(dat.trln)==0)
    tiet=sprintf('trial %i: ACCEPT',dat.trln);
  else
    tiet=sprintf('trial %i: REJECT',dat.trln);
  end;
  H=subplot(nch,1,1);
  title(tiet);

  guidata(fig_trace,dat);
  uiwait;
  if(ishandle(fig_trace)) dat=guidata(fig_trace);end;

  if(dat.stop)break;end;
end

if(~ishandle(fig_trace))
  error('Figure closed unexpectedly, no manual artifacts marked.');
else
  close(fig_trace)
end;

if(dat.FFT)
  if(ishandle(fig_spec))close(fig_spec);end;
end;

fprintf('Rejected %i trials.\n',dat.RejCount);
fprintf('Exported %i trials.\n',length(cfg.trl)-dat.RejCount);

artifact=cfg.trl(find((dat.RejMarkList)),[1,2]);

% remember the details that were used here
cfg.artfctdef.manual.artifact = artifact;

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: artifact_manual.m,v 1.14 2006/11/29 09:06:36 roboos Exp $';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here the SUBFUNCTIONS start taht implement the gui callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=nexttrial(h, eventdata, handles, varargin)
dat=guidata(h);
if(dat.trln<dat.numtrl)
  dat.trln=dat.trln + 1;
else
  dat.trln=1;
end;
guidata(h,dat);
uiresume;

function varargout=prevtrial(h, eventdata, handles, varargin)
dat=guidata(h);
if(dat.trln > 1),
  dat.trln=dat.trln - 1;
else
  dat.trln=dat.numtrl;
end;
guidata(h,dat);
uiresume;

function varargout=stop(h, eventdata, handles, varargin)
dat=guidata(h);
dat.stop=1;
guidata(h,dat);
uiresume;

function varargout=reject(h, eventdata, handles, varargin)
dat=guidata(h);
if(dat.RejMarkList(dat.trln)==0)
  dat.RejCount=dat.RejCount+1;
  dat.RejMarkList(dat.trln)=1;
  if(~dat.ResRejMarkList(dat.trln))
    dat.ResRejCount=dat.ResRejCount+1;
    dat.ResRejMarkList(dat.trln)=1;
    fprintf('Rejected trial %i!\n',dat.trln);
  else
    fprintf('Rejected trial %i again!\n',dat.trln);
  end;
else
  fprintf('Trial %i is already marked for rejection.\n',dat.trln);
end;
guidata(h,dat);
uiresume;


function varargout=reject_next(h, eventdata, handles, varargin)
dat=guidata(h);
if(dat.RejMarkList(dat.trln)==0)
  dat.RejCount=dat.RejCount+1;
  dat.RejMarkList(dat.trln)=1;
  if(~dat.ResRejMarkList(dat.trln))
    dat.ResRejCount=dat.ResRejCount+1;
    dat.ResRejMarkList(dat.trln)=1;
    fprintf('Rejected trial %i!\n',dat.trln);
  else
    fprintf('Rejected trial %i again!\n',dat.trln);
  end;
else
  fprintf('Trial %i is already marked for rejection.\n',dat.trln);
end;
if(dat.trln < dat.numtrl)
  dat.trln=dat.trln + 1;
else
  dat.trln=1;
end;
guidata(h,dat);
uiresume;

function varargout=undo(h, eventdata, handles, varargin)
dat=guidata(h);
if(dat.RejMarkList(dat.trln)>0)
  dat.RejCount=dat.RejCount-1;
  dat.RejMarkList(dat.trln)=0;
  fprintf('Rejected trial %i recovered.\n',dat.trln);
end;
guidata(h,dat);
uiresume;

function varargout=undoAll(h, eventdata, handles, varargin)
dat=guidata(h);
dat.RejCount=0;
dat.RejMarkList=zeros(dat.numtrl,1);
fprintf('All rejected trials recovered.\n');
guidata(h,dat);
uiresume;

function varargout=redoAll(h, eventdata, handles, varargin)
dat=guidata(h);
dat.RejMarkList=dat.ResRejMarkList;
dat.RejCount=dat.ResRejCount;
fprintf('All recovered trials rejected again.\n');
guidata(h,dat);
uiresume;

function varargout=nntrial(h, eventdata, handles, varargin)
dat=guidata(h);
if(dat.trln<(dat.numtrl-10))
  dat.trln=dat.trln + 10;
else
  dat.trln=1;
end;
guidata(h,dat);
uiresume;


function varargout=pptrial(h, eventdata, handles, varargin)
dat=guidata(h);
if(dat.trln > 10),
  dat.trln=dat.trln - 10;
else
  dat.trln=dat.numtrl;
end;
guidata(h,dat);
uiresume;

function keypress(h, eventdata, handles, varargin)
dat=guidata(gcbf);
key=get(gcbf, 'CurrentCharacter');
if  key
  switch key
    case 28     % arrow left
      if(dat.trln > 1),
        dat.trln=dat.trln - 1;
      else
        dat.trln=dat.numtrl;
      end;
    case 29     % arrow right
      if(dat.trln<dat.numtrl)
        dat.trln=dat.trln + 1;
      else
        dat.trln=1;
      end;
      %    case 30     % arrow up
      %      dat.zscale=data.zscale*1.5;
      %    case 31     % arrow down
      %      dat.zscale=data.zscale/1.5;
      %    case 'a'
      %      dat.reject(data.trlop)=0;
      %    case 'r'
      %      dat.reject(data.trlop)=1;
    case {'a', 'A'}
      if(dat.RejMarkList(dat.trln)==1)
        dat.RejCount=dat.RejCount-1;
        dat.RejMarkList(dat.trln)=0;
        fprintf('Accepted trial %i\n',dat.trln);
      end 
      if(dat.trln < dat.numtrl)
        dat.trln=dat.trln + 1;
      else
        dat.trln=1;
      end;
    case {'0', 'r', 'R'}
      if(dat.RejMarkList(dat.trln)==0)
        dat.RejCount=dat.RejCount+1;
        dat.RejMarkList(dat.trln)=1;
        if(~dat.ResRejMarkList(dat.trln))
          dat.ResRejCount=dat.ResRejCount+1;
          dat.ResRejMarkList(dat.trln)=1;
          fprintf('Rejected trial %i! \n',dat.trln);
        else
          fprintf('Rejected trial %i again! \n',dat.trln);
        end;
      else
        fprintf('Trial %i is already marked for rejection.\n',dat.trln);
      end;
      if(dat.trln < dat.numtrl)
        dat.trln=dat.trln + 1;
      else
        dat.trln=1;
      end;

    otherwise
      show_help;
  end
end
set(gcbf, 'CurrentCharacter', '-');
guidata(gcbf, dat);
uiresume(h);

function show_help;
fprintf('You can use the buttons in the window to navigate through the trials,   \n');
fprintf('reject trials or undo and redo your rejection selection. If you want to \n');
fprintf('use the keyboard, you first have to click in the figure every time you  \n');
fprintf('pressed a button with the mouse.                                        \n');
fprintf('                                                                        \n');
fprintf('Use the left and right arrow to walk through the trials                 \n');
fprintf('Press "a" to accept, "r" to reject.                                     \n');

function thisfig(NameOfHandle);
assignin('caller','This_Name_Will_Never_be_Used_gk',NameOfHandle);
if(~evalin('caller','exist(This_Name_Will_Never_be_Used_gk)')|~ishandle(evalin('caller',NameOfHandle)))
  H=figure;
  assignin('caller',NameOfHandle,H);
else
  figure(evalin('caller',NameOfHandle));
end;
evalin('caller','clear This_Name_Will_Never_be_Used_gk;');
