function artifact_viewer(cfg,artfctdef,rejectall,z_tdata,artfctchn,flag)

% artifact_viewer provides a GUI for browsing through the data while applying automatic artifact
% detection. it presents the data on trial by trial, and ad lib. one can jump through the data.
% the figure shows two subplot, the lower one containing the transformed data, on which the 
% thresholding has been applied. the upper window contains the data of the channel (which is a 
% member of the subset, defined in xxx.artfctdef.sgn), contributing most to the most deviant 
% transformed data-point in the respective trial.
%
% Copyright (c) 2004,  F.C. Donders Centre
%

dat.trln= 1;
dat.cfg = cfg;
dat.rej = rejectall;
dat.zdat= z_tdata;
dat.art = artfctdef;
dat.chn = artfctchn;
dat.hdr = read_fcdc_header(cfg.headerfile);
dat.trl = artfctdef.trl;
dat.numtrl = size(dat.trl,1);
dat.nsmp   = dat.trl(:,2) - dat.trl(:,1) + 1;
dat.stop   = 0;
dat.afv    = cell(1,size(dat.trl,1));
for j = 1:size(dat.trl,1)
   dat.afv{j}    = rejectall(1,dat.trl(j,1):dat.trl(j,2));
   dat.afv{j}     = diff([0 dat.afv{j} 0]);
   if isempty(find(dat.afv{j}))
     dat.trialok(j) = 1;
   else
     dat.trialok(j) = 0;
   end
end

h = figure;
guidata(h,dat);
uicontrol(gcf,'units','pixels','position',[5 5 40 18],'String','stop','Callback',@stop);
uicontrol(gcf,'units','pixels','position',[50 5 25 18],'String','<','Callback',@prevtrial);
uicontrol(gcf,'units','pixels','position',[75 5 25 18],'String','>','Callback',@nexttrial);
uicontrol(gcf,'units','pixels','position',[105 5 25 18],'String','<<','Callback',@prev10trial);
uicontrol(gcf,'units','pixels','position',[130 5 25 18],'String','>>','Callback',@next10trial);
uicontrol(gcf,'units','pixels','position',[160 5 50 18],'String','<artfct','Callback',@prevartfct);
uicontrol(gcf,'units','pixels','position',[210 5 50 18],'String','artfct>','Callback',@nextartfct);

while ishandle(h),
    dat = guidata(h);
    if dat.stop == 0,
      read_and_plot(h);
      uiwait;  
    else
      break
    end         
end
if ishandle(h)
  close(h);
end

%------------
%subfunctions
%------------

function read_and_plot(h)
   dat = guidata(h);
  
   rej       = [];
   trln      = dat.trln;
   cfg       = dat.cfg;
   rejectall = dat.rej;
   z_tdata   = dat.zdat;
   artfctdef = dat.art;
   hdr       = dat.hdr;
   trl       = dat.trl;
   nsmp      = dat.nsmp;
   chn       = dat.chn(trln);
   channel   = dat.art.sgn(chn,:);
   sgnind    = match_str(dat.hdr.label,channel); 
   afv       = dat.afv{trln};

   if ~isfield(cfg, 'datatype') || ~strcmp(cfg.datatype, 'continuous')
     % datatype is unknown or not continuous, perform epoch boundary check
     iscontinuous = 0;
   else
     % do not perform epoch boundary check, usefull for pseudo-continuous data
     iscontinuous = strcmp(cfg.datatype, 'continuous');
   end 

   if ~isfield(artfctdef, 'fltpadding')
     if iscontinuous
       if ~isfield(artfctdef, 'bpfilttype')
         fltpadding = 0;
       elseif strcmp(artfctdef.bpfilttype, 'but')
         fltpadding = 20*artfctdef.bpfiltord;
       elseif strcmp(artfctdef.bpfilttype, 'fir')
         fltpadding = artfctdef.bpfiltord;
       else
         warning('unknown filter type, cannot determine filter padding');
         fltpadding = 0;
       end
     else
       warning('default is not to apply filter padding on trial based data');
       fltpadding = 0;
     end
   else
     fltpadding = round(artfctdef.fltpadding*hdr.Fs);
   end
  
   if (strcmp(flag,'EOGflt')|strcmp(flag,'jump'))
       data    = read_fcdc_data(cfg.datafile,hdr,trl(trln,1)-fltpadding, ...
           trl(trln,2)+fltpadding,sgnind,iscontinuous);
       data    = preproc(data, channel, hdr.Fs, artfctdef, [], fltpadding, fltpadding);
       data    = blc(data);
   else
       data    = read_fcdc_data(cfg.datafile,hdr,trl(trln,1),trl(trln,2),sgnind,iscontinuous);
       data    = blc(data);
   end;
   if ~isempty(find(afv==1))
       rej(:,1)= find(afv==1)';
       rej(:,2)= (find(afv==-1)-1)';
   end
   s=sprintf('trial: %d',trln);
   
   %plot data of most aberrant channel in upper subplot
   gcf;subplot(211);plot(0,0);hold on;title(s);
   if ~isempty(rej)
       for artfct=1:size(rej,1)
           xpos   = rej(artfct,1)-round(artfctdef.pretim*hdr.Fs);
           ypos   = min(min(data(:)));
           width  = round((artfctdef.psttim+artfctdef.pretim)*hdr.Fs) + ...
               rej(artfct,2)-rej(artfct,1) + 1;
           % height = diff(minmax(data(1,:)));
           height = range(data(1,:));
           rectangle('Position',[xpos ypos width height],'FaceColor','r');
           hold on;
       end;
   end;
   plot(1:nsmp(trln),data,'b');
   lowlim = min(data(1,:)) - 0.1*(abs(min(data(1,:))));
   uplim = max(data(1,:))+ 0.1*(abs(max(data(1,:))));
   axis([0 nsmp(trln) lowlim uplim]);
   hold off;
   xlabel('samples');
   ylabel('uV or Tesla');
   
   %plot the transformed data, on which the thresholding is applied, in the lower subplot
   gcf;subplot(212);plot(1:nsmp(trln),z_tdata(:,trl(trln,1):trl(trln,2)));hold on;
   plot([1 nsmp(trln)],[artfctdef.cutoff artfctdef.cutoff],'r:');
   lowlim = min(z_tdata(1,trl(trln,1):trl(trln,2))) - ...
       0.1*(abs(min(z_tdata(1,trl(trln,1):trl(trln,2)))));
   uplim = max([max(z_tdata(1,trl(trln,1):trl(trln,2)))+ ...
       0.1*(abs(max(z_tdata(1,trl(trln,1):trl(trln,2))))) artfctdef.cutoff*1.05]);
   axis([0 nsmp(trln) lowlim uplim]);
   hold off;
   xlabel('samples');
   ylabel('zscore');
   guidata(h,dat);
       
function varargout = nexttrial(h, eventdata, handles, varargin)
   dat = guidata(h);
   if dat.trln < dat.numtrl,
      dat.trln = dat.trln + 1;
   end;
   guidata(h,dat);
   uiresume;

function varargout = next10trial(h, eventdata, handles, varargin)
   dat = guidata(h);
   if dat.trln < dat.numtrl - 10,
      dat.trln = dat.trln + 10;
   else dat.trln = dat.numtrl;
   end;
   guidata(h,dat);
   uiresume;
   
function varargout = prevtrial(h, eventdata, handles, varargin)
   dat = guidata(h);
   if dat.trln > 1,
      dat.trln = dat.trln - 1;
   else dat.trln = 1;
   end;
   guidata(h,dat);
   uiresume;
   
function varargout = prev10trial(h, eventdata, handles, varargin)
   dat = guidata(h);
   if dat.trln > 10,
      dat.trln = dat.trln - 10;
   else dat.trln = 1;
   end;
   guidata(h,dat);
   uiresume;

function varargout = nextartfct(h, eventdata, handles, varargin)
   dat = guidata(h);
   artfctindx = find(dat.trialok == 0);
   sel = find(artfctindx > dat.trln);
   if ~isempty(sel)
     dat.trln = artfctindx(sel(1));
   else
     dat.trln = dat.trln;
   end
   guidata(h,dat);
   uiresume;

function varargout = prevartfct(h, eventdata, handles, varargin)
   dat = guidata(h);
   artfctindx = find(dat.trialok == 0);
   sel = find(artfctindx < dat.trln);
   if ~isempty(sel)
     dat.trln = artfctindx(sel(end));
   else
     dat.trln = dat.trln;
   end
   guidata(h,dat);
   uiresume;

function varargout = stop(h, eventdata, handles, varargin)
   dat = guidata(h);
   dat.stop = 1;
   guidata(h,dat);
   uiresume;

