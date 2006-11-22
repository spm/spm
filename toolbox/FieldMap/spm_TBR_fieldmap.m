function spm_TBR_fieldmap
% Convert Raw Siemens DICOM files from fieldmap sequence into something 
% that SPM can use
% FORMAT spm_TBR_fieldmap
% Converted files are written to the current directory
%_______________________________________________________________________
% %W% %E%
%17.10.2003: change of smooth parameter from 2 to 16 in trajrecon
%10.12.2003: smooth parameter can be entered from dialog box
%11.02.2004: update version for converting fieldmap data
fprintf(1,'TBR 11.02.2004 Version for converting field map data  \n');

spm_defaults;
global hdr

P                    = spm_get(Inf,'*','Select files');

%SliceOrder          = spm_input('Slice Order', 2, 'm',...
%                      'Ascending Slices|Descending Slices', str2mat('ascending','descending'));
SliceOrder           = 'ascending';

% Select trajectory file
[d,f,s]    = fileparts(which(mfilename));
t          = load(spm_get(1,'traj_*_????-??-??_??-??.mat','Select trajectory',d));

% Select smoothing parameter
myline1='Enter the smoothing parameter in the following dialog box';
myline2='The default value is 4.0 pixels';
myline3='Choose a lower value if ghosting occurs in the images';
myline4='Choose a higher value if there are line artefacts in the images';
uiwait(msgbox({myline1,myline2,myline3,myline4},'Smoothing Info','modal'));

myprompt={'Enter the smoothing parameter (in pixels):'};
mydef={'4.0'};
mydlgTitle=['Smoothing parameter'];
mylineNo=1;
answer=inputdlg(myprompt,mydlgTitle,mylineNo,mydef);
myzelle=answer(1);
smooth_par=str2num(myzelle{1});

hdr        = spm_dicom_headers(P);
[raw,guff] = select_fid(hdr);

% Find fieldmap data and extract it from the raw TBR data to be converted
[raw,pm_raw] = select_series(raw,cat(2,'ralf_flash_fieldmap',' ','chloe_epi_tbr',' ','ralf_epi_fieldmap_ff','ralf_epi_fieldmap_ff_shim'));

%spm_progress_bar('Init',length(hdr),['Writing Images'], 'Files written');
%for i=1:length(raw),
%	convert_fid(raw{i},t.M,SliceOrder,smooth_par);
%	spm_progress_bar('Set',i);
%end;
%spm_progress_bar('Clear');

spm_progress_bar('Init',length(pm_raw),['Processing Fieldmap'], 'Files written');
if ~isempty(pm_raw)
   % Process fieldmap data if acquired
   % Convert double echo fid
   convert_pm(pm_raw,t.M,SliceOrder,smooth_par);
end
spm_progress_bar('Clear');

return;
%_______________________________________________________________________

%_______________________________________________________________________
function [fid,other] = select_fid(hdr)
fid      = {};
other    = {};
for i=1:length(hdr),
	if ~checkfields(hdr{i},'ImageType','StartOfCSAData') |...
		~strcmp(deblank(hdr{i}.ImageType),'ORIGINAL\PRIMARY\RAW'),
		other = {other{:},hdr{i}};
	else,
		fid = {fid{:},hdr{i}};
	end;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function ok = checkfields(hdr,varargin)
ok = 1;
for i=1:(nargin-1),
	if ~isfield(hdr,varargin{i}),
		ok = 0;
		break;
	end;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function clean = strip_unwanted(dirty)
msk = find((dirty>='a'&dirty<='z') | (dirty>='A'&dirty<='Z') |...
           (dirty>='0'&dirty<='9') | dirty=='_');
clean = dirty(msk);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function vol=trajrecon(meas,M,smooth_par)
% Oliver Joseph's TBR reconstruction
%17.10.2003: change of smooth parameter from 2 to 16
%11.11.2003: smoothing parameter has become input argument

% Regridded FT!!!!
% fid ( LINES, POINTS, SLICES )!!!
%-------------------------------------------------------------------
f1d=zeros(64,size(meas,2),size(meas,3));
for line=1:size(meas,2)
 	f1d(:,line,:)=M(:,:,line)*squeeze(meas(:,line,:));
end;

% Weighted phase reference smoothing
%-------------------------------------------------------------------
x   = fftshift(-32:31)';
y   = exp(-(x.^2)*(smooth_par/64)^2);
dd  = f1d(:,65,:).*conj(f1d(:,66,:));
df  = fft(dd,[],1);
df  = df.*repmat(y,[1 1 size(meas,3)]);
dff = ifft(df,[],1);
d2  = exp(sqrt(-1)*angle(dff));

% Slice Based Phase correction %
%-------------------------------------------------------------------
f1d(:,2:2:64,:) = f1d(:,2:2:64,:).*repmat(d2,[1,32,1]);

% 2nd FT
% Need loop to prevent odd-number-of-slices fftshift problem
%-------------------------------------------------------------------
% Chloe's changes 27/02/03 
% Need an additional conj to extract correct phase direction 
% This doesn't effect the magnitude data
%for slice=1:size(f1d,3),
% 	vol(:,:,slice)=fftshift(fft(fftshift(conj(f1d(:,1:64,slice))),[],2));
%end

for slice=1:size(f1d,3),
   vol(:,:,slice)=conj(fftshift(fft(fftshift(conj(f1d(:,1:64,slice))),[],2)));
end

% Don't calculate the magnitude here because we need complex data in case
% image is field map data
% vol = abs(vol);
return;

%_______________________________________________________________________

%_______________________________________________________________________
function convert_pm(hdr,M,SliceOrder,smooth_par)
    
% Do the right thing with each type of field map

[other_hdr, fmap_hdr] = select_series(hdr,cat(2,'ralf_epi_fieldmap_ff','ralf_flash_fieldmap','ralf_epi_fieldmap_ff_shim'));

if ~isempty(fmap_hdr)
   % Sort fieldmap data into short/long echo pairs ready for conversion if require
   Fhdr=sort_fieldmap_data(fmap_hdr);
   nfmap=length(Fhdr)/2;
   TE={'short','long'};    
   % Processing each fieldmap requires converting two dicom files

   for i=1:2:length(Fhdr) % Process each fieldmap as pair of images

      % Need 4 cells to hold real and imaginary pair
      npm=1;
      V=cell(1,4);

      for j=1:2 % Process each pair of images
         index=(i-1)+j;
         hdr=Fhdr{index};

         % Sort out output filename
         %-------------------------------------------------------------------
%        fname = sprintf('f%s-%.4d.img', strip_unwanted(hdr.PatientID),hdr.SeriesNumber);
         fname = sprintf('f%s-%.4d-%.6d.img', strip_unwanted(hdr.PatientID),...
		hdr.SeriesNumber, hdr.AcquisitionNumber);

         % Order of files for the field map routines is
         % short_real, short_imag, short_mag 
         % long_real, long_imag, long_mag (are mags this really necessary?)

         TEname = sprintf('%s',TE{j});
         fname = [deblank(spm_str_manip(fname,'rt')) '-' deblank(TEname)];
        
         [V{npm},V{npm+1}] = complex2spm(hdr,M,SliceOrder,fname,smooth_par);
         npm=npm+2;  
      end      
   end
end

if ~isempty(other_hdr)
   for i=1:length(other_hdr) 
   
      hdr=other_hdr{i};

      % Sort out output filename
      %-------------------------------------------------------------------
%     fname = sprintf('f%s-%.4d.img', strip_unwanted(hdr.PatientID),hdr.SeriesNumber);
      fname = sprintf('f%s-%.4d-%.6d.img', strip_unwanted(hdr.PatientID),...
	      hdr.SeriesNumber, hdr.AcquisitionNumber);

      % Order of files for the field map 
      TEname = sprintf('TE%s',num2str(hdr.EchoTime));
      fname = [deblank(spm_str_manip(fname,'rt')) '-' deblank(TEname)];
      complex2spm(hdr,M,SliceOrder,fname,smooth_par); 
   end
end

return;

%_______________________________________________________________________

%_______________________________________________________________________
function varargout = complex2spm(hdr,M,SliceOrder,fname,smooth_par)

persistent SpacingBetweenSlices

% Dimensions and data
%-------------------------------------------------------------------
if hdr.Columns*hdr.Rows*2 ~= hdr.SizeOfCSAData,
   warning(['Cant convert "' fname '".']);
   return;
end;

nphase = 66; % <- may change
dim    = [hdr.Columns/4 nphase hdr.Rows/nphase];
fd     = fopen(hdr.Filename,'r','ieee-le');

if fd==-1,
   warning(['Cant open "' hdr.Filename '" to create "' fname '".']);
   return;
end;

fseek(fd,hdr.StartOfCSAData,-1);
cdata  = fread(fd,hdr.SizeOfCSAData/4,'float');
fclose(fd);
cdata  = reshape(cdata(1:2:end)+sqrt(-1)*cdata(2:2:end),dim);

if ~isempty(findstr('flash_fieldmap',hdr.SeriesDescription))
   volume = flashrecon(cdata);
else
   volume = trajrecon(cdata,M,smooth_par);
end
dim    = [size(volume) spm_type('int16')];

% Orientation information as for TBR above
%-------------------------------------------------------------------
analyze_to_dicom = [diag([1 -1 1]) [0 (dim(2)+1) 0]'; 0 0 0 1]*[eye(4,3) [-1 -1 -1 1]'];
if isfield(hdr,'SpacingBetweenSlices'),
   SpacingBetweenSlices = hdr.SpacingBetweenSlices;
else,
   if isempty(SpacingBetweenSlices),
      SpacingBetweenSlices = spm_input('Spacing Between Slices (mm)',1,'r','3.0',[1 1]);
   end;
end;

vox    = [hdr.PixelSpacing SpacingBetweenSlices];
pos    = hdr.ImagePositionPatient';
orient = reshape(hdr.ImageOrientationPatient,[3 2]);
orient(:,3) = null(orient');
if det(orient)<0, orient(:,3) = -orient(:,3); end;

% The image position vector is not correct. In dicom this vector points to
% the upper left corner of the image. Perhaps it is unlucky that this is
% calculated in the syngo software from the vector pointing to the center of
% the slice (keep in mind: upper left slice) with the enlarged FoV.

dicom_to_patient = [orient*diag(vox) pos ; 0 0 0 1];
truepos          = dicom_to_patient *[([hdr.Columns hdr.Rows]-dim(1:2))/2 0 1]';
dicom_to_patient = [orient*diag(vox) truepos(1:3) ; 0 0 0 1];
patient_to_tal   = diag([-1 -1 1 1]);
mat              = patient_to_tal*dicom_to_patient*analyze_to_dicom;

% Re-arrange order of slices
switch SliceOrder,
   case {'descending'},
      order = 1:dim(3);
   case {'ascending'},
      order = dim(3):-1:1;
   otherwise,
      error('Unknown slice order');
   end;

volume(:,:,order) = volume;
mat    = mat*[eye(3) [0 0 -(order(1)-1)]'; 0 0 0 1];
      
% Really need to have SliceNormalVector - but it doesn't exist.
% Also need to do something else for interleaved volumes.
%-------------------------------------------------------------------
% SliceNormalVector = read_SliceNormalVector(hdr);
% if det([reshape(hdr.ImageOrientationPatient,[3 2]) SliceNormalVector(:)])<0;
% 	volume = volume(:,:,end:-1:1);
% 	mat    = mat*[eye(3) [0 0 -(dim(3)-1)]'; 0 0 0 1];
%  end;

% Possibly useful information
%-------------------------------------------------------------------
tim = datevec(hdr.AcquisitionTime/(24*60*60));
descrip = sprintf('%s %s TR=%gms/TE=%gms/FA=%gdeg %s %d:%d:%.5g TBR',...
	     hdr.MRAcquisitionType,...
	     deblank(hdr.ScanningSequence),...
	     hdr.RepetitionTime,hdr.EchoTime,hdr.FlipAngle,...
	     datestr(hdr.AcquisitionDate),tim(4),tim(5),tim(6));
     
% Construct name and write out real part of complex data
realname=sprintf('%s','real');
realfname=[deblank(spm_str_manip(fname,'rt')) '-' deblank(realname) '.img'];
V{1}=struct('fname',realfname,'dim',dim,'mat',mat,'descrip',descrip);
V{1}=spm_write_vol(V{1},real(volume));
        
% Construct name and write out imaginary part of complex data
imagname=sprintf('%s','imag');
imagfname=[deblank(spm_str_manip(fname,'rt')) '-' deblank(imagname) '.img'];
V{2}=struct('fname',imagfname,'dim',dim,'mat',mat,'descrip',descrip);
V{2}=spm_write_vol(V{2},imag(volume));
  
% Construct name and write out magnitude of complex data
magfname=[deblank(spm_str_manip(fname,'rt')) '.img'];
Vi=struct('fname',magfname,'dim',dim,'mat',mat,'descrip',descrip);
Vi=spm_write_vol(Vi,abs(volume));
  
varargout{1}=V{1};
varargout{2}=V{2};

return;

%_______________________________________________________________________
function F = sort_fieldmap_data(hdr)

% Want to extract the structures in the correct order to 
% create a real and imaginary and a magnitude image for
% two volumes (or more?).

nfiles=length(hdr);
if nfiles==1
   error('Need at least 2 dicom files for field map conversion');
   return
end

sessions=[];
for i=1:nfiles
        sessions(i)=hdr{i}.SeriesNumber;
end;

% Sort scans into sessions
fnum=1;
nfmap=1;
while fnum<=nfiles
   nsess=find(sessions==sessions(fnum));
   for i=1:length(nsess)
      sesshdr{i}=hdr{nsess(i)};
   end
   % Split field maps up into pairs
   n=2:2:length(nsess);
   if ~isempty(n)
%     These lines are commented out so that only the last fieldmaps 
%     images are written.
%
%     for i=1:length(n)
%        F{nfmap}=sesshdr{n(i)-1};
%        F{nfmap+1}=sesshdr{n(i)};
%        nfmap=nfmap+2;
%     end
      F{nfmap}=sesshdr{n(end)-1};
      F{nfmap+1}=sesshdr{n(end)};
      nfmap=nfmap+2;
   end
   fnum=fnum+length(nsess);
end

return;
%_______________________________________________________________________

%_______________________________________________________________________
function [fid,sfid] = select_series(hdr,seriesname)
fid      = {};
sfid    = {};
for i=1:length(hdr),
	if findstr(deblank(hdr{i}.SeriesDescription),seriesname),
		sfid = {sfid{:},hdr{i}};
	else,
		fid = {fid{:},hdr{i}};
	end;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function vol=flashrecon(cdata)
% Ralf Deichman's flash recon

[ncol nphase nslc]=size(cdata);
cdata=ifftshift(fft(fftshift(cdata),[],1));

% Double conj to reverse y direction for consistency with TBR data
% in y-direction, but without changing phase direction.
vol=conj(ifftshift(fft(fftshift(conj(cdata(:,1:nphase-2,:))),[],2)));
vol=vol(ncol/4+1:ncol*3/4,:,:);
dim=size(vol);

