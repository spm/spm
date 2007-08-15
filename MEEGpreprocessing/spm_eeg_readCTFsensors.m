function [SensLoc, SensOr, CoiLoc, Label] = spm_eeg_readCTFsensors

ctffolder = spm_select(1, 'dir', 'Select the CTF folder');

pre_data = ctf_read(ctffolder, 'meg', [0 0.01]);
% --- Save coil/sensor positions and orientations for source reconstruction (in mm) ---

% - channel locations and orientations
SensLoc = [];
SensOr  = [];
for i = 1:length(pre_data.sensor.location);
    if any(pre_data.sensor.location(:,i)) & pre_data.sensor.label{i}(1) == 'M'
        SensLoc = [SensLoc; pre_data.sensor.location(:,i)'];
        SensOr  = [SensOr ; pre_data.sensor.orientation(:,i)'];
    end
end
SensLoc = 10*SensLoc; % convertion from cm to mm
if length(SensLoc) > 275
    warning(sprintf('Found more than 275 channels!\n'));
end

% - coil locations (in this order - NZ:nazion , LE: left ear , RE: right ear)
CurrentDir = pwd;
cd(pre_data.folder);
hc_files = dir('*.hc');
if isempty(hc_files)
    warning(sprintf('Impossible to find head coil file\n'));
elseif length(hc_files) > 1
    hc_file = spm_select(1, '\.hc$', 'Select head coil file');
else
    hc_file = fullfile(pre_data.folder,hc_files.name);
end
clear hc_files
for coils=1:3
    fid = fopen(hc_file,'r');
    testlines = fgetl(fid);
    t=0;
    while t==0
        if coils==1 & strmatch('measured nasion coil position relative to head (cm):',testlines);
            t=1;
            for i = 1:3 % Nazion coordinates
                UsedLine    = fgetl(fid);
                UsedLine    = fliplr(deblank(fliplr(UsedLine)));
                [A,COUNT,ERRMSG,NEXTINDEX] = sscanf(UsedLine,'%c = %f');
                if ~isempty(ERRMSG) | (COUNT ~= 2)
                    warning(sprintf('Unable to read head coil file\n'));
                else
                    NZ(i) = A(2);
                end
            end
        end
        if coils==2 & strmatch('measured left ear coil position relative to head (cm):',testlines);
            t=1;
            for i = 1:3 % Nazion coordinates
                UsedLine    = fgetl(fid);
                UsedLine    = fliplr(deblank(fliplr(UsedLine)));
                [A,COUNT,ERRMSG,NEXTINDEX] = sscanf(UsedLine,'%c = %f');
                if ~isempty(ERRMSG) | (COUNT ~= 2)
                    warning(sprintf('Unable to read head coil file\n'));
                else
                    LE(i) = A(2);
                end
            end
        end
        if coils==3 & strmatch('measured right ear coil position relative to head (cm):',testlines);
            t=1;
            for i = 1:3 % Nazion coordinates
                UsedLine    = fgetl(fid);
                UsedLine    = fliplr(deblank(fliplr(UsedLine)));
                [A,COUNT,ERRMSG,NEXTINDEX] = sscanf(UsedLine,'%c = %f');
                if ~isempty(ERRMSG) | (COUNT ~= 2)
                    warning(sprintf('Unable to read head coil file\n'));
                else
                    RE(i) = A(2);
                end
            end
        end
        testlines = fgetl(fid);
    end
    fclose(fid);
end


CoiLoc = 10*[NZ ; LE ; RE]; % convertion from cm to mm
cd(CurrentDir);


Label=pre_data.sensor.label;