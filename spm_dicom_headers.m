function hdr = spm_dicom_headers(P, essentials)
% Read header information from DICOM files
% FORMAT hdr = spm_dicom_headers(P [,essentials])
% P   - array of filenames
% essentials - if true, then only save the essential parts of the header
% hdr - cell array of headers, one element for each file.
%
% Contents of headers are approximately explained in:
% http://medical.nema.org/dicom/2001.html
%
% This code will not work for all cases of DICOM data, as DICOM is an
% extremely complicated "standard".
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dicom_headers.m 5250 2013-02-15 21:04:36Z john $

if nargin<2, essentials = false; end

dict = readdict;
j    = 0;
hdr  = {};
if size(P,1)>1, spm_progress_bar('Init',size(P,1),'Reading DICOM headers','Files complete'); end;
for i=1:size(P,1),
    tmp = readdicomfile(P(i,:),dict);
    if ~isempty(tmp),
        if essentials, tmp = spm_dicom_essentials(tmp); end
        j      = j + 1;
        hdr{j} = tmp;
    end;
    if size(P,1)>1, spm_progress_bar('Set',i); end;
end;
if size(P,1)>1, spm_progress_bar('Clear'); end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function ret = readdicomfile(P,dict)
ret = [];
P   = deblank(P);
fp  = fopen(P,'r','ieee-le');
if fp==-1, warning('spm:dicom','%s: Cant open file.',P); return; end;

fseek(fp,128,'bof');
dcm = char(fread(fp,4,'uint8')');
if ~strcmp(dcm,'DICM'),
    % Try truncated DICOM file fomat
    fseek(fp,0,'bof');
    tag.group   = fread(fp,1,'ushort');
    tag.element = fread(fp,1,'ushort');
    if isempty(tag.group) || isempty(tag.element),
        fclose(fp);
        warning('spm:dicom','%s: Truncated file.',P);
        return;
    end;
    %t          = dict.tags(tag.group+1,tag.element+1);
    if isempty(find(dict.group==tag.group & dict.element==tag.element,1)) && ~(tag.group==8 && tag.element==0),
        % entry not found in DICOM dict and not from a GE Twin+excite
        % that starts with with an 8/0 tag that I can't find any
        % documentation for.
        fclose(fp);
        warning('spm:dicom','%s: Not a DICOM file.', P);
        return;
    else
        fseek(fp,0,'bof');
    end;
end;
try
    ret = read_dicom(fp, 'il',dict);
    ret.Filename = fopen(fp);
catch problem
    fprintf('%s: Trouble reading DICOM file (%s), skipping.\n', fopen(fp), problem.message);
end
fclose(fp);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [ret,len] = read_dicom(fp, flg, dict,lim)
if nargin<4, lim=4294967295; end; % FFFFFFFF
len = 0;
ret = [];
while len<lim,
    tag = read_tag(fp,flg,dict);
    if isempty(tag), break; end

    %fprintf('(%.4X,%.4X) "%s" %d %d %s\n', tag.group, tag.element, tag.vr, tag.length, tag.le, tag.name);

    if tag.group==65534 && tag.element==57357, % FFFE,E00D ItemDelimitationItem
        break;
    end
 
    if tag.length>0,
        switch tag.name,
            case {'GroupLength'},
                % Ignore it
                fseek(fp,tag.length,'cof');
            case {'PixelData'},
                ret.StartOfPixelData = ftell(fp);
                ret.SizeOfPixelData  = tag.length;
                ret.VROfPixelData    = tag.vr;
                fseek(fp,tag.length,'cof');
            case {'CSAData'}, % raw data
                ret.StartOfCSAData = ftell(fp);
                ret.SizeOfCSAData = tag.length;
                fseek(fp,tag.length,'cof');
            case {'CSAImageHeaderInfo', 'CSASeriesHeaderInfo','Private_0029_1110','Private_0029_1120','Private_0029_1210','Private_0029_1220'},
                dat  = decode_csa(fp,tag.length);
                ret.(tag.name) = dat;
            case {'TransferSyntaxUID'},
                dat = char(fread(fp,tag.length,'uint8')');
                dat = deblank(dat);
                ret.(tag.name) = dat;
                switch dat,
                    case {'1.2.840.10008.1.2'},      % Implicit VR Little Endian
                        flg = 'il';
                    case {'1.2.840.10008.1.2.1'},    % Explicit VR Little Endian
                        flg = 'el';
                    case {'1.2.840.10008.1.2.1.99'}, % Deflated Explicit VR Little Endian
                        warning('spm:dicom','%s: Cant read Deflated Explicit VR Little Endian file.', fopen(fp));
                        flg = 'dl';
                        return;
                    case {'1.2.840.10008.1.2.2'},    % Explicit VR Big Endian
                        %warning('spm:dicom','%s: Cant read Explicit VR Big Endian file',fopen(fp));
                        flg = 'eb'; % Unused
                    case {'1.2.840.10008.1.2.4.50','1.2.840.10008.1.2.4.51','1.2.840.10008.1.2.4.70',...
                          '1.2.840.10008.1.2.4.80','1.2.840.10008.1.2.4.90','1.2.840.10008.1.2.4.91'}, % JPEG Explicit VR
                        flg = 'el';
                        %warning('spm:dicom',['Cant read JPEG Encoded file "' fopen(fp) '".']);
                    otherwise,
                        flg = 'el';
                        warning('spm:dicom','%s: Unknown Transfer Syntax UID (%s).',fopen(fp), dat);
                        return;
                end;
            otherwise,
                switch tag.vr,
                    case {'UN'},
                        % Unknown - read as char
                        dat = fread(fp,tag.length,'uint8')';
                    case {'AE', 'AS', 'CS', 'DA', 'DS', 'DT', 'IS', 'LO', 'LT',...
                            'PN', 'SH', 'ST', 'TM', 'UI', 'UT'},
                        % Character strings
                        dat = char(fread(fp,tag.length,'uint8')');

                        switch tag.vr,
                            case {'UI','ST'},
                                dat = deblank(dat);
                            case {'DS'},
                                try
                                    dat = textscan(dat,'%f','delimiter','\\')';
                                    dat = dat{1};
                                catch
                                    dat = textscan(dat,'%f','delimiter','/')';
                                    dat = dat{1};
                                end
                            case {'IS'},
                                dat = textscan(dat,'%d','delimiter','\\')';
                                dat = double(dat{1});
                            case {'DA'},
                                dat     = strrep(dat,'.',' ');
                                dat     = textscan(dat,'%4d%2d%2d');
                                [y,m,d] = deal(dat{:});
                                dat     = datenum(double(y),double(m),double(d));
                            case {'TM'},
                                if any(dat==':'),
                                    dat     = textscan(dat,'%d:%d:%f');
                                    [h,m,s] = deal(dat{:});
                                    h       = double(h);
                                    m       = double(m);
                                else
                                    dat     = textscan(dat,'%2d%2d%f');
                                    [h,m,s] = deal(dat{:});
                                    h       = double(h);
                                    m       = double(m);
                                end
                                if isempty(h), h = 0; end;
                                if isempty(m), m = 0; end;
                                if isempty(s), s = 0; end;
                                dat = s+60*(m+60*h);
                            case {'LO'},
                                dat = uscore_subst(dat);
                            otherwise,
                        end;
                    case {'OB'},
                        % dont know if this should be signed or unsigned
                        dat = fread(fp,tag.length,'uint8')';
                    case {'US', 'AT', 'OW'},
                        dat = fread(fp,tag.length/2,'uint16')';
                    case {'SS'},
                        dat = fread(fp,tag.length/2,'int16')';
                    case {'UL'},
                        dat = fread(fp,tag.length/4,'uint32')';
                    case {'SL'},
                        dat = fread(fp,tag.length/4,'int32')';
                    case {'FL'},
                        dat = fread(fp,tag.length/4,'float')';
                    case {'FD'},
                        dat = fread(fp,tag.length/8,'double')';
                    case {'SQ'},
                        [dat,len1] = read_sq(fp, flg,dict,tag.length);
                        tag.length = len1;
                    otherwise,
                        dat = '';
                        if tag.length
                            fseek(fp,tag.length,'cof');
                            warning('spm:dicom','%s: Unknown VR [%X%X] (offset=%d).', fopen(fp), tag.vr+0, ftell(fp));
                        end
                end;
                if ~isempty(tag.name),
                    ret.(tag.name) = dat;
                end;
        end;
    end;
    len = len + tag.le + tag.length;
end;
%_______________________________________________________________________

%_______________________________________________________________________
function [ret,len] = read_sq(fp, flg, dict,lim)
ret = {};
n   = 0;
len = 0;
while len<lim,
    tag.group   = fread(fp,1,'ushort');
    tag.element = fread(fp,1,'ushort');
    tag.length  = fread(fp,1,'uint');
    if isempty(tag.length), return; end; % End of file

    %fprintf('(%.4X,%.4X) %d\n', tag.group,tag.element,tag.length);
    %if tag.length==13, tag.length=10; end;

    len         = len + 8;
    if (tag.group == 65534) && (tag.element == 57344),   % FFFE/E000 Item
        [Item,len1] = read_dicom(fp, flg, dict, tag.length);
        len    = len + len1;
        if ~isempty(Item)
            n      = n + 1;
            ret{n} = Item;
        end
    elseif (tag.group == 65279) && (tag.element == 224), % FEFF/00E0 Item (Byte-swapped)
        % Byte-swapped
        [fname,perm,fmt] = fopen(fp);
        flg1 = flg;
        if flg(2)=='b',
            flg1(2) = 'l';
        else
            flg1(2) = 'b';
        end;
        [Item,len1] = read_dicom(fp, flg1, dict, tag.length);
        len    = len + len1;
        n      = n + 1;
        ret{n} = Item;
        pos    = ftell(fp);
        fclose(fp);
        fp     = fopen(fname,perm,fmt);
        fseek(fp,pos,'bof');
    elseif (tag.group == 65534) && (tag.element == 57565), % FFFE/E0DD SequenceDelimitationItem
        break;
    elseif (tag.group == 65279) && (tag.element == 56800), % FEFF/DDE0 SequenceDelimitationItem (Byte-swapped)
        % Byte-swapped
        break;
    else
        warning('spm:dicom','%s: Tag (%.4X,%.4X) is unexpected in sequence.', fopen(fp), tag.group, tag.element);
    end;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function tag = read_tag(fp,flg,dict)
tag.group   = fread(fp,1,'ushort');
tag.element = fread(fp,1,'ushort');
if isempty(tag.element), tag=[]; return; end;
if tag.group == 2, flg = 'el'; end;
%t          = dict.tags(tag.group+1,tag.element+1); % Sparse matrix representation
t           = find(dict.group==tag.group & dict.element==tag.element);
if t>0,
    tag.name = dict.values(t).name;
    tag.vr   = dict.values(t).vr{1};
else
    % Set tag.name = '' in order to restrict the fields to those
    % in the dictionary.  With a reduced dictionary, this could
    % speed things up considerably.
    % tag.name = '';
    if tag.element~=0,
        if rem(tag.group,2),
            tag.name = sprintf('Private_%.4x_%.4x',tag.group,tag.element);
        else
            tag.name = sprintf('Tag_%.4x_%.4x',tag.group,tag.element);
        end
        tag.vr   = 'UN';
    else
        tag.name = '';
        tag.vr   = 'UN';
    end
end;

if flg(2) == 'b',
    [fname,perm,fmt] = fopen(fp);
    if strcmp(fmt,'ieee-le') || strcmp(fmt,'ieee-le.l64'),
        pos = ftell(fp);
        fclose(fp);
        fp  = fopen(fname,perm,'ieee-be');
        fseek(fp,pos,'bof');
    end;
end;

if flg(1) =='e',
    tag.vr      = char(fread(fp,2,'uint8')');
    tag.le      = 6;
    switch tag.vr,
        case {'OB','OW','SQ','UN','UT'}
            if ~strcmp(tag.vr,'UN') || tag.group~=65534,
                unused = fread(fp,1,'ushort');
                tag.le = 12;
            else
                warning('spm:dicom','%s: Possible problem with %s tag (VR="%s").', fopen(fp), tag.name, tag.vr);
                tag.le = 10;
            end;
            tag.length = double(fread(fp,1,'uint'));
        case {'AE','AS','AT','CS','DA','DS','DT','FD','FL','IS','LO','LT','PN','SH','SL','SS','ST','TM','UI','UL','US'},
            tag.length = double(fread(fp,1,'ushort'));
            tag.le     = 8;
        case char([0 0])
            if (tag.group == 65534) && (tag.element == 57357)    % ItemDeliminationItem
                % at least on GE, ItemDeliminationItem does not have a
                % VR, but 4 bytes zeroes as length
                tag.vr     = 'UN';
                tag.le     = 8;
                tag.length = 0;
                unused     = fread(fp,1,'ushort'); % Should be zero
            elseif (tag.group == 65534) && (tag.element == 57565) % SequenceDelimitationItem
                tag.vr     = 'UN';
                tag.le     = 8;
                tag.length = 0;
                unused     = fread(fp,1,'ushort'); % Should be zero
            else
                warning('spm:dicom','%s: Don''t know how to handle VR of "\0\0" in %s.',fopen(fp),tag.name);
            end;
        otherwise,
            warning('spm:dicom','%s: Possible problem with %s tag (VR="%s")\n', fopen(fp), tag.name, tag.vr);
            unused = fread(fp,1,'ushort');
            tag.length = double(fread(fp,1,'uint'));
            tag.le     = 12;
    end;
else
    tag.le     =  8;
    tag.length = double(fread(fp,1,'uint'));
end;

if isempty(tag.vr) || isempty(tag.length),
    tag = [];
    return;
end;


if rem(tag.length,2),
    if tag.length==4294967295, % FFFFFFFF
        return;
    else
        warning('spm:dicom','%s: Odd numbered Value Length in %s tag (%X).', fopen(fp), tag.name, tag.length);
        tag = [];
    end;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function dict = readdict(P)
if nargin<1, P = 'spm_dicom_dict.mat'; end;
try
    dict = load(P);
catch problem
    fprintf('\nUnable to load the file "%s".\n', P);
    rethrow(problem);
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function t = decode_csa(fp,lim)
% Decode shadow information (0029,1010) and (0029,1020)
[fname,perm,fmt] = fopen(fp);
pos = ftell(fp);
if strcmp(fmt,'ieee-be') || strcmp(fmt,'ieee-be.l64'),
    fclose(fp);
    fp  = fopen(fname,perm,'ieee-le');
    fseek(fp,pos,'bof');
end;

c   = fread(fp,4,'uint8');
fseek(fp,pos,'bof');

if all(c'==[83 86 49 48]), % "SV10"
    t = decode_csa2(fp,lim);
else
    t = decode_csa1(fp,lim);
end;

if strcmp(fmt,'ieee-be') || strcmp(fmt,'ieee-be.l64'),
    fclose(fp);
    fp  = fopen(fname,perm,fmt);
end;
fseek(fp,pos+lim,'bof');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function t = decode_csa1(fp,lim)
n   = fread(fp,1,'uint32');
if isempty(n) || n>1024 || n <= 0,
    fseek(fp,lim-4,'cof');
    t = struct('name','JUNK: Don''t know how to read this damned file format');
    return;
end;
unused = fread(fp,1,'uint32')'; % Unused "M" or 77 for some reason
tot = 2*4;
for i=1:n,
    t(i).name    = fread(fp,64,'uint8')';
    msk          = find(~t(i).name)-1;
    if ~isempty(msk),
        t(i).name    = char(t(i).name(1:msk(1)));
    else
        t(i).name    = char(t(i).name);
    end;
    t(i).vm      = fread(fp,1,'int32')';
    t(i).vr      = fread(fp,4,'uint8')';
    t(i).vr      = char(t(i).vr(1:3));
    t(i).syngodt = fread(fp,1,'int32')';
    t(i).nitems  = fread(fp,1,'int32')';
    t(i).xx      = fread(fp,1,'int32')'; % 77 or 205
    tot          = tot + 64+4+4+4+4+4;
    for j=1:t(i).nitems
        % This bit is just wierd
        t(i).item(j).xx  = fread(fp,4,'int32')'; % [x x 77 x]
        len              = t(i).item(j).xx(1)-t(1).nitems;
        if len<0 || len+tot+4*4>lim,
            t(i).item(j).val = '';
            tot              = tot + 4*4;
            break;
        end;
        t(i).item(j).val = char(fread(fp,len,'uint8')');
        fread(fp,4-rem(len,4),'uint8');
        tot              = tot + 4*4+len+(4-rem(len,4));
    end;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function t = decode_csa2(fp,lim)
unused1 = fread(fp,4,'uint8'); % Unused
unused2 = fread(fp,4,'uint8'); % Unused
n    = fread(fp,1,'uint32');
opos = ftell(fp);
if isempty(n) || n>1024 || n < 0,
    fseek(fp,lim-4,'cof');
    t = struct('name','Don''t know how to read this damned file format');
    return;
end;
unused = fread(fp,1,'uint32')'; % Unused "M" or 77 for some reason
pos    = 16;
for i=1:n,
    t(i).name    = fread(fp,64,'uint8')';
    pos          = pos + 64;
    msk          = find(~t(i).name)-1;
    if ~isempty(msk),
        t(i).name    = char(t(i).name(1:msk(1)));
    else
        t(i).name    = char(t(i).name);
    end;
    t(i).vm      = fread(fp,1,'int32')';
    t(i).vr      = fread(fp,4,'uint8')';
    t(i).vr      = char(t(i).vr(1:3));
    t(i).syngodt = fread(fp,1,'int32')';
    t(i).nitems  = fread(fp,1,'int32')';
    t(i).xx      = fread(fp,1,'int32')'; % 77 or 205
    pos          = pos + 20;
    for j=1:t(i).nitems
        t(i).item(j).xx  = fread(fp,4,'int32')'; % [x x 77 x]
        pos              = pos + 16;
        len              = t(i).item(j).xx(2);
        if len>lim-pos,
            len = lim-pos;
            t(i).item(j).val = char(fread(fp,len,'uint8')');
            fread(fp,rem(4-rem(len,4),4),'uint8');
            warning('spm:dicom','%s: Problem reading Siemens CSA field.', fopen(fp));
            return;
        end
        t(i).item(j).val = char(fread(fp,len,'uint8')');
        fread(fp,rem(4-rem(len,4),4),'uint8');
    end;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function str_out = uscore_subst(str_in)
str_out = str_in;
pos = strfind(str_in,'+AF8-');
if ~isempty(pos),
    str_out(pos) = '_';
    str_out(repmat(pos,4,1)+repmat((1:4)',1,numel(pos))) = [];
end
return;
%_______________________________________________________________________

