function [this, ok] = check(this, option)
% Method that performs integrity checks of the meeg object
% and its readiness for particular purposes.
% FORMAT  this = check(this, option)
% IN
% option - 'basic' (default) - just check the essential fields
%          '3d' - check if suitable for source reconstruction
%          'dcm'    - check if suitable for DCM
%
% OUT
% ok - 1 - OK, 0- failed
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: check.m 6542 2015-09-09 11:48:34Z karl $

if nargin == 1
    option = 'basic';
end

ok = 1;

this = meeg(struct(this));

if isequal(option, 'basic')
    return;
end
        
if ~isequal(transformtype(this), 'time')
    ok = 0;
    disp('Source reconstruction and DCM only work for time domain data');
    return;
end

if this.montage.Mind~=0
    disp('Virtual montage is applied. Make sure this is what you want.')    
end

eegind    = indchantype(this, 'EEG');
megind    = indchantype(this, {'MEG'});
planarind = indchantype(this, {'MEGPLANAR'});
lfpind    = indchantype(this, 'LFP');

if  ~isempty([eegind(:); megind(:); planarind(:)])              
    if ~isempty(eegind)
        if ~isfield(this.sensors, 'eeg') || isempty(this.sensors.eeg)
            ok = 0;
            disp('EEG sensor locations are not specified');
        else
            if ~isempty(setdiff(chanlabels(this, eegind), this.sensors.eeg.label))
                ok = 0;
                disp('Not all EEG channel locations are specified');               
            end
        end
        
        if ~isempty(strmatch('unknown', units(this, eegind)))
            this = units(this, eegind, 'uV');
            
            warning_flexible('EEG channel units are missing. Assuming uV, source scaling might be wrong');
        end
    end
    
    if ~isempty([megind(:); planarind(:)])
        if ~isfield(this.sensors, 'meg') || isempty(this.sensors.meg)
            ok = 0;
            disp('MEG sensor locations are not specified');
        else
            if ~isempty(setdiff({this.channels(megind).label}, this.sensors.meg.label))                
                disp('Not all MEG sensor locations are specified');         
            end
        end
        
        if ~isempty(strmatch('unknown', units(this, megind)))
            this = units(this,  megind, 'fT');
            
             warning_flexible('MEG channel units are missing. Assuming fT, source scaling might be wrong');
        end
        
        if ~isempty(strmatch('unknown', units(this, planarind)))
            this = units(this,  planarind, 'fT/mm');
            
            warning_flexible('MEGPLANAR channel units are missing. Assuming fT/mm, source scaling might be wrong');
        end
    end
    
    if isempty(this.fiducials)
        ok = 0;
        disp('No fiducials are defined');
    end
    
    if ~isfield(this.fiducials, 'pnt') || isempty(this.fiducials.pnt)
        if ~isempty(eegind)
            % Copy EEG sensors to fiducials.
            this.fiducials.pnt = this.sensors.eeg.elecpos;
        else
            this.fiducials.pnt = sparse(0, 3);
        end
    end
    
    if ~isfield(this.fiducials, 'fid') || ...
            ~all(isfield(this.fiducials.fid, {'pnt', 'label'})) ||...
            (length(this.fiducials.fid.label) ~= size(this.fiducials.fid.pnt, 1)) || ...
            length(this.fiducials.fid.label) < 3
        ok = 0;
        disp('At least 3 fiducials with labels are required');        
    end
    
    nzlbl = {'fidnz', 'nz', 'nas', 'nasion', 'spmnas'};
    lelbl = {'fidle', 'fidt9', 'lpa', 'lear', 'earl', 'le', 'l', 't9', 'spmlpa'};
    relbl = {'fidre', 'fidt10', 'rpa', 'rear', 'earr', 're', 'r', 't10', 'spmrpa'};
    
    [sel1, nzind] = match_str(nzlbl, lower(this.fiducials.fid.label));
    if isempty(nzind)
        disp('Could not find the nasion fiducial');
    else
        nzind = nzind(1);
    end
    
    [sel1, leind] = match_str(lelbl, lower(this.fiducials.fid.label));
    if isempty(leind)
        disp('Could not find the left fiducial');
    else
        leind = leind(1);
    end
    
    [sel1, reind] = match_str(relbl, lower(this.fiducials.fid.label));
    if isempty(reind)
        disp('Could not find the right fiducial');
    else
        reind = reind(1);
    end
    
    restind = setdiff(1:length(this.fiducials.fid.label), [nzind(:)', leind(:)', reind(:)']);
    
    this.fiducials.fid.label = this.fiducials.fid.label([nzind(:)', leind(:)', reind(:)', restind(:)']);
    this.fiducials.fid.pnt = this.fiducials.fid.pnt([nzind(:)', leind(:)', reind(:)', restind(:)'], :);
    
end

if isequal(option, '3d')
    if ~ismember(modality(this), {'EEG', 'MEG', 'Multimodal'})
        ok = 0;
        disp('Unsupported modality for 3D source reconstruction');
    end    
end

if isequal(option, 'dcm')
    if strcmp(option, 'dcm')
        if ~ismember(modality(this, 0), {'EEG', 'MEG', 'MEGPLANAR', 'Multimodal', 'LFP','ILAM'})
            ok = 0;
            disp('Unsupported modality for DCM');            
        end
    end
    
    if ~isempty(lfpind) && ~isempty([eegind, megind])
        ok = 0;
        disp('DCM does not allow mixing scalp and LFP channels in the same dataset');
    end           
end