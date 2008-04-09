function grad = yokogawa2grad(hdr);

% YOKOGAWA2GRAD converts the position and weights of all coils that
% compromise a gradiometer system into a structure that can be used
% by FieldTrip.

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: yokogawa2grad.m,v $
% Revision 1.1  2006/08/31 13:32:11  roboos
% moved from fieldtrip to fileio module
%
% Revision 1.1  2005/09/06 08:54:30  roboos
% new implementations for the Yokogawa 160 channel MEG system
%

% hdr = read_yokogawa_header(filename);
hdr = hdr.orig; % use the original Yokogawa header, not the FieldTrip header

handles    = definehandles;
isgrad     = (hdr.channel_info(:,2)==handles.AxialGradioMeter);
grad.pnt   = hdr.channel_info(find(isgrad),3:5);
grad.ori   = []; % FIXME
grad.label = []; % FIXME

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this defines some usefull constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = definehandles;
handles.output = [];
handles.sqd_load_flag = false;
handles.mri_load_flag = false;
handles.NullChannel         = 0;
handles.MagnetoMeter        = 1;
handles.AxialGradioMeter    = 2;
handles.PlannerGradioMeter  = 3;
handles.RefferenceChannelMark = hex2dec('0100');
handles.RefferenceMagnetoMeter       = bitor( handles.RefferenceChannelMark, handles.MagnetoMeter );
handles.RefferenceAxialGradioMeter   = bitor( handles.RefferenceChannelMark, handles.AxialGradioMeter );
handles.RefferencePlannerGradioMeter = bitor( handles.RefferenceChannelMark, handles.PlannerGradioMeter );
handles.TriggerChannel      = -1;
handles.EegChannel          = -2;
handles.EcgChannel          = -3;
handles.EtcChannel          = -4;
handles.NonMegChannelNameLength = 32;
handles.DefaultMagnetometerSize       = (4.0/1000.0);		% ˆê•Ó4.0mm‚Ì?³•ûŒ`
handles.DefaultAxialGradioMeterSize   = (15.5/1000.0);		% ’¼Œa15.5mm‚Ì‰~ŠÂ
handles.DefaultPlannerGradioMeterSize = (12.0/1000.0);		% ˆê•Ó12.0mm‚Ì?³•ûŒ`
handles.AcqTypeContinuousRaw = 1;
handles.AcqTypeEvokedAve     = 2;
handles.AcqTypeEvokedRaw     = 3;
handles.sqd = [];
handles.sqd.selected_start  = [];
handles.sqd.selected_end    = [];
handles.sqd.axialgradiometer_ch_no      = [];
handles.sqd.axialgradiometer_ch_info    = [];
handles.sqd.axialgradiometer_data       = [];
handles.sqd.plannergradiometer_ch_no    = [];
handles.sqd.plannergradiometer_ch_info  = [];
handles.sqd.plannergradiometer_data     = [];
handles.sqd.nullchannel_ch_no   = [];
handles.sqd.nullchannel_data    = [];
handles.sqd.selected_time       = [];
handles.sqd.sample_rate         = [];
handles.sqd.sample_count        = [];
handles.sqd.pretrigger_length   = [];
handles.sqd.matching_info   = [];
handles.sqd.source_info     = [];
handles.sqd.mri_info        = [];
handles.mri                 = [];
