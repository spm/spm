% prompt for OK and activate correct figure
%--------------------------------------------------------------------------

uiwait(warndlg('Proceed with demonstration?','MFM demo'));
spm_figure('GetWin','MFM');
clf