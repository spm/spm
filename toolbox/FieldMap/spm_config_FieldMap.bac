function vdm = spm_config_fieldmap
% Configuration file for FieldMap jobs
%_______________________________________________________________________
% Copyright (C) 2006 Wellcome Department of Imaging Neuroscience

% Chloe Hutton
% $Id$
%_______________________________________________________________________
entry = inline(['struct(''type'',''entry'',''name'',name,'...
        '''tag'',tag,''strtype'',strtype,''num'',num)'],...
        'name','tag','strtype','num');

files = inline(['struct(''type'',''files'',''name'',name,'...
        '''tag'',tag,''filter'',fltr,''num'',num)'],...
        'name','tag','fltr','num');

mnu = inline(['struct(''type'',''menu'',''name'',name,'...
        '''tag'',tag,''labels'',{labels},''values'',{values})'],...
        'name','tag','labels','values');

branch = inline(['struct(''type'',''branch'',''name'',name,'...
        '''tag'',tag,''val'',{val})'],...
        'name','tag','val');

repeat = inline(['struct(''type'',''repeat'',''name'',name,''tag'',tag,'...
         '''values'',{values})'],'name','tag','values');
     
choice = inline(['struct(''type'',''choice'',''name'',name,''tag'',tag,'...
         '''values'',{values})'],'name','tag','values');
%______________________________________________________________________

addpath(fullfile(spm('dir'),'toolbox','FieldMap'));

phase = files('Phase','p','image',1);
mag   = files('Magnitude','m','image',1);
te    = entry('Echo Time','te','e',[1 1]);
ste   = branch('Short TE','ste',{te,phase,mag});
lte   = branch('Long TE','lte',{te,phase,mag});
pm    = branch('PM','pm',{ste,lte});
rl    = files('Real','r','image',1);
im    = files('Imaginary','i','image',1);
ste   = branch('Short TE','ste',{te,rl,im});
lte   = branch('Long TE','lte',{te,rl,im});
ri    = branch('RI','ri',{ste,lte});
precalc = files('Pre-Calculated','precalc','image',1);
fm    = choice('Field Map','fm',{ri,pm,precalc});
fm.prog = @run_fmap;

ebf = mnu('EPI-based field map','ebf',{'Yes','No'},{1,0}); ebf.val = {0};
ped = mnu('Phase Encode Direction','ped',{'X','Y'},{1,2}); ped.val = {2};
pol = mnu('Polarity of phase-encode blips','pol',{'+ve','-ve'},{1,-1}); pol.val = {1};
mod = mnu('Jacobian Modulation','mod',{'Yes','No'},{1,0}); mod.val = {0};

vdm = branch('Vox Disp Map','vdm',{fm,ebf,ped,pol,mod});
vdm.help = {'This toolbox does nothing yet','','There is no help either.'};

function run_fmap(job)
disp(job)
