spm fmri;
spm_jobman;
bb=findobj(0,'tag','batch_box');
c0spm=get(bb,'userdata');
c0=cfg_struct2cfg(c0spm);
cfg_util('addapp',c0)
% stop here if you don't want to create this example job
mod_cfg_id=cfg_util('listcfg',[],cfg_findspec({{'tag','imcalc'}}))
mod_job_id{1} = cfg_util('addtojob', mod_cfg_id{1})
mod_job_id{2} = cfg_util('addtojob', mod_cfg_id{1})
mod_job_id{3} = cfg_util('addtojob', mod_cfg_id{1})
cfg_util('delfromjob',mod_job_id{1})
inp = cfg_util('listmod', mod_job_id{2}, [], cfg_findspec({{'tag','input'}}));
cfg_util('setval', mod_job_id{2}, inp{1}, {'/export/T1/c1rtAA240407_HMRF.img',...
    '/export/T1/c1rtAA290604_HMRF.img',...
    '/export/T1/c1rtAB101104_HMRF.img',...
    '/export/T1/c1rtAB130704_HMRF.img',...
    '/export/T1/c1rtAB160904_HMRF.img',...
    '/export/T1/c1rtAE170505_HMRF.img'});
[u j]=cfg_util('harvest',mod_job_id{2})
[mod_job_idlist str sts dep sout] = cfg_util('showjob');
cfg_util('setval', mod_job_id{3},inp{1},sout{1})
cfg_util('harvest',mod_job_id{3});
[mod_job_idlist str sts dep sout] = cfg_util('showjob');
cfg_ui