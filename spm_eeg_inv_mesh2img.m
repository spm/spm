function Vo = spm_eeg_inv_mesh2img(vert,Vi,fname)

Vo = Vi;
Vo.fname = fname;

dat = zeros(Vo.dim);

n = length(vert);

for i = 1:n
    dat(round(vert(i,1)),round(vert(i,2)),round(vert(i,3))) = 1;
end

Vo = spm_write_vol(Vo,dat);

