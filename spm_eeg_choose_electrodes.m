function [channels,data,d1]=choose_electrodes;


D=spm_eeg_ldata;
[Cel, Cind, x, y] = spm_eeg_locate_channels(D, 128, 1);
d1=[];
for n=1:size(Cel,1)
	d1=[d1,data(Cel(n,1),Cel(n,2))];
end
figure
imagesc(data(:,:)')
hold on
for n=1:275
	text(Cel(n,1),Cel(n,2),[D.channels.name(n)])
	hold on
end
but =1;
channels=zeros(1,275);
while but == 1
	[x1,y1,but]=ginput(1);
	tmp=((Cel(:,1)-x1).^2+(Cel(:,2)-y1).^2).^0.5;
	[junk,elec]=sort(tmp);
	if but ==1
		if channels(elec(1))==0
			channels(elec(1))=1;
		else
			channels(elec(1))=0;
		end
	end
	clf
	imagesc(data(:,:)')
	hold on
	
	plot(Cel(channels==0,1),Cel(channels==0,2),'wo')
end

