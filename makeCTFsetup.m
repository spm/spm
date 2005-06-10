%%%% function to make  CTF setup file

pre_data=ctf_read;
sens=pre_data.sensor.index.all_sens;
xd=pre_data.sensor.location(1,sens);
yd=pre_data.sensor.location(2,sens);
zd=pre_data.sensor.location(3,sens);
[th,r]=cart2pol(xd,yd);
[x,xi]=min(r);
xdata=pre_data.sensor.location(1,sens)-xd(xi);
ydata=pre_data.sensor.location(2,sens)-yd(xi);
zdata=pre_data.sensor.location(3,sens)-zd(xi);
[th,r]=cart2pol(xdata,ydata);
r=(xdata.^2+ydata.^2+zdata.^2).^0.5;
[xdata,ydata]=pol2cart(th,r);
plot(xdata,ydata,'.')


Cnames=pre_data.sensor.label(sens);;
Nchannels=length(sens);
Cpos(1,:)=(xdata+abs(min(xdata).*1.2))./((max(xdata)-min(xdata))*1.2);
Cpos(2,:)=(ydata+abs(min(ydata).*1.2))./((max(ydata)-min(ydata))*1.2);
figure
plot(Cpos(1,:),Cpos(2,:),'.')

Rxy=1.5;

save('d:/spm5/EEGtemplates/CTF275_setup','Cnames','Cpos','Rxy','Nchannels');